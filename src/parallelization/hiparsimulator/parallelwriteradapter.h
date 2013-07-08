#ifndef LIBGEODECOMP_PARALLELIZATION_HIPARSIMULATOR_PARALLELWRITERADAPTER_H
#define LIBGEODECOMP_PARALLELIZATION_HIPARSIMULATOR_PARALLELWRITERADAPTER_H

#include <libgeodecomp/parallelization/hiparsimulator/patchaccepter.h>

namespace LibGeoDecomp {
namespace HiParSimulator {

template<class CELL_TYPE, class PARTITION>
class HiParSimulator;

/**
 * ParallelWriterAdapter translates the interface of a ParallelWriter
 * to a PatchAccepter, so that we can treat IO similarly to sending
 * ghost zones.
 */
template<typename GRID_TYPE, typename CELL_TYPE, typename SIMULATOR>
class ParallelWriterAdapter : public PatchAccepter<GRID_TYPE>
{
public:

    using PatchAccepter<GRID_TYPE>::checkNanoStepPut;
    using PatchAccepter<GRID_TYPE>::pushRequest;
    using PatchAccepter<GRID_TYPE>::requestedNanoSteps;

    ParallelWriterAdapter(
        SIMULATOR * sim,
        boost::shared_ptr<ParallelWriter<CELL_TYPE> > writer,
        const std::size_t firstStep,
        const std::size_t lastStep,
        Coord<CELL_TYPE::Topology::DIM> globalGridDimensions,
        bool lastCall) :
        sim(sim),
        writer(writer),
        firstNanoStep(firstStep * CELL_TYPE::nanoSteps()),
        lastNanoStep(lastStep   * CELL_TYPE::nanoSteps()),
        stride(writer->getPeriod() * CELL_TYPE::nanoSteps()),
        lastCall(lastCall),
        globalGridDimensions(globalGridDimensions)
    {
        pushRequest(firstNanoStep);
        pushRequest(lastNanoStep);
    }

    virtual void setRegion(const Region<GRID_TYPE::DIM>& region)
    {
        writer->setRegion(region);
    }

    virtual void put(
        const GRID_TYPE& grid,
        const Region<GRID_TYPE::DIM>& validRegion,
        const std::size_t nanoStep)
    {
        if (!checkNanoStepPut(nanoStep)) {
            return;
        }

        WriterEvent event = WRITER_STEP_FINISHED;
        if (nanoStep == firstNanoStep) {
            event = WRITER_INITIALIZED;
        }
        if (nanoStep == lastNanoStep) {
            event = WRITER_ALL_DONE;
        }

        writer->stepFinished(
            grid,
            validRegion,
            globalGridDimensions,
            nanoStep / CELL_TYPE::nanoSteps(),
            event,
            lastCall);
        requestedNanoSteps.erase_min();
        std::size_t nextNanoStep = nanoStep + stride;
        pushRequest(nextNanoStep);
    }

private:
    SIMULATOR * sim;
    boost::shared_ptr<ParallelWriter<CELL_TYPE> > writer;
    std::size_t firstNanoStep;
    std::size_t lastNanoStep;
    std::size_t stride;
    bool lastCall;
    Coord<CELL_TYPE::Topology::DIM> globalGridDimensions;
};

}
}

#endif
