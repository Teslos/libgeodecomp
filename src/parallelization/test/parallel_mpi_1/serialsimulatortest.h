#include <cxxtest/TestSuite.h>
#include <sstream>
#include <libgeodecomp/io/writer.h>
#include <libgeodecomp/io/memorywriter.h>
#include <libgeodecomp/io/mockinitializer.h>
#include <libgeodecomp/io/mockwriter.h>
#include <libgeodecomp/io/mocksteerer.h>
#include <libgeodecomp/io/testinitializer.h>
#include <libgeodecomp/io/testwriter.h>
#include <libgeodecomp/misc/stringops.h>
#include <libgeodecomp/misc/testcell.h>
#include <libgeodecomp/misc/testhelper.h>
#include <libgeodecomp/parallelization/serialsimulator.h>

using namespace LibGeoDecomp;

namespace LibGeoDecomp {

class SerialSimulatorTest : public CxxTest::TestSuite
{
public:
    typedef MockSteerer<TestCell<2> > SteererType;

    void setUp()
    {
        dim = Coord<2>(17, 12);
        maxSteps = 21;
        startStep = 13;

        init.reset(createInitializer());
        simulator.reset(new SerialSimulator<TestCell<2> >(createInitializer()));
    }

    void tearDown()
    {
        simulator.reset();
        init.reset();
    }

    void testInitialization()
    {
        TS_ASSERT_EQUALS(simulator->getGrid()->getDimensions().x(),
                         (unsigned)17);
        TS_ASSERT_EQUALS(simulator->getGrid()->getDimensions().y(),
                         (unsigned)12);
        TS_ASSERT_TEST_GRID(Grid<TestCell<2> >, *simulator->getGrid(), startStep * TestCell<2>::nanoSteps());
    }

    void testStep()
    {
        TS_ASSERT_EQUALS(startStep, (int)simulator->getStep());

        simulator->step();
        const Grid<TestCell<2> > *grid = simulator->getGrid();
        TS_ASSERT_TEST_GRID(Grid<TestCell<2> >, *grid,
                            (startStep + 1) * TestCell<2>::nanoSteps());
        TS_ASSERT_EQUALS(startStep + 1, (int)simulator->getStep());
    }

    void testRun()
    {
        simulator->run();
        TS_ASSERT_EQUALS(init->maxSteps(), simulator->getStep());
        TS_ASSERT_TEST_GRID(
            Grid<TestCell<2> >,
            *simulator->getGrid(),
            init->maxSteps() * TestCell<2>::nanoSteps());
    }

    // fixme: use unified test suite like this test for all simulators?
    void testWriterInvocation()
    {
        unsigned period = 4;
        SuperVector<unsigned> expectedSteps;
        SuperVector<WriterEvent> expectedEvents;
        expectedSteps << 13
                      << 16
                      << 20
                      << 21;
        expectedEvents << WRITER_INITIALIZED
                       << WRITER_STEP_FINISHED
                       << WRITER_STEP_FINISHED
                       << WRITER_ALL_DONE;

        simulator->addWriter(new TestWriter(period, expectedSteps, expectedEvents));
        simulator->run();
    }

    //fixme: move those tests converning just the abstract base class
    //Simulator to unitsimulator.h
    void testDeleteInitializer()
    {
        MockInitializer::events = "";
        {
            SerialSimulator<TestCell<2> > foo(new MockInitializer);
        }
        TS_ASSERT_EQUALS(
            MockInitializer::events,
            "created, configString: ''\ndeleted\n");
    }

    void testRegisterWriter()
    {
        MockWriter *w = new MockWriter();
        simulator->addWriter(w);
        SerialSimulator<TestCell<2> >::WriterVector writers = simulator->writers;
        TS_ASSERT_EQUALS(1, writers.size());
        TS_ASSERT_EQUALS(w, writers[0].get());
    }

    void testSerialSimulatorShouldCallBackWriter()
    {
        MockWriter *w = new MockWriter(3);
        simulator->addWriter(w);
        simulator->run();

        std::string expectedEvents = "initialized()\n";
        for (unsigned i = startStep + 2; i <= init->maxSteps(); i += 3) {
            expectedEvents +=
                "stepFinished(step=" + StringConv::itoa(i) + ")\n";
        }
        expectedEvents += "allDone()\n";

        TS_ASSERT_EQUALS(expectedEvents, w->events());
    }

    void testRunMustResetGridPriorToSimulation()
    {
        MockWriter *eventWriter1 = new MockWriter();
        MemoryWriter<TestCell<2> > *gridWriter1 =
            new MemoryWriter<TestCell<2> >();
        simulator->addWriter(eventWriter1);
        simulator->addWriter(gridWriter1);

        simulator->run();
        std::string events1 = eventWriter1->events();
        std::vector<Grid<TestCell<2> > > grids1 = gridWriter1->getGrids();

        MockWriter *eventWriter2 = new MockWriter();
        MemoryWriter<TestCell<2> > *gridWriter2 =
            new MemoryWriter<TestCell<2> >();
        simulator->addWriter(eventWriter2);
        simulator->addWriter(gridWriter2);
        simulator->run();
        std::string events2 = eventWriter2->events();
        std::vector<Grid<TestCell<2> > > grids2 = gridWriter2->getGrids();

        TS_ASSERT_EQUALS(events1, events2);
        TS_ASSERT_EQUALS(grids1, grids2);
    }

    typedef Grid<TestCell<3>, TestCell<3>::Topology> Grid3D;

    void test3D()
    {
        SerialSimulator<TestCell<3> > sim(new TestInitializer<TestCell<3> >());
        TS_ASSERT_TEST_GRID(Grid3D, *sim.getGrid(), 0);

        sim.step();
        TS_ASSERT_TEST_GRID(Grid3D, *sim.getGrid(),
                            TestCell<3>::nanoSteps());

        sim.nanoStep(0);
        TS_ASSERT_TEST_GRID(Grid3D, *sim.getGrid(),
                            TestCell<3>::nanoSteps() + 1);

        sim.run();
        TS_ASSERT_TEST_GRID(Grid3D, *sim.getGrid(),
                            21 * TestCell<3>::nanoSteps());
    }

    void testSteererCallback()
    {
        std::stringstream events;
        simulator->addSteerer(new SteererType(5, &events));
        std::stringstream expected;
        expected << "created, period = 5\n";
        TS_ASSERT_EQUALS(events.str(), expected.str());

        simulator->run();
        int i = startStep;
        if (i % 5) {
            i += 5 - (i % 5);
        }
        for (; i < maxSteps; i += 5) {
            expected << "nextStep(" << i << ")\n";
        }
        expected << "deleted\n";

        simulator.reset();
        TS_ASSERT_EQUALS(events.str(), expected.str());
    }

private:
    boost::shared_ptr<SerialSimulator<TestCell<2> > > simulator;
    boost::shared_ptr<Initializer<TestCell<2> > > init;
    unsigned maxSteps;
    unsigned startStep;
    Coord<2> dim;

    Initializer<TestCell<2> > *createInitializer()
    {
        return new TestInitializer<TestCell<2> >(dim, maxSteps, startStep);
    }
};

}
