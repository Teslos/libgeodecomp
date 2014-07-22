// vim: noai:ts=4:sw=4:expandtab
#include <libgeodecomp/misc/simplexoptimizer.h>
#include <libgeodecomp/misc/test/unit/patternoptimizertest.h>
#include <libgeodecomp/io/logger.h>

using namespace LibGeoDecomp;

namespace LibGeoDecomp {

class SimplexOptimizerTest :  public CxxTest::TestSuite
{
public:
    void testHimmelblau()
    {
        SimulationParameters params;
        params.addParameter("x", -500, 500);
        params.addParameter("y", -500, 500);
        params["x"].setValue(800);
        params["y"].setValue(400);

        PatternOptimizerTest::HimmelblauFunction eval;
        SimplexOptimizer optimizer(params);
        params = optimizer(20, eval);

        // for now it's OK if the optimum is not realized, but we
        // should get really close to it:
        TS_ASSERT(optimizer.fitness > (0.999 * eval.getGlobalMax()));

        LOG(Logger::INFO, "Test Himmelblau with simplexOptimizer" << std::endl
            << "x: " << (int) params["x"] << std::endl
            << "y: " << (int) params["y"] << std::endl
            << "Calls: " << eval.getCalls() <<  std::endl);
    }

    void test5D()
    {
        SimulationParameters params5;
        params5.addParameter("v", -60, 60);
        params5.addParameter("w", -20, 40);
        params5.addParameter("x", -10, 10);
        params5.addParameter("y", 0, 2);
        params5.addParameter("z", -10, 40);

        SimplexOptimizer optimizer5(params5);
        PatternOptimizerTest::FiveDimFunction eval5;
        params5 = optimizer5(20 ,eval5);

        // for now it's OK if the optimum is not realized, but we
        // should get really close to it:
        TS_ASSERT(optimizer5.fitness > (0.997 * eval5.getGlobalMax()));

        LOG(Logger::INFO,  "Test 5, five dimensions: "
                << std::endl<< "fitness: " << optimizer5.fitness
                <<  " calls: " << eval5.getCalls() << std::endl
                << "v:"   << (int)params5["v"]
                << " w: " << (int)params5["w"]
                << " x: " << (int)params5["x"]
                << " y: " << (int)params5["y"]
                << " z: " << (int)params5["z"]
                << std::endl);
    }
};
} // namespace LibGeoDecomp
