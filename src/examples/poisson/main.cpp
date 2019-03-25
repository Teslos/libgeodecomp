//
// Created by ivt on 03.10.2017.
//
#include<iostream>
#include<libgeodecomp.h>
#include<libgeodecomp/communication/typemaps.h>
#include<libgeodecomp/io/bovwriter.h>
#include<libgeodecomp/io/tracingwriter.h>
#include<libgeodecomp/io/simpleinitializer.h>
#include<libgeodecomp/parallelization/hiparsimulator.h>
#include<libgeodecomp/loadbalancer/noopbalancer.h>
//#include<libgeodecomp/parallelization/stripingsimulator.h>

#ifdef PCUDA
#include <libgeodecomp/parallelization/cudasimulator.h>
#endif
//#define PMPI 1
#define PSERIAL 1
using namespace std;
using namespace LibGeoDecomp;

const int nx=101;
const int ny=101;

int outputFrequency = 1000;
int numSteps = 500;
const int nit = 50;

const double dx = 0.05;
const double mu = 0.1;
const double dt = 0.001;
const double rho = 1.;
const double nu = .1;

// model to solve
#include "poisson.hpp"

//MPI_Datatype Cell::MPIDataType;
class CellInitializer:public SimpleInitializer<Cell>
{
public:
    CellInitializer(Coord<2>dim, int numSteps): SimpleInitializer<Cell>(dim, numSteps)
    {}
    virtual void grid(GridBase<Cell,2> *ret)
    {
        CoordBox<2> box = ret->boundingBox();
        Coord<2> size = gridDimensions();
        /**
         * Set starting values for cells
         *
         */
	for (int y = 0; y < size.y(); ++y) {
            for (int x = 0; x < size.x(); ++x) {
                // get coordinates of cells
                Coord<2> c(x,y);
		Cell::State s = Cell::LIQUID;
                if (box.inBounds(c)) {
			ret->set(c, Cell(0.,0.,0.,s));
		}

            }
        }

        // set the boundary conditions according to cavity flow problem
        for(int x = 0; x < size.x(); x++) {
            Coord<2> cS(x, 0);
            Coord<2> cN(x, size.x()-1);
            if(box.inBounds(cS)){
		Cell::State s = Cell::SOUTH_NOSLIP;    
                ret->set( cS, Cell(0.,0.,0.,s));
            }

            if(box.inBounds(cN)){
		Cell::State s = Cell::NORTH_BOUND;    
                ret->set( cN, Cell(1.,0.,0.,s));
            }
        }

        for(int y = 0; y < size.y(); y++) {
            Coord<2> cE(0, y);
            Coord<2> cW(size.y()-1, y);
            if(box.inBounds(cE)) {
		Cell::State s = Cell::EAST_NOSLIP;
                ret->set( cE, Cell(0.,0.,0.,s));
            }
            if(box.inBounds(cW)) {
		Cell::State s = Cell::WEST_NOSLIP;    
                ret->set( cW, Cell(0.,0.,0.,s));
            }
        }
    }
};

/**
 * Sets and runs simulations
 */
void runSimulations()
{
#ifdef PMPI
    MPI_Aint displacements[] = { 0 };
    MPI_Datatype memberTypes[] = { MPI_CHAR };
    int lengths[] = { sizeof(Cell) };
    MPI_Type_create_struct(1, lengths, displacements, memberTypes, &Cell::MPIDataType);
    MPI_Type_commit(&Cell::MPIDataType);
#endif
    CellInitializer *init = new CellInitializer(Coord<2>(nx,ny),numSteps);

    int outputFrequency = 100;
    
#ifdef PMPI
    HiParSimulator<Cell, RecursiveBisectionPartition<2> > sim(
            init,
            MPILayer().rank() ? 0 : new NoOpBalancer(),
            numSteps,
            1);
/*     
    StripingSimulator<Cell> sim(
            init,
            MPILayer().rank() ? 0 : new NoOpBalancer());

*/
    /**
     * we can add different simulators for now we have only MPI simulator
     */
    if (MPILayer().rank() == 0)
        sim.addWriter(new TracingWriter<Cell>(outputFrequency, init->maxSteps()));
#endif
#ifdef PSERIAL
    SerialSimulator<Cell> sim(init);

    sim.addWriter(new TracingWriter<Cell>(outputFrequency, init->maxSteps()));
    sim.addWriter(
            new SerialBOVWriter<Cell>(
                    Selector<Cell>(&Cell::p, "p"),
                    "p",
                    outputFrequency));

    sim.addWriter(
            new SerialBOVWriter<Cell>(
                    Selector<Cell>(&Cell::v, "v"),
                    "v",
                    outputFrequency));
    sim.addWriter(
            new SerialBOVWriter<Cell>(
                    Selector<Cell>(&Cell::u, "u"),
                    "u",
                    outputFrequency));

#endif
#ifdef PMPI
    sim.addWriter(
            new BOVWriter<Cell>(
                    Selector<Cell>(&Cell::p, "p"),
                    "p",
                    outputFrequency));

    sim.addWriter(
            new BOVWriter<Cell>(
                    Selector<Cell>(&Cell::v, "v"),
                    "v",
                    outputFrequency));
    sim.addWriter(
            new BOVWriter<Cell>(
                    Selector<Cell>(&Cell::u, "u"),
                    "u",
                    outputFrequency));
#endif

    // running simulation
    sim.run();
}


int main(int argc, char *argv[])
{
#ifdef PMPI
    MPI_Init(&argc, &argv);
#endif
    std::cout << "*** poisson2d solving ***" << std::endl;
    runSimulations();

#ifdef PMPI
    MPI_Finalize();
#endif
    return 0;
}
