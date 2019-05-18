// phase field model

#include <libgeodecomp.h>
#define PSERIAL 1
#ifdef PCUDA
#include <libgeodecomp/parallelization/cudasimulator.h>
#endif

#ifdef PMPI
#include <libgeodecomp/parallelization/hiparsimulator.h>
#include <libgeodecomp/loadbalancer/tracingbalancer.h>
#include <libgeodecomp/loadbalancer/noopbalancer.h>
#include <libgeodecomp/geometry/partitions/recursivebisectionpartition.h>
#include <libgeodecomp/io/bovwriter.h>
#endif

#ifdef PMGPU
#include <libgeodecomp/parallelization/hiparsimulator/cudastepper.h>
#endif

#include <iostream>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace LibGeoDecomp;

class RigidMotions;
// constants

const int nf = 12;
//const int nx = 1152/2;
//const int ny = 1152/2;
const int nx = 256;
const int ny = 256;
const int nz = 1;

int outputFrequency = 1000;
int numSteps = 100000;

const double dt = 0.001;
// model

//#include "allen_cahn.hpp"
//#include "cahn_hilliard.hpp"
#include "wang.hpp"

// init
vector<double> Cell::staticData = {0};
class CellInitializer : public SimpleInitializer<Cell>
{
public:
    CellInitializer() : SimpleInitializer<Cell>(Coord<3>(nx, ny, nz), numSteps)
    {}

 	void readparticles(	vector< vector<double> > &particles )
	{
		using namespace boost;
		using namespace std;
		cout << "Reading particles" << endl;
		string data("26particles.txt");
		ifstream in(data.c_str());
		if (!in) {
			cerr << "Can't open the input file " << data << endl;
			exit(1);
		}


		typedef tokenizer<escaped_list_separator<char> > Tokenizer;
		string line;
		vector<string> vec;
		int l = 0;
		while(getline(in,line))
		{
			Tokenizer tok(line);
			vec.assign(tok.begin(),tok.end());
			//vector <double> v = {atof(vec[5].c_str()),atof(vec[6].c_str()),atof(vec[8].c_str())};
			particles[l][0] = atof(vec[0].c_str());
			particles[l][1] = atof(vec[1].c_str());
			particles[l][2] = atof(vec[2].c_str());
			cout << "vector of coor:" << particles[l][0] << " " << particles[l][1] << ',' << particles[l][2]<<  endl;
			l++;
		}

	}
		   	
    virtual void grid(GridBase<Cell, 3> *ret)
    {
	vector< vector<double> > particles(30, vector<double>(3, 0.0));	
        boost::mt19937 rng;
        boost::normal_distribution<> nd(0.0, 0.5);
        boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > var_nor(rng, nd);
        
        CoordBox<3> rect = ret->boundingBox();
	// read particles from external file
        readparticles(particles);
        double size_img = 1152/2;	
        for (int z = 0; z < nz; ++z) 
            for (int y = 0; y < ny; ++y) 
                for (int x = 0; x < nx; ++x) 
                {
                    Coord<3> c(x, y, z);
                    vector<double> v(nf);
                    
                    //for(int i = 0; i < nf; ++i)
                    //    v[i] = var_nor();
                    
                    //if (rect.inBounds(c))
                    //    ret->set(c, Cell(v));
                        
                    if (rect.inBounds(c))
                    {
						// sets the particles inside domain
                        for(int i = 0; i < nf; ++i)
                        {
							/*
                            double dx = double(x - nx/2)/double(nx);
                            double dy = double(y - ny/2)/double(ny);
                            double dz = double(z - nz/2)/double(nz);
                            double r2 = dx*dx + dy*dy + dz*dz;

                            // reduced the internal pore from 0.04 to 0.02 
                            if(r2 > 0.02)
                              v[i] = var_nor() + 0.12;
                            else
                              v[i] = 0;
                      		*/
							double xc = particles[i][0]/size_img;
							double yc = particles[i][1]/size_img;
							double rc = particles[i][2]/size_img;	
							//cout << "xc: " << xc << "yc: " << yc << "rc: " << rc << endl;
							double dx = double(x - xc*nx)/double(nx);
							double dy = double(y - yc*ny)/double(ny);
							double dz = double(z - nz/2)/double(nz);
							dz = 0.0;
							double r2 = dx*dx + dy*dy + dz*dz;
							// scale radius of particles
							//double scale = 5.0/160.0;
							//rc *= scale;
							if (r2 < rc*rc) {
								   //v[i] = var_nor() + 0.12;
								   v[i] = 2.0/3.0 * M_PI ;
                            	   				   ret->set(c, Cell(v));
							} else {	   
								   v[i] = 0.0;
							}	   
                            //ret->set(c, Cell(v));
                        }
                    }   
                }
    }
};
//
// Steerer class
// this class is used to run the rigid motion 
// of the sintering particles.
//
//
class RigidMotions : public Steerer<Cell> 
{
	public:
		using Steerer<Cell>::CoordType;
		using Steerer<Cell>::GridType;
		using Steerer<Cell>::Topology;
	RigidMotions(const unsigned ioPeriod): Steerer<Cell>(ioPeriod), vol(30),xc(30),yc(30) 
	{}
	void nextStep(
		GridType *grid,
		const Region<Topology::DIM>& validRegion,
		const CoordType& globalDimensions,
		unsigned step,
		SteererEvent event,
		std::size_t rank,
		bool lastCall,
		SteererFeedback *feedback)
	{
		Coord<3> dim = grid->dimensions();
		//cout << "Dim: " << "("<< dim.x() <<","<< dim.y()<<"," << dim.z() <<")"<< endl;
		for (int i = 0; i < nf; ++i) vol[i] = 0.0; 
		// calculate volume of the particles
		// and center of gravity for the particles
		for (int z = 0; z < nz; ++z)
           	for (int y = 0; y < ny; ++y)
              	for (int x = 0; x < nx; ++x)
                {
                    Coord<3> c(x,y,z);
                    Cell cell = grid->get(c);
                    for(int i = 0; i < nf; ++i) {
                        vol[i] += cell.eta[i];
                        xc[i]  += x * cell.eta[i]; 
                        yc[i]  += y * cell.eta[i];
                    }	
                }   
		for (int i = 0; i < nf; ++i) {
			//cout << "Volume of particle first eval." << i << ":" << vol[i] << endl;
			xc[i] /= vol[i];
			yc[i] /= vol[i];
		}


		// write the current value of volume for particles
		for (int i = 0; i < nf; ++i) {
			//xc[i] /= vol[i];
			//yc[i] /= vol[i];
			//cout << "Volume of particle " << i << ":" << vol[i] << endl;
			//cout << "(xc,yc): " << xc[i] << "," << yc[i] << endl;
		}
        feedback->setStaticData(vol);
													
	}

	private:
	vector<double> vol;
	vector<double> xc, yc;
};
// run
    
void runSimulation()
{
#ifdef PMPI
    MPI_Aint displacements[] = { 0 };
    MPI_Datatype memberTypes[] = { MPI_CHAR };
    int lengths[] = { sizeof(Cell) };
    MPI_Type_create_struct(1, lengths, displacements, memberTypes, &Cell::MPIDataType);
    MPI_Type_commit(&Cell::MPIDataType);
#endif
    
    // init

    CellInitializer *init = new CellInitializer();
    
    // simulators
    
#ifdef PSERIAL
    SerialSimulator<Cell> sim(init);
#endif

#ifdef PCUDA
    CudaSimulator<Cell> sim(init, Coord<3>(64, 8, 1));
#endif

#ifdef PMPI
#ifdef PMGPU
    HiParSimulator::HiParSimulator
        <Cell, RecursiveBisectionPartition<3>, HiParSimulator::CUDAStepper<Cell> > sim(
            init,
            MPILayer().rank() ? 0 : new NoOpBalancer(),
            numSteps,
            1);
#else
    HiParSimulator<Cell, RecursiveBisectionPartition<3> > sim(
        init,
        MPILayer().rank() ? 0 : new NoOpBalancer(),
        numSteps,
        1);
#endif
#endif

    // writers
    
#ifdef PMPI
    if (MPILayer().rank() == 0)
        sim.addWriter(new TracingWriter<Cell>(outputFrequency, init->maxSteps()));
        
    sim.addWriter(
        new BOVWriter<Cell>(
            Selector<Cell>(&Cell::phi, "phi"),
            "phi",
            outputFrequency));
            
    sim.addWriter(
        new BOVWriter<Cell>(
            Selector<Cell>(&Cell::rho, "rho"),
            "rho",
            outputFrequency));
            
#else
    sim.addWriter(new TracingWriter<Cell>(outputFrequency, init->maxSteps()));
    
    sim.addWriter(
        new SerialBOVWriter<Cell>(
            Selector<Cell>(&Cell::phi, "phi"),
            "phi",
            outputFrequency));
            
    sim.addWriter(
        new SerialBOVWriter<Cell>(
            Selector<Cell>(&Cell::rho, "rho"),
            "rho",
            outputFrequency));
            
#endif
    cout << "running simulation" << endl;
    sim.addSteerer(new RigidMotions(10));
    sim.run();
}

// main

int main(int argc, char **argv)
{
    //simParamsHost.initParams(argc, argv);
    //cudaSetDevice(simParamsHost.cudaDevice);
#ifdef PMPI
    MPI_Init(&argc, &argv);
#endif    

    cout << "** xphase **" << endl;
    
    runSimulation();
    
#ifdef PMPI
    MPI_Finalize();
#endif

    return 0;
}
