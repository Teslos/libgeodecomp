// model 2d poisson equation
//#include<libgeodecomp.h>
#include<libgeodecomp/io/simpleinitializer.h>
#include<libgeodecomp/io/ppmwriter.h>
#include<libgeodecomp/io/simplecellplotter.h>
#include<libgeodecomp/parallelization/serialsimulator.h>


using namespace LibGeoDecomp;
//const double dx = 0.05;
//const double mu = 0.1;

class Cell {
public:
       // static MPI_Datatype MPIDataType;
	friend void runSimulations();

#ifdef PMPI
	static MPI_Datatype MPIDataType;
#endif

	class API :
			public APITraits::HasStencil<Stencils::VonNeumann<2, 1> >,
			public APITraits::HasCubeTopology<2>,
			public APITraits::HasCustomMPIDataType <Cell>,
			public APITraits::HasNanoSteps<2> {};

	enum State {LIQUID, WEST_NOSLIP, EAST_NOSLIP, SOUTH_NOSLIP, NORTH_BOUND};
//	inline explicit Cell() {};

	// initialize cells for begin simulation to zero values
	inline explicit Cell(double u = 0.0, double v = 0.0, double p = 0.0, State s = LIQUID) : state(s),
																							 u(u),v(v),p(p)
	{}

	/**
	 * defines the central difference operator for u in x-direction
	 * @param hood
	 */
	template<typename COORD_MAP>
	double u_xcent(const COORD_MAP &hood )
	{
		return (hood[Coord<2>(1,0)].u - hood[Coord<2>(-1, 0)].u)/(2.*dx);
	}

	/**
	 * defines the central difference operator for u in y-direction
	 * @param hood
	 */
	template<typename COORD_MAP>
	double u_ycent(const COORD_MAP &hood)
	{
		return (hood[Coord<2>(0, 1)].u - hood[Coord<2>(0, -1)].u)/(2.*dx);
	}
	/**
	 * defines central difference operator for v in y-direction
	 * @param hood
	 */
	template<typename COORD_MAP>
	double v_ycent(const COORD_MAP &hood)
	{
		return (hood[Coord<2>(0, 1)].v - hood[Coord<2>(0, -1)].v)/(2.*dx);
	}
	/**
	 * defines central difference operator for v in x-direction
	 * @param hood
	 */
	template<typename COORD_MAP>
	double v_xcent(const COORD_MAP& hood)
	{
		return (hood[Coord<2>(1, 0)].v - hood[Coord<2>(-1, 0)].v)/(2.*dx);
	}

	/**
	 * calculates the pressure poisson equation
	 * nit- is pseudo time variable for ensuring
	 * that Poisson calculation gives a divergence-free field.
	 * @param hood
	 * @return  a new pressure
	 */
    template<typename COORD_MAP>
	double pressure_poisson(const COORD_MAP &hood)
	{
		for(int it = 0; it < nit; ++it){
			p = 1./4.*(hood[Coord<2>(1, 0)].p + hood[Coord<2>(-1, 0)].p
					+hood[Coord<2>(0, 1)].p + hood[Coord<2>(0, -1)].p)
				-rho*dx*dx/4.*(1./dt * (  u_xcent(hood) + v_ycent(hood) )
							   -u_xcent(hood)* u_xcent(hood)
					           -2.*u_ycent(hood)*v_xcent(hood)
							   -v_ycent(hood)*v_ycent(hood));
			//hood[Coord<2>(0, 0)].p = p;  // update the stencil with the new value
				// set boundary condition for pressure
		switch(state) {
			case WEST_NOSLIP:
				updateWestPress(hood);
				break;
			case EAST_NOSLIP:
				updateEastPress(hood);
				break;
			case SOUTH_NOSLIP:
				updateSouthPress(hood);
				break;
			case NORTH_BOUND:
				updateNorthPress(hood);
				break;
		}

		}
		return p;
	}

	template<typename COORD_MAP>
	double old_u(const COORD_MAP &hood, int i) {
		return hood[Coord<2>(0, 0)].u[i];
	}

	template<typename COORD_MAP>
	double old_v(const COORD_MAP &hood, int i) {
		return hood[Coord<2>(0, 0)].v[i];
	}

	template<typename COORD_MAP>
	double old_p(const COORD_MAP &hood, int i) {
		return hood[Coord<2>(0, 0)].p[i];
	}

	/**
     * Defines laplace operator for u-momentum of fluid flow
     *
     */
	template<typename COORD_MAP>
	double lap_uvel(const COORD_MAP &hood) {
		return hood[Coord<2>(0, -1)].u
			   + hood[Coord<2>(-1, 0)].u
			   + hood[Coord<2>(1, 0)].u
			   + hood[Coord<2>(0, 1)].u
			   - 4 * hood[Coord<2>(0, 0)].u;
	}

	/**
     * Defines laplace operator for v-momentum of fluid flow
     *
     */
	template<typename COORD_MAP>
	double lap_vvel(const COORD_MAP &hood) {
		return hood[Coord<2>(0, -1)].v
			   + hood[Coord<2>(-1, 0)].v
			   + hood[Coord<2>(1, 0)].v
			   + hood[Coord<2>(0, 1)].v
			   - 4 * hood[Coord<2>(0, 0)].v;
	}

    /**
     *  Defines the backward Euler operator for u-velocity
     */
	template<typename COORD_MAP>
	double u_back(const COORD_MAP &hood) {
		return hood[Coord<2>(0, 0)].u - hood[Coord<2>(-1, 0)].u;
	}

    /**
     *  Defines the backward Euler operator for v-velocity
     */
	template<typename COORD_MAP>
	double v_back(const COORD_MAP &hood) {
		return hood[Coord<2>(0, 0)].v - hood[Coord<2>(0, -1)].v;
	}
    /**
     * Defines the backward Euler operator for p- pressure
     * @tparam COORD_MAP
     * @param hood
     * @return
     */
    template<typename COORD_MAP>
    double p_xback(const COORD_MAP& hood)
    {
        return hood[Coord<2>(0, 0)].p - hood[Coord<2>(-1, 0)].p;
    }
    template<typename COORD_MAP>
    double p_yback(const COORD_MAP& hood)
    {
        return hood[Coord<2>(0, 0)].p - hood[Coord<2>(0, -1)].p;
    }
    template<typename COORD_MAP>
    void updateWestNoSlip(const COORD_MAP &hood)
    {
        u = 0.;
        v = 0.;
    }

    template<typename COORD_MAP>
    void updateEastNoSlip(const COORD_MAP &hood)
    {
        u = 0.;
        v = 0.;
    }

    template<typename COORD_MAP>
    void updateSouthNoSlip(const COORD_MAP &hood)
    {
       u = 0.;
       v = 0.;
    }

    template <typename COORD_MAP>
    void updateNorthBound(const COORD_MAP &hood)
    {
       u = 1.0;
       v = 0.;
    }

    template <typename COORD_MAP>
    void updateWestPress(const COORD_MAP &hood)
    {
	    p = hood[Coord<2>(-1,0)].p;  // dp/dx = 0 at x = W
    }

    template <typename COORD_MAP>
    void updateEastPress(const COORD_MAP &hood)
    {
	    p = hood[Coord<2>(1,0)].p; // dp/dx = 0 at x = E
    }

    template <typename COORD_MAP>
    void updateNorthPress(const COORD_MAP &hood)
    {
	    p = 0.;    // p = 0 at y = N
    }

    template <typename COORD_MAP>
    void updateSouthPress(const COORD_MAP &hood)
    {
	    p = hood[Coord<2>(0, 1)].p; // dp/dy = 0 at y = S
    } 

	/**
	 * updates the solution using the stencil
	 */
	template<typename COORD_MAP>
	void update(const COORD_MAP& hood, const unsigned& nanostep)
	{

        // solve pressure poisson equation
        p = pressure_poisson(hood);
	
        u = hood[Coord<2>(0, 0)].u;
        v = hood[Coord<2>(0, 0)].v;

	// solve for u-momentum
        u = u - u*dt/dx*u_back<COORD_MAP>(hood) - v*dt/dx*u_back<COORD_MAP>(hood)
                -dt/(2.*rho*dx)*p_xback(hood)+mu*dt/(dx*dx)*lap_uvel(hood);
        // solve for v-momentum
        v = v - u*dt/dx*v_back<COORD_MAP>(hood) - v*dt/dx*v_back<COORD_MAP>(hood)
                -dt/(2.*rho*dx)*p_yback<COORD_MAP>(hood)+mu*dt/(dx*dx)*lap_vvel<COORD_MAP>(hood);

        // set boundary conditions
        switch (state) {
            case WEST_NOSLIP:
                updateWestNoSlip(hood);
                break;
            case EAST_NOSLIP:
                updateEastNoSlip(hood);
                break;
            case SOUTH_NOSLIP:
                updateSouthNoSlip(hood);
                break;
            case NORTH_BOUND:
                updateNorthBound(hood);
                break;


        }
	}


	double u;
	double p;
	double v;
	State state;
};





		



