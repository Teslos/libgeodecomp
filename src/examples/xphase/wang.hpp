const double A      = 16.0;
const double B      =  1.0;
const double D      =  0.1;
const double L      = 10.0;
const double k_rho  = 50.0; 
const double k_eta  =  5.0;
const double mr     = 500.0;
const double mt     = 1.0;
const double kappa  = 100.0;
const double rho_o  = 0.9816;


// model allen cahn

class Cell
{
public:
        friend void runSimulation();
#ifdef PMPI
        static MPI_Datatype MPIDataType;
#endif
        
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >,
        public APITraits::HasTorusTopology<3>,
        public APITraits::HasCustomMPIDataType<Cell>,
        public APITraits::HasNanoSteps<2>
    {};

    inline explicit Cell() {}
    
    inline explicit Cell(vector<double> v)
    {
        rho = 0.0;   // initialize mass density
        for(int i = 0; i < nf; ++i)
        {
            eta[i] = v[i];  // vector of crystalline orientations
            //rho += eta[i];
			rho = 1.0; // initialize mass density to one ( ivt )
        }
    }
    
    template<typename MAP>
    inline __host__ __device__
    double old_eta(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0,  0>()].eta[i];
    }
        
    template<typename MAP>
    inline __host__ __device__
    double old_chi(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0,  0>()].chi[i];
    }

    template<typename MAP>
    inline __host__ __device__
    double old_rho(const MAP& hood)
    {
        return hood[FixedCoord< 0,  0,  0>()].rho;
    }
   
    /**
     * Defines laplace operator for eta- multicomponent order parameter
     * that is zero in liquid and has value {0,1,0,....,0} for complete solid rho=1.
     */ 
    template<typename MAP>
    inline __host__ __device__
    double lap_eta(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0, -1>()].eta[i]
             + hood[FixedCoord< 0, -1,  0>()].eta[i]
             + hood[FixedCoord<-1,  0,  0>()].eta[i]
             + hood[FixedCoord< 1,  0,  0>()].eta[i]
             + hood[FixedCoord< 0,  1,  0>()].eta[i]
             + hood[FixedCoord< 0,  0,  1>()].eta[i]
         - 6 * hood[FixedCoord< 0,  0,  0>()].eta[i];
    }

	/**
	 * Defines nabla operator for eta- multicomponent order parameter
	 * that is zero in liquid and has value {0,1,0,....,0} for complete solid rho=1.
	 */
	template<typename MAP>
	inline __host__ __device__
	double nabla_eta(const MAP& hood, int i) 
	{
			
	}	
    
    template<typename MAP>
    inline __host__ __device__
    double lap_psi(const MAP& hood)
    {
        return hood[FixedCoord< 0,  0, -1>()].psi
             + hood[FixedCoord< 0, -1,  0>()].psi
             + hood[FixedCoord<-1,  0,  0>()].psi
             + hood[FixedCoord< 1,  0,  0>()].psi
             + hood[FixedCoord< 0,  1,  0>()].psi
             + hood[FixedCoord< 0,  0,  1>()].psi
         - 6 * hood[FixedCoord< 0,  0,  0>()].psi;
    }

    /**
     * Defines laplace operator for rho order parameter that is mass density
     * 0 in liquid and 1 in solid.
     */
    template<typename MAP>
    inline __host__ __device__
    double lap_rho(const MAP& hood)
    {
        return hood[FixedCoord< 0,  0, -1>()].rho
             + hood[FixedCoord< 0, -1,  0>()].rho
             + hood[FixedCoord<-1,  0,  0>()].rho
             + hood[FixedCoord< 1,  0,  0>()].rho
             + hood[FixedCoord< 0,  1,  0>()].rho
             + hood[FixedCoord< 0,  0,  1>()].rho
         - 6 * hood[FixedCoord< 0,  0,  0>()].rho;
    }

    /**
     * updates the solution using the stencil
     *
     */
    template<typename MAP>
    __host__ __device__
    void update(const MAP& hood, const unsigned& nanoStep)
    {
        if(nanoStep == 0)
        {
            phi   = 0.0;
            eta3  = 0.0;
            // update the particles 
            for(int i = 0; i < nf; ++i)
            {
                double oe = old_eta<MAP>(hood, i); // take previous time step eta
                
                phi  += oe * oe;
                eta3 += oe * oe * oe;
            }
            
            double orho = old_rho<MAP>(hood); // take previous time step rho
            
            //phi += orho * orho;
                
            for(int i = 0; i < nf; ++i)
            {
                double oe = old_eta<MAP>(hood, i);
                
                chi[i]  = - k_eta * lap_eta<MAP>(hood, i) 
                          + 12.0 * B * (1 - orho) * oe
                          - 12.0 * B * (2 - orho) * oe * oe
                          + 12.0 * B * oe * phi;
            } 
            // changes to original Wang paper 
	        // -4.0 * B * eta3 is +4*B*eta3 and additional +2.0*B*orho term 
            psi = - k_rho * lap_rho<MAP>(hood)
                  + 2.0 * A * orho 
                  - 6.0 * A * orho * orho
                  + 4.0 * A * orho * orho * orho
		  		  + 2.0 * B * orho
                  - 6.0 * B * phi
                  + 4.0 * B * eta3;
                  

        }
        else if(nanoStep == 1)
        {
            //phi = hood[FixedCoord< 0, 0, 0>()].phi;
            
            phi   = 0.0;
            
            for(int i = 0; i < nf; ++i)
            {
                eta[i] += - dt * old_chi<MAP>(hood, i);
                
                phi  += eta[i] * eta[i];
            }

            rho += D * dt * lap_psi<MAP>(hood);
        }
        
    }
    
    template<typename MAP>
    __host__ __device__
    void updateCUDA(const MAP& hood, const unsigned& nanoStep)
    {
        update(hood, nanoStep);
    }
    
    double rho, psi, phi, eta[nf], chi[nf], vr[nf], vt[nf], eta3;
};

#ifdef PMPI
MPI_Datatype Cell::MPIDataType;
#endif


