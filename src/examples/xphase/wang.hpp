#include<vector>
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
const double gc     = 0.05; 


// model allen cahn

class Cell
{
public:
        friend void runSimulation();
        friend class RigidMotions;
#ifdef PMPI
        static MPI_Datatype MPIDataType;
#endif
    static vector<double> staticData;        
    class API :
        public APITraits::HasFixedCoordsOnlyUpdate,
        public APITraits::HasStencil<Stencils::VonNeumann<3, 1> >,
        public APITraits::HasTorusTopology<3>,
        public APITraits::HasCustomMPIDataType<Cell>,
        public APITraits::HasNanoSteps<2>,
        public APITraits::HasStaticData<vector<double>>
    {};

    inline explicit Cell() {}
    
    inline explicit Cell(vector<double> v)
    {
        rho = 0.0;   // initialize mass density
        //vel_sumx = vel_sumy = vel_sumz = 0.0; // initialize velocities 
        for(int i = 0; i < nf; ++i)
        {
            eta[i] = v[i];  // vector of crystalline orientations
            //rho += eta[i];
		    rho = 1.0; // initialize mass density to one ( ivt )
        }
		// initialize forces and velocities to zero
		memset(f, 0, sizeof(f[0][0][0]) *nf*nf*3 );
        memset(vel, 0, sizeof(vel[0][0]) *nf*3);
    }
    
	template<typename MAP>
	inline __host__ __device__
	double grainb(const MAP& hood, int i, int j)
	{
		return hood[FixedCoord< 0, 0, 0 >()].eta[i] * hood[FixedCoord< 0, 0, 0 >()].eta[j];
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
	double* nabla_eta(const MAP& hood, int i) 
	{
       static double v[3];
       v[0] =  hood[FixedCoord< 1, 0, 0>()].eta[i] - hood[FixedCoord< -1, 0, 0 >()].eta[i];      
       v[1] =  hood[FixedCoord< 0, 1, 0>()].eta[i] - hood[FixedCoord< 0, -1 ,0 >()].eta[i];
       v[2] =  hood[FixedCoord< 0, 0, 1>()].eta[i] - hood[FixedCoord< 0, 0, -1 >()].eta[i];
	   return v;
	}	
    /**
     * Defines nabla operator for velocity of the particle i
     * \partial_x v_x 
     * \partial_y v_y
     * \partial_z v_z
     * @WARNING is this correct maybe should divide by 2 because deltax is unit. 
     */
    template<typename MAP>
    inline __host__ __device__
    double* nabla_vel(const MAP& hood, int i)
    {
        static double gradv[3];
        gradv[0] = hood[FixedCoord< 1, 0, 0 >()].vel[i][0] - hood[FixedCoord< -1, 0, 0 >()].vel[i][0];
        gradv[1] = hood[FixedCoord< 0, 1, 0 >()].vel[i][1] - hood[FixedCoord< 0, -1, 0 >()].vel[i][1];
        gradv[2] = hood[FixedCoord< 0, 0, 1 >()].vel[i][2] - hood[FixedCoord< 0, 0, -1 >()].vel[i][2];
        return gradv;
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

            // set the forces to zero
            memset(f, 0, sizeof(double)*nf*nf*3);
            // update the particles 
            for(int i = 0; i < nf; ++i)
            {
                double oe = old_eta<MAP>(hood, i); // take previous time step eta
                
                phi  += oe * oe;
                eta3 += oe * oe * oe;
            }
            
            double orho = old_rho<MAP>(hood); // take previous time step rho

            // here calculate the forces coefficients to account
			// for rigid motion of the particles
			for (int i = 0; i < nf; ++i) 
			{
			   for (int j = i+1; j < nf; ++j) 
               {
                   double gij = grainb<MAP>(hood,i,j);    
                   if( gij > gc ) 
                   {
                       // @ factor 0.5 comes from the central differences for gradient operator
                       // added template <MAP> designator--ivt 15.05.19
                       f[i][j][0] = kappa * (orho - rho_o) * 0.5*(nabla_eta(hood,i)[0] - nabla_eta(hood,j)[0]); 
                       f[i][j][1] = kappa * (orho - rho_o) * 0.5*(nabla_eta(hood,i)[1] - nabla_eta(hood,j)[1]); 
                       f[i][j][2] = kappa * (orho - rho_o) * 0.5*(nabla_eta(hood,i)[2] - nabla_eta(hood,j)[2]); 
		       		   //cout << "Forces on (" << i <<","<<j<<"):"<<"("<<f[i][j][0] << ","<<f[i][j][1]<<","<<f[i][j][2]<<")"<<endl;
                       f[j][i][0]=-f[i][j][0];
                       f[j][i][1]=-f[i][j][1];
                       f[j][i][2]=-f[i][j][2];
                   }    
               }
            }

            memset(sumf, 0.0, sizeof(sumf[0][0])*nf*3); // @@ CLEAN not used anymore

            // reset the global velocity field to zero
            //vel_sumx = 0.0;
            //vel_sumy = 0.0;
            //vel_sumz = 0.0;
            // forces on the particle i
            for(int i = 0; i < nf; ++i)
            {
                double sumfx = 0.0;
                double sumfy = 0.0;
                double sumfz = 0.0;
                for(int j=0; j < nf; ++j)
                {
                    if(i != j) 
                    {
                        /*
                        sumf[i][0] += f[i][j][0];
                        sumf[i][1] += f[i][j][1];
                        sumf[i][2] += f[i][j][2];*/
                        sumfx += f[i][j][0];
                        sumfy += f[i][j][1];
                        sumfz += f[i][j][2];
                    }
                }
                if (sumfx > 0.0 || sumfy > 0.0 || sumfz > 0.0)
                    //cout << "Sumf "<<i<<" (x,y,z):" << sumfx << "," << sumfy << "," << sumfz << endl; 
                /*vel[i][0] = mt * sumf[i][0] / staticData[i] * eta[i];
                vel[i][1] = mt * sumf[i][1] / staticData[i] * eta[i];
                vel[i][2] = mt * sumf[i][2] / staticData[i] * eta[i];*/
                //cout << "Volume of the particle " << i << ":" << staticData[i] << endl; 
                if (staticData[i] > 0.0) {
                    vel[i][0] = mt * sumfx / staticData[i] * eta[i];
                    vel[i][1] = mt * sumfy / staticData[i] * eta[i];
                    vel[i][2] = mt * sumfz / staticData[i] * eta[i];

                    // here we add all particles velocities to get
                    // global velocities for rho field.
                    /*
                    vel_sumx += vel[i][0];
                    vel_sumy += vel[i][1];
                    vel_sumz += vel[i][2];
                    cout << "vel "<<i<<" (x,y,z):" << vel_sumx << "," << vel_sumy << "," << vel_sumz << endl; 
                    */
                }
            }

            //phi += orho * orho;
                
            for(int i = 0; i < nf; ++i)
            {
                double oe = old_eta<MAP>(hood, i);
                
                chi[i]  = - k_eta * lap_eta<MAP>(hood, i) 
                    // adding advection term to model rigid motion
                    // of the particles -- ivt 15.05.19
                    // 0.5 factor is comming from the central differences used
                    // for nabla_eta and nabla_vel operators.
                          +0.5*(nabla_eta(hood,i)[0] * vel[i][0] + nabla_eta(hood,i)[1] * vel[i][1] 
                                  + nabla_eta(hood,i)[2] * vel[i][2] + oe * nabla_vel(hood,i)[0] + oe * nabla_vel(hood,i)[1] + oe * nabla_vel(hood,i)[2])
                          + 12.0 * B * (1 - orho) * oe
                          - 12.0 * B * (2 - orho) * oe * oe
                          + 12.0 * B * oe * phi;
            } 
            // changes to original Wang paper 
	        // -4.0 * B * eta3 is +4*B*eta3 and additional +2.0*B*orho term 
            // @@ NOTE added the new advection term for the velocity--ivt 19.05.19
            double velx = 0.0;
            double vely = 0.0;
            double velz = 0.0;
            for (int i = 0; i < nf; ++i) {
                velx += vel[i][0];
                vely += vel[i][1];
                velz += vel[i][2];
            }
            psi = - k_rho * lap_rho<MAP>(hood)
                  - (orho * velx + orho * vely + orho * velz) 
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
    double f[nf][nf][3]; 
    double sumf[nf][3]; 
    double vel[nf][3]; 
    //double vel_sumx, vel_sumy, vel_sumz;
    double rho, psi, phi, eta[nf], chi[nf], vr[nf], vt[nf], eta3;
};

#ifdef PMPI
MPI_Datatype Cell::MPIDataType;
#endif


