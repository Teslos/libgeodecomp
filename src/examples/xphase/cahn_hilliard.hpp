// model cahn hilliard

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
        for(int i = 0; i < nf; ++i)
            eta[i] = v[i];
    }
    
    template<typename MAP>
    inline __host__ __device__
    double old_eta(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0,  0>()].eta[i];
    }
        
    template<typename MAP>
    inline __host__ __device__
    double old_psi(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0,  0>()].psi[i];
    }
    
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
    
    template<typename MAP>
    inline __host__ __device__
    double lap_psi(const MAP& hood, int i)
    {
        return hood[FixedCoord< 0,  0, -1>()].psi[i]
             + hood[FixedCoord< 0, -1,  0>()].psi[i]
             + hood[FixedCoord<-1,  0,  0>()].psi[i]
             + hood[FixedCoord< 1,  0,  0>()].psi[i]
             + hood[FixedCoord< 0,  1,  0>()].psi[i]
             + hood[FixedCoord< 0,  0,  1>()].psi[i]
         - 6 * hood[FixedCoord< 0,  0,  0>()].psi[i];
    }
    
    template<typename MAP>
    inline __host__ __device__
    double delta_eta(const MAP& hood, int i)
    {
        return lap_eta<MAP>(hood, i) 
             + old_eta<MAP>(hood, i)
             + old_eta<MAP>(hood, i) * old_eta<MAP>(hood, i) * old_eta<MAP>(hood, i)
       - 2.0 * old_eta<MAP>(hood, i) * phi;
    }

    template<typename MAP>
    __host__ __device__
    void update(const MAP& hood, const unsigned& nanoStep)
    {
        if(nanoStep == 0)
        {
            phi = 0.0;
            
            for(int i = 0; i < nf; ++i)
                phi += old_eta<MAP>(hood, i) * old_eta<MAP>(hood, i);
                
            for(int i = 0; i < nf; ++i)
                psi[i] = delta_eta<MAP>(hood, i);
        }
        else if(nanoStep == 1)
        {
            phi = hood[FixedCoord< 0, 0, 0>()].phi;
            
            for(int i = 0; i < nf; ++i)
                eta[i] -= dt * lap_psi<MAP>(hood, i);
        }
        
    }
    
    template<typename MAP>
    __host__ __device__
    void updateCUDA(const MAP& hood, const unsigned& nanoStep)
    {
        update(hood, nanoStep);
    }
    
    double phi, eta[nf], psi[nf];
};

#ifdef PMPI
MPI_Datatype Cell::MPIDataType;
#endif


