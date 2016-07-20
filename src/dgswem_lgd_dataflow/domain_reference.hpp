#ifndef DOMAIN_REFERENCE 
#define DOMAIN_REFERENCE

#include "../fname.h"
#include "../fortran_declarations.hpp"

std::vector<int> neighboringDomainIDs(void *size, void *dg, void *global)
{
    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    if (!size) {
        return std::vector<int>();
    }

    FNAME(get_neighbors_fort)(&size,
			      &dg,
			      &global,
			      neighbors_fort,
			      &numneighbors_fort);
    std::vector<int> ret(numneighbors_fort);
    std::copy(neighbors_fort, neighbors_fort + numneighbors_fort, ret.begin());

    return ret;
}

class DomainReference
{
public:
    const static int NANO_STEPS = 1;

    class API :
        public LibGeoDecomp::APITraits::HasUnstructuredTopology,
        public LibGeoDecomp::APITraits::HasNanoSteps<NANO_STEPS>,
	public LibGeoDecomp::APITraits::HasCustomMessageType<std::vector<double> >
    {};

    class FortranPointerWrapper
    {
    public:
        explicit
        FortranPointerWrapper(void *size, void *global, void *dg, void *nodalattr) :
            size(size),
	    global(global),
	    dg(dg),
	    nodalattr(nodalattr)
        {}
      
	~FortranPointerWrapper()
	{
	    if (size) {
		//	      std::cout << "CPP: about to call term_fort" << std::endl;
		FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    }
	}
      

	void *size;
	void *global;
	void *dg;
	void *nodalattr;

    };

    explicit
    DomainReference(int id = 0, 
		    void *size = 0, 
		    void *global = 0, 
		    void *dg = 0, 
		    void *nodalattr = 0, 
		    const std::vector<int>& neighbors_here = std::vector<int>()) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors(neighbors_here),
        id(id),
        timestep(1),
	rkstep(1)
    {}
	
    template<typename HOOD>
    void update(HOOD& hood, int nanoStep, int step)
    {		
        int globalNanoStep = step * NANO_STEPS + nanoStep;

	
	//	if (id == 0)
	//hpx::cout << "globalNanoStep = " << globalNanoStep << " timestep = " << timestep << std::endl;
	hpx::cout << "globalNanoStep = " << globalNanoStep << " timestep = " << timestep << " id = " << id << std::endl;
	

	// Retrieve variables 
	if (timestep > 1) {
	    // FNAME(cpp_vars_from_fort)(&domainWrapper->size,&rkstep,&timestep);
	    
	    // Place incoming buffers into our cell 
	    // Loop over neighbors
	    for (auto&& neighbor: neighbors) {
		std::vector<double> buffer_vector = hood[neighbor];
		int volume; // We don't really use this right now
		double buffer[MAX_BUFFER_SIZE];

		/* DEBUGGING
		   std::cout << "buffer_vector.size() = " << buffer_vector.size() << "\n";
		*/
		
		/* Set correct rk indices: This is required to place incoming data
		   in the correct places in the arrays */	
		int rkindex;       
		if (rkstep == 1) rkindex = 1;
		if (rkstep == 2) rkindex = 2;	
		
		// Make sure the buffer isn't empty, which it never should be
		if (buffer_vector.size() == 0) {
		    throw std::logic_error("Empty incoming buffer!");
		}	
		// Make sure buffer isn't larger than MAX_BUFFER_SIZE
		if (buffer_vector.size() > MAX_BUFFER_SIZE) {
		    throw std::logic_error("buffer_vector > MAX_BUFFER_SIZE!");
		}

		// Load buffer into c-style array
		for (int i=0; i<buffer_vector.size(); i++) {
		    buffer[i] = buffer_vector[i];
		}

		rkindex = 0;
		/*
		  if (id==0)
		    hpx::cout << "hpx_put_elems" << std::endl;
		*/
		FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
					  &neighbor,
					  &volume,
					  buffer,
					  &rkindex);
			 
	    } // End loop over neighbors
	    
	// Call dg_timestep_advance before first RK hydro step
	    if (rkstep == 1) { 

		int timestep_minus = timestep-1; // Because we are
		// actually calling advance for the *previous*
		// timestep.
		/*
		if (id==0)
		    hpx::cout << "dg_timestep_advance" << std::endl;
		*/
		FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
						&domainWrapper->dg,
						&domainWrapper->global,
						&domainWrapper->nodalattr,
						&timestep_minus
						);
	    }
	}
	    
	// Call dg_hydro_timestep_fort
	/*
		if (id==0)
		    hpx::cout << "dg_hydro_timestep, rkstep = " << rkstep << std::endl;
	*/
	FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
				      &domainWrapper->dg,
				      &domainWrapper->global,
				      &domainWrapper->nodalattr,
				      &timestep,
				      &rkstep
				      );

	// Loop over neighbors
	for (auto&& neighbor: neighbors) {
	    std::vector<double> send_buffer;
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];
	    int rkindex;
	    if (rkstep == 1) rkindex = 2;
	    if (rkstep == 2) rkindex = 3;
	   
	    /*
		if (id==0)
		    hpx::cout << "hpx_get_elems" << std::endl;
	    */
	    FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
				      &neighbor, // pointer to neighbor to send to
				      &volume,
				      buffer,
				      &rkindex);

	    for (int i=0; i<volume; i++) {
		send_buffer.push_back(buffer[i]);
	    }
	    
	    hood.send(neighbor, send_buffer);
	}
	    
	if (rkstep == 2) {
	    ++timestep;
	    --rkstep;
	} else {
	    ++rkstep;
	}
     
	// FNAME(cpp_vars_to_fort)(&domainWrapper->size,&rkstep,&timestep);

    } // End of update()

    
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	ar & id & timestep & rkstep; // DO THIS
    }

private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    int id;
    int timestep;
    int rkstep;
    std::vector<int> neighbors;
};

LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(DomainReference)
#endif
