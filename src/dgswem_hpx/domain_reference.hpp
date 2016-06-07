#ifndef DOMAIN_REFERENCE 
#define DOMAIN_REFERENCE

#include "../fname.h"
#include "../fortran_declarations.hpp"

//Defining settings for SELL-C-q
static const std::size_t MATRICES = 1;
static const int C = 4;
static const int SIGMA = 1;

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
    class API :
        public LibGeoDecomp::APITraits::HasUnstructuredTopology,
	public LibGeoDecomp::APITraits::HasSellType<double>,
        public LibGeoDecomp::APITraits::HasSellMatrices<MATRICES>,
        public LibGeoDecomp::APITraits::HasSellC<C>,
        public LibGeoDecomp::APITraits::HasSellSigma<SIGMA>
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
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0, const std::vector<int>& neighbors = std::vector<int>() ) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighbors),
        id(id),
        timestep(0),
	rkstep(1),
	update_step(true),
	exchange_step(false),
	advance_step(false)
    {}

    template<typename HOOD>
    void update(const HOOD& hood, int nanoStep)
    {
      
	*this = hood[hood.index()];

	std::cout << "LGD update timestep =  " << timestep << " id = " << id << " ";
	
	if (update_step) {
	    std::cout << "update_step" << std::endl;
	} else if (exchange_step) {
	    std::cout << "exchange_step" << std::endl;
	} else if (advance_step) {
	    std::cout << "advance_step" << std::endl;
	} else {
	    std::cout << "LGD update, error! timestep = " << timestep << std::endl;
	}
	
	if (timestep != 0) {
	    
	    if (update_step) {
		
		// Clear output buffer map
		output_buffer.clear();
		
		std::cout << "updating (domain_id = " << id
			  << ", timestep = " << timestep
			  << ", rkstep = " << rkstep
			  << ", locality = " << hpx::find_here()
			  << ")...\n";
	      
	      
		std::cout << "domainWrapper->size = " << domainWrapper->size << std::endl;
		if (domainWrapper->size == 0) {
		    throw std::logic_error("bad domain pointer!");
		}
		
		FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
					      &domainWrapper->dg,
					      &domainWrapper->global,
					      &domainWrapper->nodalattr,
					      &timestep,
					      &rkstep
					      );
		
		// Fill output buffer
		//Loop over neighbor domains
		for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
		    int neighbor_here = neighbors_here[neighbor];
		    int volume;
		    double buffer[MAX_BUFFER_SIZE];
		    int rkstep_fort = 0;
		    FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
					      &neighbor_here, // pointer to neighbor to send to
					      &volume,
					      buffer,
					      &rkstep_fort);
		    
		    // Pack buffer into std vector
		    std::vector<double> buffer_vector;
		    for (int i=0; i<volume; i++) {
			buffer_vector.push_back(buffer[i]);
		    }
		    
		    std::pair<std::_Rb_tree_iterator<std::pair<const int, std::vector<double> > >, bool> outval = output_buffer.insert(std::map<int, std::vector<double> >::value_type(neighbor_here, buffer_vector));
		    if (std::get<1>(outval)) {
			//std::cout << "Insert successful!" << std::endl;
		    } else {
			//FIXME: add exception handling here
			std::cout << "Insert not successful!" << std::endl;
		    }
		    
		}   
		update_step = false;
		exchange_step = true;
		advance_step = false;
	    } else if (exchange_step) {
		std::cout << "exchange step stuff being called" << std::endl;
		
		// Boundary exchange
		//Loop over neighbors
		for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
		    int neighbor_here = neighbors_here[neighbor];
		    int volume;
		    double buffer[MAX_BUFFER_SIZE];
		    
		    // Unpack buffer from neighbor
		    std::vector<double> buffer_vector = hood[neighbor_here].output_buffer.at(id);
		    for (int i=0; i<buffer_vector.size(); i++) {
			buffer[i] = buffer_vector[i];
		    }
		    
		    // Put elements into our own subdomain	
		    int rkstep_fort = 0;
		    FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
					      &neighbor_here,
					      &volume,
					      buffer,
					      &rkstep_fort);
		    
		}// end loop over neighbors
		
		std::cout << "about to set exchange_step to false" << std::endl;
		exchange_step = false;
		
		if (rkstep == 2) {
		    advance_step = true;
		    update_step = false;
		    rkstep = 1;
		} else {
		    ++rkstep;
		    update_step = true;
		    advance_step = false;
		}
		
		
		} else if (advance_step) {	
		    std::cout << "advancing domain " << id << " at timestep " << timestep <<std::endl;
		    FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
						    &domainWrapper->dg,
						    &domainWrapper->global,
						    &domainWrapper->nodalattr,
						    &timestep
						    );
		    ++timestep;
		    
		    update_step = true;
		    exchange_step = false;
		    advance_step = false;
		} else {
		    std::cout << "Something went wrong!" << std::endl;
		}
		
	    } else {
		++timestep;				
	    }
	}

	// Serialization function
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	ar & neighbors_here & output_buffer & id & timestep & rkstep & update_step & exchange_step & advance_step;
    }

private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    std::vector<int> neighbors_here;
    std::map<int,std::vector<double> > output_buffer;
    int id;
    int timestep;
    int rkstep;
    bool update_step;
    bool exchange_step;
    bool advance_step;
    };


LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(DomainReference)
#endif
