#include <iostream>
#include <vector>
#include <cstdlib>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

//#include "fname.h"
//#include "fortran_declarations.hpp" // This includes parameters defined

#define MAX_DOMAIN_NEIGHBORS 20
#define MAX_BUFFER_SIZE 100

using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

std::vector<int> neighboringDomainIDs(void *size, void *dg, void *global)
{
    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    if (!size) {
        return std::vector<int>();
    }
    
    //FNAME(get_neighbors_fort)(&size,
    //    &dg,
    //	&global,
    //neighbors_fort,
    //&numneighbors_fort);
//std::vector<int> ret(numneighbors_fort);

//std::copy(neighbors_fort, neighbors_fort + numneighbors_fort, ret.begin());

    std::vector<int> ret { 1,2,3 };
    return ret;

}


class DomainReference
{
public:
    class FortranPointerWrapper
    {
    public:
	explicit
	FortranPointerWrapper(void *size, void *global, void *dg, void *nodalattr) :
	    size(size),
	    global(global),
	    dg(dg),
	    nodalattr(nodalattr)
	{
	    std::cout << "Calling FortranPointerWrapper Constructor" << std::endl;
	}


	~FortranPointerWrapper()
	{
	    if (size) {
		//	      std::cout << "CPP: about to call term_fort" << std::endl;
		//		FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    }
	}
	
	void *size;
	void *global;
	void *dg;
	void *nodalattr;
           
    };

    explicit
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0, bool busywork = false) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighboringDomainIDs(size,dg,global)),
        id(id),
        timestep(0),
	rkstep(1),
	update_step(true),
	exchange_step(false),
	advance_step(false),
	busywork(busywork)
  {
      std::cout << "Calling Domain Reference Constructor" << std::endl;
  }

    hpx::future<std::map<int,std::vector<double> > > update(std::vector<std::map<int,std::vector<double> > >)
    {
	//Put futures from 
	//Loop over neighbors
	for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
	    int neighbor_here = neighbors_here[neighbor];
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];
	    
	    // Unpack buffer from neighbor
	    //const std::vector<double>*  buffer_vector = &hood[neighbor_here].output_buffer.at(id);
	    //std::vector<double> buffer_vector = &hood[neighbor_here].output_buffer.at(id);
	    //FIX		    std::vector<double> buffer_vector = hood[neighbor_here].output_buffer.at(id);
	    //std::vector<double> buffer_vector = output_buffer.at(neighbor_here);
	    	    
	    /* FIX
	       // Unpack Buffer Vector
	       for (int i=0; i<buffer_vector.size(); i++) {
	       buffer[i] = buffer_vector[i];
	       //std::cout << "buffer[" << i << "] = " << buffer[i] << " "; 
	       }
	    */
	    
	    // Put elements into our own subdomain
	    /*
	      FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
	      &neighbor_here,
	      &volume,
	      buffer);
	    */
	    
	    // Get outgoing boundarys from the neighbors	
	    //		    std::cout << "CPP: about to call hpx_get_elems_fort" << std::endl;		   
	    /*
	      FNAME(hpx_get_elems_fort)(&hood[neighbor_here].domainWrapper->dg,
	      &id,
	      &volume,
	      buffer);
	    */
	    
	    // Put those arrays inside current domain
	    //		    std::cout << "CPP: about to call hpx_put_elems_fort" << std::endl;
	}// end loop over neighbors
	
	
	// Clear output buffer map
	output_buffer.clear();
	
	// Hydro timestep
	/*
	  FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
	  &domainWrapper->dg,
	  &domainWrapper->global,
	  &domainWrapper->nodalattr,
	  &timestep,
	  &rkstep
	  );
	*/
	
	
	// Fill output buffer
	//Loop over neighbor domains
	for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
	    int neighbor_here = neighbors_here[neighbor];
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];
	    /*
	      FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
	      &neighbor_here, // pointer to neighbor to send to
	      &volume,
	      buffer);
	    */
	    
	    // Pack buffer into std vector
	    /* FIX
	       std::vector<double> buffer_vector;
	       for (int i=0; i<volume; i++) {
	       //std::cout << "buffer[" << i << "] = " << buffer[i] << " ";
	       buffer_vector.push_back(buffer[i]);		   
	       }
	    */
	    
	    // Insert pair into map
	    /* FIX
	       std::pair<std::_Rb_tree_iterator<std::pair<const int, std::vector<double> > >, bool> outval = output_buffer.insert(std::map<int, std::vector<double> >::value_type(neighbor_here, buffer_vector));
	       if (std::get<1>(outval)) {
	       //std::cout << "Insert successful!" << std::endl;
	       } else {
	       std::cout << "Insert not successful!" << std::endl;
	       }
	    */
	    
		}
	
	
	
	//std::cout << "output_buffer.size() = " << output_buffer.size() << std::endl;
	
	if (rkstep == 2) {	
	    //std::cout << "advancing domain " << id << " at timestep " << timestep <<std::endl;
	    //	      std::cout << "CPP: about to call dg_timestep_advance_fort" << std::endl;
	    /*
	      FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
	      &domainWrapper->dg,
	      &domainWrapper->global,
	      &domainWrapper->nodalattr,
	      &timestep
	      );
	    */
	    ++timestep;
	    
	}
	
    }

private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    //FortranPointerWrapper domainWrapper;
    std::vector<int> neighbors_here;
    std::map<int,std::vector<double> > output_buffer;
    int id;
    int timestep;
    int rkstep;
    bool update_step;
    bool exchange_step;
    bool advance_step;
    bool busywork;
};
  
struct stepper
{

    // Partition type
    typedef hpx::shared_future<DomainReference> partition;
    
    // Data for one time step
    typedef std::vector<DomainReference> space;

    hpx::future<space> do_work(int total_rksteps, int n_domains)
    {
	using hpx::dataflow;
	using hpx::util::unwrapped; // What is this?

	// U[t][i] is state of domain i at time t
	std::vector<space> U(2);
	//	for (space& s : U) // Is this basically a for each statement?
	//	    s.resize(n_domains);

	// Initialize Domains
	for (int i=1; i<n_domains; i++) {
	    std::cout << "about to initialize domain i=" << i << std::endl;
	    DomainReference dr_here = DomainReference(i);
	    U[0].push_back(dr_here);
	    U[1].push_back(dr_here);
	}

	auto Op = unwrapped(&DomainReference::update);
	
	// Timestep loop
	for (int substep = 0; substep != total_rksteps; ++substep) {
	    space const& current = U[substep % 2];
	    space& next = U[(substep + 1) % 2];
	    
	    // Domain loop
	    for (int i=1; i<n_domains; i++) {

		//pack futures into an array
		//loop over neighbors  (current[i] is a DomainReference)
		current[i].update(

		next[i] = dataflow(hpx::launch::async, Op,
				   
				   );
	    }
	return hpx::when_all(U[1]);
    }

};


int hpx_main(variables_map & vm)
{

    // number of domains read from config file
    int n_domains;
    //    FNAME(hpx_read_n_domains)(&n_domains);
    n_domains = 10;

    int n_timesteps = vm["n_timesteps"].as<std::size_t>();
    bool busywork = vm["busywork"].as<bool>();

    int n_rksteps = 2;
    int total_rksteps = n_timesteps*(n_rksteps*2+1)+1;
   
    //Timestepping loop goes here?
    stepper step;
    //hpx::future<stepper::space> result = step.do_work(total_rksteps, n_domains);
    int result = step.do_work(total_rksteps, n_domains);
    /*
    stepper::space solution = result.get();
    hpx::wait_all(solution);
    */
    return hpx::finalize();
}

int main(int argc, char **argv)
{
    
  // Parse command line options
  options_description
    desc_commandline("usage: " HPX_APPLICATION_STRING " [options]");

  desc_commandline.add_options()
    (
            "n_timesteps"
	    , value<std::size_t>()->default_value(100)
	    , "Number of timesteps"
     )
    (
            "busywork"
	    , value<bool>()->default_value(false)
	    , "Option to use busywork instead of calling fortran subroutines"
     )
    ;

    // We want HPX to run hpx_main() on all localities to avoid the
    // initial overhead caused by broadcasting the work from one to
    // all other localities:
    std::vector<std::string> config(1, "hpx.run_hpx_main!=1");


    // fixme: I'm not sure that config is being used by hpx_main
    return hpx::init(desc_commandline,argc,argv,config);

}
