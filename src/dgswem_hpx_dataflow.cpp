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

    void update() //FIXME stuff goes here
    {
	
      bool only_busywork = this->busywork;
      
      if (only_busywork) {
	
	// Adding some busy work FIXME
	volatile double a;
	volatile double b;
	volatile double c;
	a = 23.2394908;
	std::cout << "timestep = " << timestep << " domainID = " << id << std::endl;
	for (int i=0; i<1000000; i++) {
	  //		for (int j=0; j<10000000; j++) {
	  b = 2.23423;
	  c = a/b;
	  c = a*(b+0.000001);
	  a = c;		  
	  //		}

	}
	timestep++;
      
      } else {
	
	//      std::cout << "CPP: LGD update" << std::endl;
	if (timestep != 0) {
	  
	  if (update_step) {

	      // Clear output buffer map
	      output_buffer.clear();
	    
	      //std::cout << "********* update step ************" << std::endl;

	      /*
	      std::cout << "updating (domain_id = " << id
			  << ", timestep = " << timestep
			  << ", rkstep = " << rkstep
			  << ")...\n";
	      */

		//                 std::cout << "CPP: about to call dg_hydro_timestep_fort" << std::endl;
	      /*
	      FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
					      &domainWrapper->dg,
					      &domainWrapper->global,
					      &domainWrapper->nodalattr,
					      &timestep,
					      &rkstep
					      );
	      */
		//std::cout << "a=" << a << std::endl;
		
		
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
		    //std::cout << "volume = " << volume << std::endl;

		    //std::cout << "******************************************************" << std::endl;

		    //std::map<int,std::vector<double> > output_buffer;
		    //std::map<int,std::vector<double> > test_map;

		    
		    //std::cout << "Buffer from " << id << ", for " << neighbor_here << std::endl;
		    //std::cout << "volume = " << volume << std::endl;
		    // Pack buffer into std vector
		    /* FIX
		    std::vector<double> buffer_vector;
		    for (int i=0; i<volume; i++) {
			//std::cout << "buffer[" << i << "] = " << buffer[i] << " ";
			buffer_vector.push_back(buffer[i]);
			
		    }
		    */
		    //std::cout << std::endl;

		    /*
		    std::cout << "Buffer Vector:" << std::endl;
		    std::cout << "buffer_vector.size() = " << buffer_vector.size() << std::endl;
		    for (int i=0; i<buffer_vector.size(); i++) {
			std::cout << "buffer_vector[" << i << "] = " << buffer_vector[i] << " ";
		    }
		    std::cout << std::endl;
		    */

		    /* FIX
		    std::pair<std::_Rb_tree_iterator<std::pair<const int, std::vector<double> > >, bool> outval = output_buffer.insert(std::map<int, std::vector<double> >::value_type(neighbor_here, buffer_vector));
		    if (std::get<1>(outval)) {
			//std::cout << "Insert successful!" << std::endl;
		    } else {
			std::cout << "Insert not successful!" << std::endl;
		    }
		    */

		    //std::cout << "Vector retrieved from map:" << std::endl;
		    //std::vector<double> temp = output_buffer.at(neighbor_here);
		    /*
		    std::cout << "temp.size() = " << temp.size() << std::endl;
		    for (int i=0; i<temp.size(); i++) {	
			std::cout << "temp[" << i << "] = " << temp[i] << " "; 
		    }
		    std::cout << std::endl;
		    std::cout << "******************************************************" << std::endl;
		    */

		    //std::cout << "output_buffer.at(neighbor_here).size() = " << output_buffer.at(neighbor_here).size() << std::endl;
		    //std::cout << "temp.size() = " << temp.size() << std::endl;
		}

		

		//std::cout << "output_buffer.size() = " << output_buffer.size() << std::endl;

		update_step = false;
		exchange_step = true;
		advance_step = false;
	    } else if (exchange_step) {

		// Boundary exchange
		//Loop over neighbors
		for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
		    int neighbor_here = neighbors_here[neighbor];
		    int volume;
		    double buffer[MAX_BUFFER_SIZE];

		    // Unpack buffer from neighbor
		    //std::cout << "id here:" << id << ", unpacking buffer from " << neighbor_here << std::endl;
		    //const std::vector<double>*  buffer_vector = &hood[neighbor_here].output_buffer.at(id);
		    //std::vector<double> buffer_vector = &hood[neighbor_here].output_buffer.at(id);
		    //FIX		    std::vector<double> buffer_vector = hood[neighbor_here].output_buffer.at(id);
		    //std::vector<double> buffer_vector = output_buffer.at(neighbor_here);
		    //for (int i=0; i<buffer_vector->size(); i++) {
		    /* FIX
		    for (int i=0; i<buffer_vector.size(); i++) {
			buffer[i] = buffer_vector[i];
			//std::cout << "buffer[" << i << "] = " << buffer[i] << " "; 
		    }
		    */
		    //std::cout << std::endl;

		    // Put elements into our own subdomain
		    /*
		    FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
					      &neighbor_here,
					      &volume,
					      buffer);
		    */
		    //std::cout << "domain " << id << " is exchanging with " << neighbor_here << " at timestep " << timestep << std::endl;

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
		
		update_step = true;
		exchange_step = false;
		advance_step = false;
	    } else {
		std::cout << "Something went wrong!" << std::endl;
	    }
	    	    
	} else {
	    ++timestep;
	}

      } // else busywork

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

    //hpx::future<space> do_work(int total_rksteps, int n_domains)
    int do_work(int total_rksteps, int n_domains)
    {
	using hpx::dataflow;
	using hpx::util::unwrapped; // What is this?

	std::vector<space> U(2);
	//	for (space& s : U) // Is this basically a for each statement?
	//	    s.resize(n_domains);

	// Initialize Domains
	for (int i=1; i<n_domains; i++) {
	    std::cout << "about to initialize domain i=" << i << std::endl;
	    DomainReference dr_here = DomainReference(i);
	    U[1].push_back(dr_here);
	    U[2].push_back(dr_here);
	}
	
	//return hpx::when_all(U[1]);
	return 1;
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
