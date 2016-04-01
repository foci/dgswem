#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cassert>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>

#include <boost/format.hpp>
#include <boost/cstdint.hpp>

#include "fname.h"
#include "fortran_declarations.hpp" // This includes parameters defined

using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

typedef std::map<int,std::vector<double> > buffer_map;

boost::uint64_t t;
//std::ofstream fortran_calls

std::vector<int> neighboringDomainIDs(void *size, void *dg, void *global)
{
    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    
    //hpx::cout << "calling neighboringDomainIDs" << std::endl;
    
    if (!size) {
#ifdef VERBOSE
      hpx::cout << "no size, so returning empty vector" << std::endl; // DEBUG
#endif
      return std::vector<int>();
    }
    
#ifdef VERBOSE
    hpx::cout << "Calling get_neighbors_fort" << std::endl;
#endif
    FNAME(get_neighbors_fort)(&size,
			      &dg,
			      &global,
			      neighbors_fort,
			      &numneighbors_fort);
    std::vector<int> ret(numneighbors_fort);

    std::copy(neighbors_fort, neighbors_fort + numneighbors_fort, ret.begin());

    return ret;

} // end neighboringDomainIDs


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
#ifdef VERBOSE
	  hpx::cout << "Calling FortranPointerWrapper Constructor" << std::endl; // DEBUG
#endif
	}


	~FortranPointerWrapper()
	{
	    if (size) {
#ifdef VERBOSE
		hpx::cout << "CPP: about to call term_fort" << std::endl; // DEBUG
#endif
		FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    }
	}
	
	void *size;
	void *global;
	void *dg;
	void *nodalattr;
	
    }; //end FortranPointerWrapper
    
    
    explicit
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighboringDomainIDs(size,dg,global)),
        id(id),
        timestep(0),
	rkstep(1)
    {
#ifdef VERBOSE
	hpx::cout << "Calling DomainReference Constructor" << std::endl; // DEBUG
#endif
    }

    
    ~DomainReference () {
#ifdef VERBOSE
	hpx::cout << "Calling DomainReference destructor" << std::endl;
#endif
    }
    
    buffer_map update(std::vector<hpx::shared_future<buffer_map> > incoming_maps_futures)
    {

	// Retrieve state variables (timestep, rkstep) in fortran data structure
#ifdef VERBOSE
	std::cout << "About to call cpp_vars_from_fort" << std::endl;
#endif

	FNAME(cpp_vars_from_fort)(&domainWrapper->size,&rkstep,&timestep);
#ifdef VERBOSE
	std::cout << "after calling cpp_vars_from_fort, timestep = " << timestep << " rkstep = " << rkstep << std::endl << std::flush;
	std::cout << "Update step, timestep = " << timestep << " rkstep = " << rkstep << " domain = " << id << " size = " << &domainWrapper->size << std::endl << std::flush;
#endif	

	if (timestep == 0) {
	    timestep++;
	    
	    FNAME(cpp_vars_to_fort)(&domainWrapper->size,&rkstep,&timestep);		
	    
	    return output_buffer;
	} // end if timestep == 0
	

#ifdef VERBOSE
	//hpx::cout << "update call, timestep = " << timestep << " rkstep = " << rkstep << " domain = " << id << std::endl << std::flush; // DEBUG
	//hpx::cout << "incoming_maps.size() = " << incoming_maps.size() << std::endl << std::flush; // DEBUG
	//hpx::cout << "neighbors_here.size() = " << neighbors_here.size() << std::endl << std::flush; // DEBUG
#endif
	
	//#############  Place ghost zones from neighboring cells into this domain  ################
	// Check size of incoming buffer, should be equal to number of neighbors
	if ( (timestep !=1 && rkstep == 1) || (rkstep == 2) )
	    {
		if (incoming_maps_futures.size() == neighbors_here.size())
		    {
			//hpx::cout << "trying to put ghost zones in domain " << id << std::endl << std::flush; // DEBUG
			//Loop over neighbors
			for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
			    int neighbor_here = neighbors_here[neighbor];
			    int volume;
			    double buffer[MAX_BUFFER_SIZE];
			    
			    // put rk indices
			    int rkindex;
			    if (rkstep == 1) rkindex = 1;
			    if (rkstep == 2) rkindex = 2;
			    
			    
			    // Unpack buffer from neighbor
			    //hpx::cout << "Unpacking buffer from neighbor: " << neighbor_here << " incoming_maps.size() = " << incoming_maps.size() << std::endl << std::flush;
			    buffer_map buffer_map_here = incoming_maps_futures[neighbor].get(); // This needs to be "neighbor" and NOT "neighbor_here"
			    // because this is a vector not a map!
			    
			    //hpx::cout << "pair from buffer: id =" << id << " buffer_map_here.size() = " << buffer_map_here.size() << std::endl << std::flush;
			    if (buffer_map_here.size() != 0) {  // Make sure the buffer isn't empty (like it always will be at the first timestep)
				std::vector<double> buffer_vector = buffer_map_here.at(id);
				// Unpack Buffer Vector
				for (int i=0; i<buffer_vector.size(); i++) {
				    buffer[i] = buffer_vector[i];
				    // hpx::cout << "buffer[" << i << "] = " << buffer[i] << " "; // DEBUG
				}
				
				// Put elements into our own subdomain
				/*
				  fortran_calls << "calling hpx_put_elems_fort, timestep = " << timestep << " rkstep = "
				  << rkstep << " domain = "<< id << " neighbor = " << neighbor_here << " rkindex = " << rkindex << " " ;
				  
				  fortran_calls << "buffer = ";
				  for (int i=0; i<MAX_BUFFER_SIZE; i++) {
				  fortran_calls << buffer[i] << " ";
				  }
				  fortran_calls << std::endl << std::flush;
				
				  fortran_calls << "buffer_vector.size() = " << buffer_vector.size() << std::endl << std::flush;
				*/
   				FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
							  &neighbor_here,
							  &volume,
							  buffer,
							  &rkindex);
			    } else {
				hpx::cout << "Buffer map is empty!" << std::endl << std::flush;
			    }
			}// end loop over neighbors
		    }
		else {
		    hpx::cout << "size of incoming buffer not equal to number of neighbors";
		    hpx::cout << " skipping placing boundary zones " << std::endl << std::flush;
		}
	    } //end if timestep !=1
	//  #########################################################################################
	
	
	
	// ######################### Hydro timestep ########################
	//	fortran_calls << "calling dg_hydro_timestep_fort, timestep = " << timestep << " rkstep = "
	//  << rkstep << " domain = "<< id << std::endl << std::flush;
	// ASSERT that rkstep != 0
	assert(rkstep != 0);
#ifdef VERBOSE
	std::cout << "calling hydro timestep, rkstep = " << rkstep << " timestep = " << timestep << std::endl << std::flush;
#endif
	FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
				      &domainWrapper->dg,
				      &domainWrapper->global,
				      &domainWrapper->nodalattr,
				      &timestep,
				      &rkstep
				      );
	// #################################################################
	
	
	
	// ###################### Pack ghost zones to transfer to neigbors ####################
	// Clear output buffer map
	output_buffer.clear();
	
	// Fill output buffer
	//Loop over neighbor domains
	for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
	    int neighbor_here = neighbors_here[neighbor];
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];
	    
	    // get rk indices
	    int rkindex;
	    if (rkstep == 1) rkindex = 2;
	    if (rkstep == 2) rkindex = 3;
	
	    //	    fortran_calls << "calling hpx_get_elems_fort, timestep = " << timestep << " rkstep = "
	    //		  << rkstep << " domain = "<< id << " neighbor = " << neighbor_here << " rkindex = " << rkindex << " ";    


	    //hpx::cout << "cpp: rkindex = " << rkindex << std::endl << std::flush; // DEBUG
	    FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
				      &neighbor_here, // pointer to neighbor to send to
				      &volume,
				      buffer,
				      &rkindex);
	    //hpx::cout << "volume = " << volume << std::endl << std::flush; // DEBUG
	     //DEBUG	     
	    //fortran_calls << "volume = " << volume << std::endl << std::flush; // DEBUG
	     /*
	    fortran_calls << "buffer = ";
	    for (int i=0; i<volume; i++) {
		fortran_calls << buffer[i] << " ";
	    }
	    fortran_calls << std::endl << std::flush;
	     */
	    
	    // Pack buffer into std vector
	    
	    std::vector<double> buffer_vector;
	    for (int i=0; i<volume; i++) {
		//hpx::cout << "buffer[" << i << "] = " << buffer[i] << " "; //DEBUG
		buffer_vector.push_back(buffer[i]);
	    }
	    
	    // hpx::cout << "after calling hpx_get_elems: buffer_vector.size() = " << buffer_vector.size() << std::endl << std::flush; // DEBUG
	    
	    //hpx::cout << "pair inserted into buffer:  neighbor_here = " << neighbor_here << " buffer_vector.size() = " << buffer_vector.size() << std::endl << std::flush; //DEBUG

	    // Insert pair into map
	    std::pair<std::_Rb_tree_iterator<std::pair<const int, std::vector<double> > >, bool> outval = output_buffer.insert(std::map<int, std::vector<double> >::value_type(neighbor_here, buffer_vector));
	    if (std::get<1>(outval)) {
		// hpx::cout << "Insert successful!" << std::endl << std::flush; // DEBUG
	    } else {
		hpx::cout << "Insert not successful!" << std::endl << std::flush; // DEBUG - EXCEPTION
	    }
	    
	} // end loop over neighbors
	// ######################################################################################
	
	
	
	//hpx::cout << "output_buffer.size() = " << output_buffer.size() << std::endl << std::flush;
	
	
	//##################### Advance domain at end of 2nd RK step ############################
	
	if (rkstep == 2) {	
	    //hpx::cout << "advancing domain " << id << " at timestep " << timestep <<std::endl << std::flush;
	    //	      hpx::cout << "CPP: about to call dg_timestep_advance_fort" << std::endl << std::flush;
	    
	    //fortran_calls << "calling dg_timestep_advance_fort, timestep = " << timestep << " domain = "<< id << std::endl << std::flush;			
	    FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
					    &domainWrapper->dg,
					    &domainWrapper->global,
					    &domainWrapper->nodalattr,
					    &timestep
					    );
	    
	} // end if rkstep == 2	
	
	
	if (rkstep == 2) {		    
	    ++timestep;
	    --rkstep;
	} else {
	    ++rkstep;
	} // end if else rkstep == 2
	
	//########################################################################################

	FNAME(cpp_vars_to_fort)(&domainWrapper->size,&rkstep,&timestep);
	
	// Return output buffer
	//	return hpx::make_ready_future(output_buffer);
	return output_buffer;
	
    } // End update()
    
    
    std::vector<int> neighbors_here;
    buffer_map output_buffer;
    
private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    int id;
    int timestep;
    int rkstep;
}; // End class DomainReference

struct stepper
{
    // Domains
    std::vector<DomainReference> domains;

    // Partition type
    typedef hpx::shared_future<buffer_map> partition; //uncomment me to futurize
    
    // Data for one time step
    typedef std::vector<partition> space; //uncomment me to futurize

    
    //    typedef std::vector<buffer_map> space; //comment me out to futurize

    //void do_work(int total_substeps, int n_domains) // comment me out to futurize
    //    void do_work(int total_substeps, int n_domains)
    hpx::future<space> do_work(int total_substeps, int n_domains) // uncomment me to futurize
    {
	using hpx::dataflow;
	//using hpx::util::unwrapped; 

	// U[t][i] is state of domain i at time t
	std::vector<space> U(2);
	for (space& s : U)
	    s.resize(n_domains);

	//domains.resize(n_domains); // resize to number of domains

	// Initialize Domains
	for (int i=0; i<n_domains; i++) {
#ifdef VERBOSE
	  hpx::cout << "about to initialize domain i=" << i << std::endl << std::flush;
#endif
      	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    //fortran_calls << "calling dgswem_init_fort, domain = "<< i <<std::endl << std::flush;	    
	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &i);

	    /* Don't need to do this
	    // Get list of neighbors for this domain
	    std::vector<int> neighbors_here = neighboringDomainIDs(size, dg, global);      

	    hpx::cout << "neighbors_here.size() = " << neighbors_here.size() << std::endl << std::flush;
	    if (neighbors_here.size() == 0) {
	    }
	    */
	    DomainReference dr_here = DomainReference(i, size, global, dg, nodalattr);
	    domains.push_back(dr_here);
	    //domains[i] = dr_here;

	    buffer_map empty_buffer;
	    //U[0].push_back(empty_buffer); // comment me out to futurize
	    //U[0].push_back(hpx::make_ready_future(empty_buffer)); // uncomment me to futurize
	    U[0][i]=hpx::make_ready_future(empty_buffer); // uncomment me to futurize
	} //end domain initialize

#ifdef VERBOSE
	hpx::cout << "***** Done initializing domains ******" << std::endl << std::flush;
#endif

	// Start time variable
	t = hpx::util::high_resolution_clock::now();    
	// Reset counters
	hpx::reset_active_counters();

	// Timestep loop
	for (int substep = 0; substep != total_substeps; ++substep) {
	    space const& current = U[substep % 2];
	    space& next = U[(substep + 1) % 2];
#ifdef VERBOSE	    
	    hpx::cout << "Timestep loop: substep = " << substep << std::endl << std::flush;
#endif
	    //fortran_calls << "substep = " << substep << std::endl << std::flush;
	    // Domain loop
	    for (int i=0; i<n_domains; i++) {
		//hpx::cout << "Domain loop: substep = " << substep << " domain = " << i << std::endl << std::flush;
		//hpx::cout << "current[i].size() = " << current[i].size();
		
		//pack futures into an array
		std::vector<hpx::shared_future<buffer_map> > output_buffer_future_vector; //uncomment me to futurize
		//std::vector<buffer_map> output_buffer_future_vector; //comment me out to futurize
		//loop over neighbors current[ ] is a future of a buffer_map
		for (int neighbor=0; neighbor<domains[i].neighbors_here.size(); neighbor++) {
		    int neighbor_here = domains[i].neighbors_here[neighbor];
		    //hpx::cout << "packing domains, domain = " << i << " neighbor = " << neighbor << " current[neighbor_here].size() = " << current[neighbor_here].size() << std::endl << std::flush; // DEBUG
		    output_buffer_future_vector.push_back(current[neighbor_here]);
		    //output_buffer_future_vector.push_back(current[neighbor_here]);
		}
		//output_buffer_future_vector.push_back(current[i]);

		//Define Op

		// Passing by reference:
#ifdef VERBOSE
		hpx::cout << "&domains[" << i << "] = " << &domains[i] << std::endl << std::flush;
#endif
		auto Op = boost::bind(&DomainReference::update,&domains[i],_1);

		// Passing by copying:    and should this be outside of the loops? I don't think it can be 
		//auto Op = boost::bind(&DomainReference::update,domains[i],_1);
		
		//hpx::cout << "Calling dataflow" << std::endl << std::flush;
		next[i] = dataflow(hpx::launch::async, Op,
				   output_buffer_future_vector,
				   current[i]
				   );
	
		//hpx::wait_all(next);
	    }
	} //end timestep loop
#ifdef VERBOSE	
	hpx::cout << " &&&&&&& Done building tree &&&&&& " << std::endl << std::flush;
#endif
	return hpx::when_all(U[total_substeps % 2]); // uncomment to futurize
	//	hpx::wait_all(U[total_substeps % 2]); // uncomment to futurize
	//hpx::wait_all(U[total_substeps]); // uncomment to futurize
	
    }// end do_work()

};


int hpx_main(variables_map & vm)
{

    // number of domains read from config file
    int n_domains;
    FNAME(hpx_read_n_domains)(&n_domains);

    int n_timesteps = vm["n_timesteps"].as<std::size_t>();
    bool busywork = vm["busywork"].as<bool>();

    int n_rksteps = 2;
    int total_substeps = n_timesteps*(n_rksteps)+1; 
   
    stepper step;


    // Uncomment to futurize
    hpx::future<stepper::space> result = step.do_work(total_substeps, n_domains);

    //stepper::space result = step.do_work(total_substeps, n_domains);
    stepper::space solution = result.get();
    hpx::wait_all(solution);
    
    boost::uint64_t elapsed = hpx::util::high_resolution_clock::now() - t;    

    boost::uint64_t const num_os_threads = hpx::get_os_thread_count();

    // Close debug output file
    //   fortran_calls.close();

    std::string const threads_str = boost::str(boost::format("%lu,") % num_os_threads);
    hpx::cout << ( boost::format("%-21s %.14g\n")
		   % threads_str % (elapsed / 1e9) ) << std::flush;

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
  //std::vector<std::string> config(1, "hpx.run_hpx_main!=1");


  //return hpx::init(desc_commandline,argc,argv,config);
    return hpx::init(desc_commandline,argc,argv);

}
