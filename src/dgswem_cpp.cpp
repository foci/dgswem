#include <iostream>
#include <fstream>
#include "fname.h"
#include <vector>
#include <cstdlib>

#include "fortran_declarations.hpp" // This includes parameters defined

int main(
	 int argc
	 , char* argv[]
	 )
{

    int n_timesteps = 100;  

    std::vector <std::string> args;
    int n_timesteps_cmdline;
    
    std::ofstream fortran_calls;
    fortran_calls.open("fortran_calls.txt");


    for (int i=1; i<argc; ++i) {
	std::string arg = argv[i];
	if((arg == "-h") || (arg == "--help")) {
	    std::cerr << "Usage: " << argv[0] 
		      << " --n_timesteps <number of timesteps>" << std::endl;
	    return 1;
	} else if  ((arg == "-n") || (arg == "--n_timesteps")) {		
	    if (i+1 < argc) { // make sure there's one more argument
		
		//std::cout << "arg is " << std::atoi(argv[i+1]) << std::endl;
		const char *str = argv[++i];
		//std::cout << "argv[i++] = " << argv[i++] << std::endl;
		//std::cout << "str = " << str << std::endl;
		n_timesteps_cmdline = std::atoi(str);
		n_timesteps = n_timesteps_cmdline;
	    } else { // no argument
		std::cerr << "--n_timesteps requires one argument" << std::endl;
		return 1;
	    } 
	} else {
	    std::cerr << "Unknown argument" << argv[i] <<std::endl;
	}
    }


  std::vector<void *> sizes;
  std::vector<void *> dgs;
  std::vector<void *> globals;
  std::vector<void *> nodalattrs;

  std::vector<int> ids;

  std::vector<int> numneighbors;
  std::vector<std::vector<int> > neighbors;

  int n_domains;
  int n_rksteps = 2;
  

  FNAME(hpx_read_n_domains)(&n_domains);

  std::cout << "c++: n_domains = " << n_domains << std::endl;

  // Create vectors of fortran pointers and ids
  for (int i=0; i<n_domains; i++) {
    void *size = NULL;
    void *dg = NULL;
    void *global = NULL;
    void *nodalattr = NULL;
    ids.push_back(i);
    sizes.push_back(size);
    dgs.push_back(dg);
    globals.push_back(global);
    nodalattrs.push_back(nodalattr);
  }


 
  std::cout << "n_domains = " << n_domains << std::endl;

  // Initialize all domains
  for(int i=0; i<n_domains; i++) {
    std::cout << "initializing domain " << i << std::endl;
    fortran_calls << "calling dgswem_init_fort, domain = "<< i <<std::endl;	    
    FNAME(dgswem_init_fort)(&sizes[i],
			    &dgs[i],
			    &globals[i],
			    &nodalattrs[i],
			    &ids[i]
			    );
    
    //Get a list of neighbors for each domain

    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    std::vector<int> neighbors_here;
    fortran_calls << "Calling get_neighbors_fort" << std::endl;
    FNAME(get_neighbors_fort)(&sizes[i],
			      &dgs[i],
			      &globals[i],
			      neighbors_fort,
			      &numneighbors_fort);

    numneighbors.push_back(numneighbors_fort);
    
    for (int j=0; j<numneighbors_fort; j++) {
      std::cout << neighbors_fort[j] << " ";
      neighbors_here.push_back(neighbors_fort[j]);
    }
    std::cout << std::endl;

    neighbors.push_back(neighbors_here);

    std::cout << "Done initializing " << i << std::endl;
    std::cout << "n_domains = " << n_domains << std::endl;
  } // End loop over domains


  // Print out some information about the domains and their neighbors
  std::cout << "*** Grid Information ***" << std::endl;
  std::cout << "ids.size() = " << ids.size() << std::endl;
  std::cout << "numneighbors.size() = " << numneighbors.size() << std::endl;
  for (int domain=0; domain<numneighbors.size(); domain++) {
    std::cout << "numneighbors[" << domain << "] = " << numneighbors[domain] << std::endl;
    std::cout << "neighbors: ";
    std::vector<int> neighbors_here = neighbors[domain];
    for (int neighbor=0; neighbor<numneighbors[domain]; neighbor++) {      
      std::cout << neighbors_here[neighbor] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "*** End Grid Information ***" << std::endl;

  std::cout << "c++: Starting timestep loop: n_timesteps = " << n_timesteps << " n_rksteps = " << n_rksteps << std::endl;
  
  
  // Start timestepping loop
  for (int timestep=1; timestep<=n_timesteps; timestep++) {
    
    std::cout << "starting timestep loop, timestep = " << timestep << std::endl;

    for (int rkstep=1; rkstep<=n_rksteps; rkstep++) {
/*     std::cout << "#################################################################" << std::endl;
     std::cout << " timestep loop, timestep = " << timestep << std::endl;     
     std::cout << " rk loop, rkstep = " << rkstep << std::endl;
     std::cout << "#################################################################" << std::endl;   */  
      
      for (int j=0; j<ids.size(); j++) {
	//std::cout << "j=" << j << std::endl;
	        std::cout << "updating (domain_id = " << ids[j]
		  << ", timestep = " << timestep
		  << ", rkstep = " << rkstep
		  << ")...\n";
		fortran_calls << "calling dg_hydro_timestep_fort, timestep = " << timestep << " rkstep = "
			      << rkstep << " domain = "<< j << std::endl;
		FNAME(dg_hydro_timestep_fort)(&sizes[j],
				&dgs[j],
				&globals[j],
				&nodalattrs[j],
				&timestep,
				&rkstep
				);
      } // End loop over domains
      

      // Element Boundary exchange
      // Loop over domains   
      for (int domain=0; domain<ids.size(); domain++) {
	std::vector<int> neighbors_here = neighbors[domain];

	//Loop over neighbors
	for (int neighbor=0; neighbor<numneighbors[domain]; neighbor++) {	
	  int neighbor_here = neighbors_here[neighbor];
	  int volume;
	  double buffer[MAX_BUFFER_SIZE];

	  std::cout << "domain " << domain << " is exchanging with " << neighbor_here << " at timestep " << timestep << std::endl;
	  
	  // Get outgoing boundarys from the neighbors
	  fortran_calls << "calling hpx_get_elems_fort, timestep = " << timestep << " rkstep = "
			<< rkstep << " domain = "<< domain << " neighbor = " << neighbor_here << std::endl;
	  int rkindex = 0;
	  FNAME(hpx_get_elems_fort)(&dgs[neighbor_here],
				    &domain,
				    &volume,
				    buffer,
				    &rkindex);
	  // DEBUG
	  //	  fortran_calls << "volume = " << volume << std::endl; // DEBUG
	  /*
	  fortran_calls << "buffer = ";
	  for (int i=0; i<volume; i++) {
	      fortran_calls << buffer[i] << " ";
	  }
	  fortran_calls << std::endl;
	  */
	  // Put those arrays inside current domain 
	  fortran_calls << "calling hpx_put_elems_fort, timestep = " << timestep << " rkstep = "
			<< rkstep << " domain = "<< domain << " neighbor = " << neighbor_here << std::endl;	  
	  /*
	  fortran_calls << "buffer = ";
	  for (int i=0; i<MAX_BUFFER_SIZE; i++) {
	      fortran_calls << buffer[i] << " ";
	  }
	  fortran_calls << std::endl;	  
	  */ 	  //fortran_calls << "buffer_vector.size() = " << buffer_vector.size() << std::endl;
	  FNAME(hpx_put_elems_fort)(&dgs[domain],
				    &neighbor_here,
				    &volume,
				    buffer,
				    &rkindex);	
	  
	  
	}// end loop over neighbors
	
      }// end loop over domains

      //end boundary exchange


      // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Slopelimiter_part1 &&&&&&&&&&&&&&&&&
      for (int domain=0; domain<ids.size(); domain++) {
	  fortran_calls << "calling SL part 1, timestep = " << timestep << " rkstep = " << rkstep 
			<< " domain = " << domain << std::endl;
      }
      // &&&&&&&&&&&&&&&&&&&&&&& Slopelimiter boundary exchange &&&&&&&&&&&&&&
      for (int domain=0; domain<ids.size(); domain++) {
      fortran_calls << "calling SL boundary exchange GET, timestep = " << timestep << " rkstep = " << rkstep 
		   << " domain = " << domain << std::endl;
      }

      for (int domain=0; domain<ids.size(); domain++) {
      fortran_calls << "calling SL boundary exchange PUT, timestep = " << timestep << " rkstep = " << rkstep 
		   << " domain = " << domain << std::endl;
      }
      // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Slopelimiter_part2 &&&&&&&&&&&&&&&&&
      fortran_calls << "calling SL part 2, timestep = " << timestep << " rkstep = " << rkstep 
		   << " domain = " << domain << std::endl;


      // &&&&&&&&&&&&&&&&&&&&&&&& wetting and drying &&&&&&&&&&&&&&&&&&&&&&&&&
      for (int j=0; j<ids.size(); j++) {
	//std::cout << "j=" << j << std::endl;
	        std::cout << "calling wetdry() (domain_id = " << ids[j]
		  << ", timestep = " << timestep
		  << ", rkstep = " << rkstep
		  << ")...\n";
		fortran_calls << "calling dg_wetdry_timestep_fort, timestep = " << timestep << " rkstep = "
			      << rkstep << " domain = "<< j << std::endl;
		FNAME(wetdry_fort)(&dgs[j],
				   &globals[j]);
      } // End loop over domains

      // &&&&&&&&&&&&&&&&&&&&& end wetting and drying &&&&&&&&&&&&&&&&&&&&&&&
       
      }// end loop over domains
      //return 0;
      
    } // end rkstep loop
    
    for (int domain=0; domain<ids.size(); domain++) {
	std::cout << "advancing domain " << domain << " at timestep " << timestep <<std::endl;  
	fortran_calls << "calling dg_timestep_advance_fort, timestep = " << timestep << " domain = "<< domain << std::endl;
	FNAME(dg_timestep_advance_fort)(&sizes[domain],
				      &dgs[domain],
				      &globals[domain],
				      &nodalattrs[domain],
				      &timestep
				      );
    }

     if ( timestep > 2) {
      return 0;  // stop after one timestep for debugging
     }


  } // End timestep loop


  
  
  // Destroy all domains
  //std::cout << "Destroying domains" << std::endl;
  for (int domain=0; domain<ids.size(); domain++) {
      fortran_calls << "calling term_fort" << std::endl;
      FNAME(term_fort)(&sizes[domain],
		       &dgs[domain],
		       &globals[domain],
		       &nodalattrs[domain]
		       );
  }
  
  fortran_calls.close();
  return 0;  
}
