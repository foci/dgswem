#include <iostream>
#include "fname.h"
#include <vector>

#ifdef HPX
#include <hpx/hpx.hpp>
#include <hpx/hpx_start.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/parallel_for_each.hpp>
#endif

#include "fortran_declarations.hpp" // This includes parameters defined

int hpx_main(
	 int argc
	 , char* argv[]
	 )
{
  std::vector<void *> sizes;
  std::vector<void *> dgs;
  std::vector<void *> globals;
  std::vector<void *> nodalattrs;

  std::vector<int> ids;

  std::vector<int> numneighbors;
  std::vector<std::vector<int> > neighbors;
  
#ifdef HPX
  std::vector<hpx::future<void> > inits;
#endif

  //  int n_timesteps = 4000;
  int n_timesteps = 86400;
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
#ifdef HPX
    inits.push_back(hpx::async(FNAME(dgswem_init_fort),
			       &sizes[i],
			       &dgs[i],
			       &globals[i],
			       &nodalattrs[i],
			       &ids[i]
			       ));
#else
    FNAME(dgswem_init_fort)(&sizes[i],
			    &dgs[i],
			    &globals[i],
			    &nodalattrs[i],
			    &ids[i]
			    );
#endif
    
    //Get a list of neighbors for each domain

    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    std::vector<int> neighbors_here;
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

#ifdef HPX
  wait_all(inits);
#endif

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
      
#ifdef HPX
      std::vector<hpx::future<void> > updates;    
#endif
      for (int j=0; j<ids.size(); j++) {
	//std::cout << "j=" << j << std::endl;
#ifdef HPX
	updates.push_back(hpx::async(FNAME(dg_hydro_timestep_fort),
				     &sizes[j],
				     &dgs[j],
				     &globals[j],
				     &nodalattrs[j],
				     &timestep,
				     &rkstep
				     ));
#else
	std::cout << "about to call dg_hydro_timestep_fort" << std::endl;
	std::cout << "sizes["<<j<<"] = " << sizes[j] << std::endl;
	FNAME(dg_hydro_timestep_fort)(&sizes[j],
				&dgs[j],
				&globals[j],
				&nodalattrs[j],
				&timestep,
				&rkstep
				);
#endif
      } // End loop over domains
#ifdef HPX
      wait_all(updates);
#endif
      


      /*
      // Boundary exchange
      // Loop over domains   
      for (int domain=0; domain<ids.size(); domain++) {
	std::vector<int> neighbors_here = neighbors[domain];
	
        //std::cout << "#################################################################" << std::endl;
        //std::cout << " domain loop, domain = " << domain << std::endl;     
        //std::cout << "#################################################################" << std::endl;             
        
	//Loop over neighbors
	for (int neighbor=0; neighbor<numneighbors[domain]; neighbor++) {	
	  int neighbor_here = neighbors_here[neighbor];
	  int volume;
	  double buffer[MAX_BUFFER_SIZE];
	  
	  // Get outgoing boundarys from the neighbors
	  FNAME(hpx_get_elems_fort)(&dgs[neighbor_here],
				    &domain,
				    &volume,
				    buffer);
	  
	  // Put those arrays inside current domain
	  FNAME(hpx_put_elems_fort)(&dgs[domain],
				    &neighbor_here,
				    &volume,
				    buffer);	
         
//          FNAME(hpx_swap_elems_fort)(&dgs[domain],
//                                     &dgs[neighbor_here]);
	  
	  
	  
	}// end loop over neighbors
	
      }// end loop over domains
      
      */
      //return 0;
      
    } // end rkstep loop
    
    for (int domain=0; domain<ids.size(); domain++) {
      FNAME(dg_timestep_advance_fort)(&sizes[domain],
				      &dgs[domain],
				      &globals[domain],
				      &nodalattrs[domain],
				      &timestep
				      );
    }

//     if ( timestep > 2) {
//      return 0;  // stop after one timestep for debugging
//     }


  } // End timestep loop


  
  
  // Destroy all domains
  std::cout << "Destroying domains" << std::endl;
  for (int domain=0; domain<ids.size(); domain++) {
      FNAME(term_fort)(&sizes[domain],
		       &dgs[domain],
		       &globals[domain],
		       &nodalattrs[domain]
		       );
  }
  
  
#ifdef HPX
  return hpx::finalize();
#else
  return 0;
#endif
  
} // End hpx_main


int main(
	 int argc
	 , char* argv[]
	 )
{
  // Initialize and run HPX
  
#ifdef HPX
  hpx::init(argc,argv);
#else  
  hpx_main(argc,argv);
#endif
  return 0;
}
