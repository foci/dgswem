#include <iostream>
#include "fname.h"
#include <vector>

#include <hpx/hpx.hpp>
#include <hpx/hpx_start.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/parallel_for_each.hpp>

#include "fortran_declarations.hpp"

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

  int n_timesteps;
  int n_domains;
  int n_rksteps;
  

  // This is hacky and needs to be changed
  {
    int dummy_id = 0;
    void *size = NULL;
    void *dg = NULL;
    void *global = NULL;
    void *nodalattr = NULL;
    // Initialize dummy domain so we can get n_domains
    FNAME(dgswem_init_fort)(&size,
			    &dg,
			    &global,
			    &nodalattr,
			    &n_timesteps,
			    &n_domains,
			    &dummy_id,
			    &n_rksteps
			    );
  }

  std::cout << "n_timesteps = " << n_timesteps << std::endl;
  std::cout << "n_rksteps = " << n_rksteps << std::endl;
  std::cout << "n_domains = " << n_domains << std::endl;


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


 

  // Initialize all domains
  for(int i=0; i<n_domains; i++) {
#ifdef HPX
    inits.push_back(hpx::async(FNAME(dgswem_init_fort),
			       &sizes[i],
			       &dgs[i],
			       &globals[i],
			       &nodalattrs[i],
			       &n_timesteps,
			       &n_domains,
			       &ids[i],
			       &n_rksteps
			       ));
#else
    FNAME(dgswem_init_fort)(&sizes[i],
			    &dgs[i],
			    &globals[i],
			    &nodalattrs[i],
			    &n_timesteps,
			    &n_domains,
			    &ids[i],
			    &n_rksteps
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

  }

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

      std::cout << "starting rk loop, rkstep = " << rkstep << std::endl;
      
#ifdef HPX
      std::vector<hpx::future<void> > updates;    
#endif
      for (int j=0; j<ids.size(); j++) {
	std::cout << "j=" << j << std::endl;
#ifdef HPX
	updates.push_back(hpx::async(FNAME(dg_timestep_fort),
				     &sizes[j],
				     &dgs[j],
				     &globals[j],
				     &nodalattrs[j],
				     &timestep,
				     &rkstep
				     ));
#else
	FNAME(dg_timestep_fort)(&sizes[j],
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

      // Boundary exchange
      // Loop over domains   
      for (int domain=0; domain<ids.size(); domain++) {
	std::vector<int> neighbors_here = neighbors[domain];
	
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
	  
	  
	  
	}// end loop over neighbors
	
      }// end loop over domains
      
    } // end rkstep loop
    
    return 0;  // stop after one timestep for debugging


  } // End timestep loop
  
  
  return hpx::finalize();
  
} // End hpx_main


int main(
	 int argc
	 , char* argv[]
	 )
{
  // Initialize and run HPX
  
  hpx::init(argc,argv);

  return 0;
}