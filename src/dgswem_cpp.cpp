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

#define MAX_DOMAIN_NEIGHBORS 10

extern"C" {
  void FNAME(dgswem_init_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* n_timesteps,
			       int* n_domains,
			       int* id);
  void FNAME(dg_timestep_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* timestep);	
  void FNAME(get_neighbors_fort)(void** sizes,
				 void** dg,
				 void** global,
				 int neighbors[MAX_DOMAIN_NEIGHBORS],
				 int* num_neighbors);				 
}


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

  n_domains = 4;

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
			       &ids[i]
			       ));
#else
    FNAME(dgswem_init_fort)(&sizes[i],
			    &dgs[i],
			    &globals[i],
			    &nodalattrs[i],
			    &n_timesteps,
			    &n_domains,
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
    
    std::cout << "domain =" << i << std::endl;
    std::cout << "num_neighbors = " << numneighbors_fort << std::endl;
    std::cout << "neighbors: ";
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


  std::cout << "c++: n_timesteps = " << n_timesteps << std::endl;

  
  // Start timestepping loop
  for (int timestep=1; timestep<=n_timesteps; timestep++) {
    
#ifdef HPX
    std::vector<hpx::future<void> > updates;    
#endif
    for (int j=0; j<ids.size(); j++) {
      #ifdef HPX
      updates.push_back(hpx::async(FNAME(dg_timestep_fort),
				   &sizes[j],
				   &dgs[j],
				   &globals[j],
				   &nodalattrs[j],
				   &timestep
				   ));
#else
      FNAME(dg_timestep_fort)(&sizes[j],
			      &dgs[j],
			      &globals[j],
			      &nodalattrs[j],
			      &timestep
			      );
#endif
    } // End loop over domains
#ifdef HPX
    wait_all(updates);
#endif

  } // End timestep loop
  
  
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