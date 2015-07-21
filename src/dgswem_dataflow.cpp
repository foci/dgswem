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

typedef std::map<int, double[MAX_BUFFER_SIZE]>  sendbuffers;
typedef std::map<int, double[MAX_BUFFER_SIZE]>  recvbuffers;

struct dgswem_domain
{
  void *size = NULL;
  void *dg = NULL;
  void *global = NULL;
  void *nodalattr = NULL;
  int id = 0;
  int numneighbors;
  std::vector<int> neighbors;
  
  //I would like to get rid of these
  int n_timesteps;
  int n_domains;
  int n_rksteps;

  dgswem_domain(int i) // Constructor initializes the domain
  {
    id = i;
    FNAME(dgswem_init_fort)(&size,
			    &dg,
			    &global,
			    &nodalattr,
			    &n_timesteps,
			    &n_domains,
			    &id,
			    &n_rksteps
			    );
    //Get a list of neighbors for each domain

    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    //    std::vector<int> neighbors_here;
    FNAME(get_neighbors_fort)(&size,
			      &dg,
			      &global,
			      neighbors_fort,
			      &numneighbors_fort);

    numneighbors = (numneighbors_fort);    
    for (int j=0; j<numneighbors_fort; j++) {
      neighbors.push_back(neighbors_fort[j]);
    }
  };

  sendbuffers update(int timestep, int rkstep, recvbuffers)
  {
    FNAME(dg_timestep_fort)(&size,
			    &dg,
			    &global,
			    &nodalattr,
			    &timestep,
			    &rkstep
			    );
    sendbuffers sendbuffers_here; 
    
    //Loop over neighbors
    for (int neighbor=0; neighbor<numneighbors; neighbor++) {	
      int neighbor_here = neighbors[neighbor];
      int volume;
      double buffer[MAX_BUFFER_SIZE];
      
      // Get outgoing boundarys from the neighbors
      FNAME(hpx_get_elems_fort)(&dgs[neighbor_here],
				&domain,
				&volume,
				buffer);
      sendbuffers.insert(neighbor_here,buffer);      
    };


  }
  
};
/*
struct stepper 
{

  typedef hpx::shared_future<domain> domain_future;

  typedef std::vector<domain_future> space;
*/
  

int hpx_main(
	 int argc
	 , char* argv[]
	 )
{

  int n_timesteps;
  int n_domains;
  int n_rksteps;

  
  std::vector<dgswem_domain> domains;

  FNAME(hpx_read_n_domains)(&n_domains);

  std::cout << "n_timesteps = " << n_timesteps << std::endl;

  for (int i=0; i<n_domains; i++) {
    dgswem_domain domain(i); // This initializes the domain
    domains.push_back(domain);
  }

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
  
  /*  
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
  
  */
  
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
