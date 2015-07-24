#include <iostream>
#include "fname.h"
#include <vector>
#include <map>
#include <cstring>

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

#include "fortran_declarations.hpp"
/*
struct buffer
{
  double this_buffer[MAX_BUFFER_SIZE];
};
*/

typedef std::map<int, std::vector<double> >  buffers;

struct dgswem_domain
{
  void *size = NULL;
  void *dg = NULL;
  void *global = NULL;
  void *nodalattr = NULL;
  int id = 0;
  int numneighbors;
  std::vector<int> neighbors;

  buffers sendbuffers;
  
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

  // update returns outgoing buffers map, indexed with domain id
  void update(int timestep, int rkstep, buffers recvbuffers)
  {

    if (timestep > 1) {
      // Insert buffers from recvbuffers
      // Loop over neighbors
      for (int neighbor=0; neighbor<numneighbors; neighbor++) {	
	int neighbor_here = neighbors[neighbor];
	int volume;
	double recvbuffer_array[MAX_BUFFER_SIZE];
	std::vector<double> recvbuffer_vector;
	
	recvbuffer_vector = recvbuffers[neighbor_here];
	
	// Put recvbuffer into c-array
	for (int i=0; i<recvbuffers[neighbor_here].size(); i++) {
	  recvbuffer_array[i]=recvbuffer_vector[i];
	};
	
	FNAME(hpx_put_elems_fort)(&dg,
				  &neighbor_here,
				  &volume,
				  recvbuffer_array);
      };
    }

    // Actual update
    FNAME(dg_timestep_fort)(&size,
			    &dg,
			    &global,
			    &nodalattr,
			    &timestep,
			    &rkstep
			    );
    buffers sendbuffers_here; 
    
    //Loop over neighbors
    for (int neighbor=0; neighbor<numneighbors; neighbor++) {	
      int neighbor_here = neighbors[neighbor];
      int volume;
      double sendbuffer_array[MAX_BUFFER_SIZE];
      std::vector<double> sendbuffer_vector;


      FNAME(hpx_get_elems_fort)(&dg,
				&neighbor_here,
				&volume,
				sendbuffer_array);

      std::cout << "volume =" << volume << std::endl;
     
      // Put sendbuffer into c-array
      for (int i=0; i<volume; i++) {
	sendbuffer_array[i]=sendbuffer_vector[i];
      };

      sendbuffers_here[neighbor_here] = sendbuffer_vector;
    };
    
    this->sendbuffers=sendbuffers_here;

  }
  
};  

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

  n_timesteps = domains[0].n_timesteps;
  n_rksteps = domains[0].n_rksteps;

  std::cout << "c++: Starting timestep loop: n_timesteps = " << n_timesteps << " n_rksteps = " << n_rksteps << std::endl;
    
  // Start timestepping loop
  for (int timestep=1; timestep<=n_timesteps; timestep++) {
    
    std::cout << "starting timestep loop, timestep = " << timestep << std::endl;

    for (int rkstep=1; rkstep<=n_rksteps; rkstep++) {
      
      std::cout << "starting rk loop, rkstep = " << rkstep << std::endl;


      // Loop over domains   
      for (int domain=0; domain<domains.size(); domain++) {

	buffers recvbuffers;
	// Assemble recvbuffers
	// Loop over neighbors
	for (int neighbor=0; neighbor<domains[domain].numneighbors; neighbor++) {
	  int neighbor_here = domains[domain].neighbors[neighbor];
	  recvbuffers[neighbor_here]=domains[neighbor_here].sendbuffers[domain];
	}

	domains[domain].update(timestep,rkstep,recvbuffers);
	
      }// end loop over domains

      
      
    } // end rkstep loop
    
    //return 0;  // stop after one timestep for debugging


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


/*
no match for call to ‘
(std::map<int, std::vector<double> >::mapped_type {aka std::vector<double>}) 
(std::map<int, std::vector<double> >::mapped_type&)’
*/
