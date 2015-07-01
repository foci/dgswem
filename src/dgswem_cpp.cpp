#include <iostream>
#include "fname.h"

extern"C" {
  void FNAME(dgswem_init_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* n_timesteps);
  void FNAME(dg_timestep_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* timestep);	
}		       

int main(
	 int argc
	 , char* argv[]
	 )
{

  void *sizes = NULL;
  void *dg = NULL;
  void *global = NULL;
  void *nodalattr = NULL;
  int n_timesteps;

  FNAME(dgswem_init_fort)(&sizes,&dg,&global,&nodalattr,&n_timesteps);

  std::cout << "n_timesteps = " << n_timesteps << std::endl;

  for (int timestep=1; timestep<=n_timesteps; timestep++) {
    FNAME(dg_timestep_fort)(&sizes,&dg,&global,&nodalattr,&timestep);
  }
}
