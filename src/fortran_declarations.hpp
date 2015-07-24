#define MAX_DOMAIN_NEIGHBORS 10
#define MAX_BUFFER_SIZE 1000

extern"C" {
  void FNAME(dgswem_init_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* n_domains,
			       int* id,
			       int* n_rksteps);
  void FNAME(dg_hydro_timestep_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* timestep,
			       int* rkstep);	
  void FNAME(dg_timestep_advance_fort)(void** sizes,
			       void** dg,
			       void** global,
			       void** nodalattr,
			       int* timestep,
			       int* rkstep);	
  void FNAME(get_neighbors_fort)(void** sizes,
				 void** dg,
				 void** global,
				 int neighbors[MAX_DOMAIN_NEIGHBORS],
				 int* num_neighbors);	
  void FNAME(hpx_get_elems_fort)(void** dg,
				 int* neighbor,
				 int* volume,
				 double* sendbuf);
  void FNAME(hpx_put_elems_fort)(void** dg,
				 int* neighbor,
				 int* volume,
				 double* recvbuf);
  void FNAME(hpx_read_n_domains)(int* n_domains);
  void FNAME(hpx_swap_elems_fort)(void** dg_domain,
                                  void** dg_neighbor);  
}
