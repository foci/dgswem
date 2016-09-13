#define MAX_DOMAIN_NEIGHBORS 100
#define MAX_BUFFER_SIZE 10000
#define MAX_CELL_SIZE 20

extern"C" {
    void FNAME(dgswem_init_fort)(void** size,
				 void** dg,
				 void** global,
				 void** nodalattr,
				 int* id);
    void FNAME(dg_hydro_timestep_fort)(void** size,
				       void** dg,
				       void** global,
				       void** nodalattr,
				       int* timestep,
				       int* rkstep);	
   void FNAME(slopelimiter_partA_fort)(void** size,
				       void** dg,
				       void** global);
   void FNAME(slopelimiter_partB_fort)(void** size,
				       void** dg,
				       void** global);
    void FNAME(dg_timestep_advance_fort)(void** size,
					 void** dg,
					 void** global,
					 void** nodalattr,
					 int* timestep
					 );	
    void FNAME(get_neighbors_fort)(void** size,
				   void** dg,
				   void** global,
				   int neighbors[MAX_DOMAIN_NEIGHBORS],
				   int* num_neighbors);	
    void FNAME(hpx_get_elems_fort)(void** dg,
				   int* neighbor,
				   int* volume,
				   double* sendbuf,
				   int* rkindex);
    void FNAME(hpx_put_elems_fort)(void** dg,
				   int* neighbor,
				   int* volume,
				   double* recvbuf,
				   int* rkindex);
    void FNAME(hpx_get_nodes_fort)(void** dg,
				   int* neighbor,
				   int* volume,
				   double* sendbuf); 
    void FNAME(hpx_put_nodes_fort)(void** dg,
				   int* neighbor,
				   int* volume,
				   double* recvbuf);    
    void FNAME(hpx_read_n_domains)(int* n_domains);
    void FNAME(hpx_swap_elems_fort)(void** dg_domain,
				    void** dg_neighbor);
    void FNAME(lgd_yield_subdomain_coord)(void** global,
					  double *coord);
    void FNAME(term_fort)(void** size,
			  void** dg,
			  void** global,
			  void** nodalattr);			  
    void FNAME(cpp_vars_from_fort)(void** size,				   
				   int* rkstep,
				   int* timestep);
    void FNAME(cpp_vars_to_fort)(void** size,
				   int* rkstep,
				   int* timestep);
    void FNAME(check_c_ptr)(void** cptr,
			    int* retval);
}
