#ifndef DOMAIN_REFERENCE 
#define DOMAIN_REFERENCE

#include <hpx/config.hpp>
#include <hpx/util/itt_notify.hpp>
#include <hpx/include/performance_counters.hpp>

#include "../fname.h"
#include "../fortran_declarations.hpp"

std::vector<int> neighboringDomainIDs(void *size, void *dg, void *global)
{
    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    if (!size) {
        return std::vector<int>();
    }

    
    FNAME(get_neighbors_fort)(&size,
			      &dg,
			      &global,
			      neighbors_fort,
			      &numneighbors_fort);
    std::vector<int> ret(numneighbors_fort);
    std::copy(neighbors_fort, neighbors_fort + numneighbors_fort, ret.begin());

    return ret;
}

class DomainReference
{
public:
    const static int NANO_STEPS = 1;

    class API :
        public LibGeoDecomp::APITraits::HasUnstructuredTopology,
        public LibGeoDecomp::APITraits::HasNanoSteps<NANO_STEPS>,
	public LibGeoDecomp::APITraits::HasCustomMessageType<std::vector<double> >
    {};

    class FortranPointerWrapper
    {
    public:
        explicit
        FortranPointerWrapper(void *size, void *global, void *dg, void *nodalattr) :
            size(size),
	    global(global),
	    dg(dg),
	    nodalattr(nodalattr)
        {}
      
	~FortranPointerWrapper()
	{
	    if (size) {
	      FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    }
	}
      

	void *size;
	void *global;
	void *dg;
	void *nodalattr;

    };

    explicit
    DomainReference(int id = 0, 
		    void *size = 0, 
		    void *global = 0, 
		    void *dg = 0, 
		    void *nodalattr = 0, 
		    const std::vector<int>& neighbors_here = std::vector<int>()) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        id(id),
        timestep(1),
	rkstep(1),
	slopelimiter(false),
	substep(1),
        neighbors(neighbors_here),
	avg_slopelim_update(0.0),
	avg_rk_update(0.0)
    {
    }

#if defined(HPX_HAVE_ITTNOTIFY)
    char* get_function_name()
    {
	char* name = new char[120];
	std::strcpy(name, "update#");
	char buffer[33];
	//itoa(id, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", id);
	std::strcat(name, buffer);
	std::strcat(name, "#");
	//itoa(timestep, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", timestep);
	std::strcat(name, buffer);
	std::strcat(name, "#");
	//itoa(substep, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", substep);
	std::strcat(name, buffer);
	return name;
    }

    char* get_comm_name()
    {
	char* name = new char[120];
	std::strcpy(name, "comm#");
	char buffer[33];
	//itoa(id, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", id);
	std::strcat(name, buffer);
	std::strcat(name, "#");
	//itoa(timestep, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", timestep);
	std::strcat(name, buffer);
	std::strcat(name, "#");
	//itoa(substep, buffer, 10);
	snprintf(buffer, sizeof(buffer), "%d", substep);
	std::strcat(name, buffer);
	return name;
    }
#endif
    
  // //HPX Performance counter
  // boost::atomic<std::int64_t> counter_rk;
  // boost::atomic<std::int64_t> counter_slopelim;
  // boost::atomic<std::int64_t> rk_time_int;
  // boost::atomic<std::int64_t> slopelim_time_int;

  // //  std::int64_t rk_update_perf_data(double rk_time)
  // std::int64_t rk_update_perf_data(bool reset)
  // {
  //   std::int64_t counter = ++counter_rk;
  //   rk_time_int+=rk_time*1000000000;
  //   std::int64_t rk_time_int_now = rk_time_int;

  //   // return average time in nanoseconds
  //   return rk_time_int_now/counter;
  // }

  // //  std::int64_t slopelim_update_perf_data(double slopelim_time)
  // std::int64_t slopelim_update_perf_data(bool reset)
  // {
  //   std::int64_t counter = ++counter_rk;
  //   slopelim_time_int+=slopelim_time*1000000000;
  //   std::int64_t slopelim_time_int_now = slopelim_time_int;

  //   // return average time in nanoseconds
  //   return slopelim_time_int_now/counter;
  // }



  // void register_counter_type()
  // {
  //   hpx::performance_counters::install_counter_type(
  // 						    "/dgswem/rk_update_time",
  // 						    &DomainReference::rk_update_perf_data,
  // 						    "returns average time of rk update step"
  // 						    "ns"
  // 						    );
  //   hpx::performance_counters::install_counter_type(
  // 						    "/dgswem/slopelim_update_time",
  // 						    &DomainReference::slopelim_update_perf_data,
  // 						    "returns average time of rk update step"
  // 						    "ns"
  // 						    );
  // }
    
    template<typename HOOD, typename EVENT>
    void update(HOOD& hood, const EVENT& event)
    {

	
        int globalNanoStep = event.step() * NANO_STEPS + event.nanoStep();
      	
	if (timestep == 0) {
	  timestep++;
	  // hpx::cout << "skipping timestep 0" << std::endl;
	  return;
	} // skip timestep = 0

	switch(substep) {
	case 1: rkstep = 1;
	  slopelimiter = false;
	  break;
	case 2: rkstep = 1;
	  slopelimiter = true;
	  break;
	case 3: rkstep = 2;
	  slopelimiter = false;
	  break;
	case 4: rkstep = 2;
	  slopelimiter = true;
	  break;
	default:
	  hpx::cout << "ERROR!" << std::endl;
	}

		   


	/*if (id == 0) {
	  hpx::cout << " * globalNanoStep = " << globalNanoStep << " timestep = " << timestep << " substep = " << substep << 
	    " rkstep = " << rkstep << " slopelimiter = " << slopelimiter << std::endl;
	      
		if (timestep % 100 == 0)
		{
		}
	      
		}*/
	
	if (slopelimiter) {
	  
	  { // Open "update" code block (slopelimiter)
#if defined(HPX_HAVE_ITTNOTIFY)
	    std::unique_ptr<char> name(get_function_name());
	    static hpx::util::itt::domain d(hpx::get_thread_name().data());
	    //hpx::util::itt::task t(d, hpx::util::itt::string_handle(name.get()));
	    hpx::util::itt::task t(d, hpx::util::string_handle(name.get()));
#endif

	    hpx::util::high_resolution_timer slopelim_update_timer;
	    
	    //  if (id==0)
	  //  hpx::cout << "put element boundaries, timestep = " << timestep << std::endl;
	// Place incoming buffers into our cell 
	// Loop over neighbors
	  for (auto&& neighbor: neighbors) {
	    std::vector<double> buffer_vector = hood[neighbor];
	    int volume; // We don't really use this right now
	    double buffer[MAX_BUFFER_SIZE];
	    
	    // Make sure the buffer isn't empty, which it never should be
	    if (buffer_vector.size() == 0) {
	      throw std::logic_error("Empty incoming buffer!");
	    }	
	    // Make sure buffer isn't larger than MAX_BUFFER_SIZE
	    if (buffer_vector.size() > MAX_BUFFER_SIZE) {
	      throw std::logic_error("buffer_vector > MAX_BUFFER_SIZE!");
	    }
	    
	    // Load buffer into c-style array
	    for (int i=0; i<buffer_vector.size(); i++) {
	      buffer[i] = buffer_vector[i];
	    }
	    

	    int rkindex = 0;
	    
	    FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
				      &neighbor,
				      &volume,
				      buffer,
				      &rkindex);
	  
	  } // End loop over neighbors  
	  
	  //if (id==0)
	    //  hpx::cout << "call slopelimiter_parta, timestep = " << timestep << std::endl;
	  FNAME(slopelimiter_parta_fort)(&domainWrapper->size,
					 &domainWrapper->dg,
					 &domainWrapper->global);
	  
	  double elapsed = slopelim_update_timer.elapsed();
	  double total = avg_slopelim_update*total_slopelim_updates;
	  total_slopelim_updates++;
	  avg_slopelim_update = (elapsed + total)/total_slopelim_updates;

	}// end "update" code block (slopelimiter)


	{// begin "comm" code block (slopelimiter)

#if defined(HPX_HAVE_ITTNOTIFY)
	  std::unique_ptr<char> name(get_comm_name());
	  static hpx::util::itt::domain d(hpx::get_thread_name().data());
	  //hpx::util::itt::task t(d, hpx::util::itt::string_handle(name.get()));
	  hpx::util::itt::task t(d, hpx::util::string_handle(name.get()));
#endif
	  //if (id==0)
	    // hpx::cout << "get node boundaries, timestep = " << timestep << std::endl;
	  for (auto&& neighbor: neighbors) {
	    std::vector<double> send_buffer;
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];

	    
	    FNAME(hpx_get_nodes_fort)(&domainWrapper->dg, //pointer to current domain
				      &neighbor, // pointer to neighbor to send to
				      &volume,
				      buffer);

	    for (int i=0; i<volume; i++) {
	      send_buffer.push_back(buffer[i]);
	    }
	    
	    hood.send(neighbor, send_buffer);
	  } // End get node boundaries
	}// end "comm" code block (slopelimiter)
	  
	} else { // RK hydro step


	  {//begin "update" code block (rk step)
	    
#if defined(HPX_HAVE_ITTNOTIFY)
	    std::unique_ptr<char> name(get_function_name());
	    static hpx::util::itt::domain d(hpx::get_thread_name().data());
	    //hpx::util::itt::task t(d, hpx::util::itt::string_handle(name.get()));
	    hpx::util::itt::task t(d, hpx::util::string_handle(name.get()));
#endif

	    
	    hpx::util::high_resolution_timer rk_update_timer;
	 
	  int timestep_here;
	  if (rkstep == 1) {
	    timestep_here=timestep-1;
	  } else {
	    timestep_here=timestep;
	  }


	  if ( (timestep > 1) || (rkstep == 2) ){	  

	    // Put node boundaries
	    //if (id==0)
	      //hpx::cout << "put node boundaries, timestep = " << timestep_here << std::endl;	      	      
	    
	    for (auto&& neighbor: neighbors) {
	      std::vector<double> buffer_vector = hood[neighbor];
	      int volume; // We don't really use this right now
	      double buffer[MAX_BUFFER_SIZE];
	      
	      // Make sure the buffer isn't empty, which it never should be
	      if (buffer_vector.size() == 0) {
		throw std::logic_error("Empty incoming buffer!");
	      }	
	      // Make sure buffer isn't larger than MAX_BUFFER_SIZE
	      if (buffer_vector.size() > MAX_BUFFER_SIZE) {
		throw std::logic_error("buffer_vector > MAX_BUFFER_SIZE!");
	      }
	      
	      // Load buffer into c-style array
	      for (int i=0; i<buffer_vector.size(); i++) {
		buffer[i] = buffer_vector[i];
	      }
	      
	      FNAME(hpx_put_nodes_fort)(&domainWrapper->dg,
					&neighbor,
					&volume,
					buffer);
	      
	    } // End put node boundaries
	    
	    
	    // Call slopelimiter_part2
	    //if (id==0)
	      //hpx::cout << "Call slopelimiter_partb (with w/d), timestep = " << timestep_here << std::endl;	      	      
	    FNAME(slopelimiter_partb_fort)(&domainWrapper->size,
					   &domainWrapper->dg,
					   &domainWrapper->global);

	    	    	    
	    // if rkstep == 1, call dg_timestep_advance
	    if (rkstep == 1) {
	      //if (id==0)
		//hpx::cout << "dg_timestep_advance, timestep = " << timestep_here << std::endl;
	      
	      
		FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
						&domainWrapper->dg,
						&domainWrapper->global,
						&domainWrapper->nodalattr,
						&timestep_here
						);
	      
	    }
	    
	  } // end if timestep > 1 || rkstep == 2 
	  
	  // if (id==0)
	  //  hpx::cout << "dg_hydro_timestep, rkstep = " << rkstep << " timestep = " << timestep << std::endl;
	  // Do hydro
	  
	  FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
					&domainWrapper->dg,
					&domainWrapper->global,
					&domainWrapper->nodalattr,
					&timestep,
					&rkstep
					);
	  
	  double elapsed = rk_update_timer.elapsed();
	  double total = avg_rk_update*total_rk_updates;
	  total_rk_updates++;
	  avg_rk_update = (elapsed + total)/total_rk_updates;
	  
	}//end "update" code block (rk step)

	  {// begin "comm" code block (rk step)	  
	  
#if defined(HPX_HAVE_ITTNOTIFY)
	    std::unique_ptr<char> name(get_function_name());
	    static hpx::util::itt::domain d(hpx::get_thread_name().data());
	    //	    hpx::util::itt::task t(d, hpx::util::itt::string_handle(name.get()));
	    hpx::util::itt::task t(d, hpx::util::string_handle(name.get()));
#endif


	  //	  if (id==0)
	  //  hpx::cout << "Get element boundaries, timestep = " << timestep << std::endl;
	  
	  // Get element boundaries
	  for (auto&& neighbor: neighbors) {
	    std::vector<double> send_buffer;
	    int volume;
	    double buffer[MAX_BUFFER_SIZE];
	    int rkindex;
	    //if (rkstep == 1) rkindex = 2;
	    //if (rkstep == 2) rkindex = 3;
	    rkindex = 0;

	    
	    FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
				      &neighbor, // pointer to neighbor to send to
				      &volume,
				      buffer,
				      &rkindex);
	    
	    for (int i=0; i<volume; i++) {
	      send_buffer.push_back(buffer[i]);
	    }
	    
	    hood.send(neighbor, send_buffer);
	  } // End get element boundaries

	} // end slopelimiter == false	(RK step)
	
	  // Control flow
	if (substep == 4) {
	  substep = 1;
	  timestep++;
	} else {
	  substep++;
	}
	
	}// end "comm" code block (rk step)
	
    } // End of update() ------------------------------------------------------------------
  
  
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, unsigned)
  {
    ar & id & timestep & rkstep; // DO THIS
  }
  
private:
  boost::shared_ptr<FortranPointerWrapper> domainWrapper;
  int id;
  int timestep;
  int rkstep;
  bool slopelimiter;
  int substep;
  std::vector<int> neighbors;


  int total_rk_updates;
  double avg_rk_update;
  int total_slopelim_updates;
  double avg_slopelim_update;
};

LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(DomainReference)
#endif
