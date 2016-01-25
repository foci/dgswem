#include <iostream>
#include <limits>
#include <vector>
#include <libgeodecomp.h>
#include <libgeodecomp/io/logger.h>
#include <libgeodecomp/storage/unstructuredgrid.h>
#include <libgeodecomp/misc/apitraits.h>

// Needed by HPX
#include <libgeodecomp/parallelization/hpxsimulator.h>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/iostreams.hpp>
#include <libgeodecomp/loadbalancer/tracingbalancer.h>
#include <libgeodecomp/geometry/partitions/recursivebisectionpartition.h>

using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

#include "fname.h"
#include "fortran_declarations.hpp"

//Defining settings for SELL-C-q
static const std::size_t MATRICES = 1;
static const int C = 4;
static const int SIGMA = 1;

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
    class API :
        public LibGeoDecomp::APITraits::HasUnstructuredTopology,
	public LibGeoDecomp::APITraits::HasSellType<double>,
        public LibGeoDecomp::APITraits::HasSellMatrices<MATRICES>,
        public LibGeoDecomp::APITraits::HasSellC<C>,
        public LibGeoDecomp::APITraits::HasSellSigma<SIGMA>
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
	      //	      std::cout << "CPP: about to call term_fort" << std::endl;
		FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    }
	}
      
      
	void *size;
	void *global;
	void *dg;
	void *nodalattr;

    };

    explicit
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0, bool busywork = false) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighboringDomainIDs(size,dg,global)),
        id(id),
        timestep(0),
	rkstep(1),
	update_step(true),
	exchange_step(false),
	advance_step(false),
	busywork(busywork)
  {}

    template<typename HOOD>
    void update(const HOOD& hood, int nanoStep)
    {
      
      *this = hood[hood.index()];

      bool only_busywork = this->busywork;
      
      if (only_busywork) {
	
	// Adding some busy work FIXME
	volatile double a;
	volatile double b;
	volatile double c;
	a = 23.2394908;
	std::cout << "timestep = " << timestep << " domainID = " << id << std::endl;
	for (int i=0; i<1000000; i++) {
	  //		for (int j=0; j<10000000; j++) {
	  b = 2.23423;
	  c = a/b;
	  c = a*(b+0.000001);
	  a = c;		  
	  //		}

	}
	timestep++;
      
      } else {
	
	//      std::cout << "CPP: LGD update" << std::endl;
	if (timestep != 0) {
	  
	  if (update_step) {

	      // Clear output buffer map
	      output_buffer.clear();
	    
	      //std::cout << "********* update step ************" << std::endl;

	      /*
	      std::cout << "updating (domain_id = " << id
			  << ", timestep = " << timestep
			  << ", rkstep = " << rkstep
			  << ")...\n";
	      */

		//                 std::cout << "CPP: about to call dg_hydro_timestep_fort" << std::endl;
		FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
					      &domainWrapper->dg,
					      &domainWrapper->global,
					      &domainWrapper->nodalattr,
					      &timestep,
					      &rkstep
					      );
	      
		//std::cout << "a=" << a << std::endl;
		
		
		// Fill output buffer
		//Loop over neighbor domains
		for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
		    int neighbor_here = neighbors_here[neighbor];
		    int volume;
		    double buffer[MAX_BUFFER_SIZE];
		    FNAME(hpx_get_elems_fort)(&domainWrapper->dg, //pointer to current domain
					      &neighbor_here, // pointer to neighbor to send to
					      &volume,
					      buffer);
		    
		    //std::cout << "volume = " << volume << std::endl;

		    //std::cout << "******************************************************" << std::endl;

		    //std::map<int,std::vector<double> > output_buffer;
		    //std::map<int,std::vector<double> > test_map;

		    
		    //std::cout << "Buffer from " << id << ", for " << neighbor_here << std::endl;
		    //std::cout << "volume = " << volume << std::endl;
		    // Pack buffer into std vector
		    std::vector<double> buffer_vector;
		    for (int i=0; i<volume; i++) {
			//std::cout << "buffer[" << i << "] = " << buffer[i] << " ";
			buffer_vector.push_back(buffer[i]);
			
		    }
		    //std::cout << std::endl;

		    /*
		    std::cout << "Buffer Vector:" << std::endl;
		    std::cout << "buffer_vector.size() = " << buffer_vector.size() << std::endl;
		    for (int i=0; i<buffer_vector.size(); i++) {
			std::cout << "buffer_vector[" << i << "] = " << buffer_vector[i] << " ";
		    }
		    std::cout << std::endl;
		    */
		    		   
		    std::pair<std::_Rb_tree_iterator<std::pair<const int, std::vector<double> > >, bool> outval = output_buffer.insert(std::map<int, std::vector<double> >::value_type(neighbor_here, buffer_vector));
		    if (std::get<1>(outval)) {
			//std::cout << "Insert successful!" << std::endl;
		    } else {
			std::cout << "Insert not successful!" << std::endl;
		    }

		    //std::cout << "Vector retrieved from map:" << std::endl;
		    //std::vector<double> temp = output_buffer.at(neighbor_here);
		    /*
		    std::cout << "temp.size() = " << temp.size() << std::endl;
		    for (int i=0; i<temp.size(); i++) {	
			std::cout << "temp[" << i << "] = " << temp[i] << " "; 
		    }
		    std::cout << std::endl;
		    std::cout << "******************************************************" << std::endl;
		    */

		    //std::cout << "output_buffer.at(neighbor_here).size() = " << output_buffer.at(neighbor_here).size() << std::endl;
		    //std::cout << "temp.size() = " << temp.size() << std::endl;
		}

		

		//std::cout << "output_buffer.size() = " << output_buffer.size() << std::endl;

		update_step = false;
		exchange_step = true;
		advance_step = false;
	    } else if (exchange_step) {

		// Boundary exchange
		//Loop over neighbors
		for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
		    int neighbor_here = neighbors_here[neighbor];
		    int volume;
		    double buffer[MAX_BUFFER_SIZE];

		    // Unpack buffer from neighbor
		    //std::cout << "id here:" << id << ", unpacking buffer from " << neighbor_here << std::endl;
		    //const std::vector<double>*  buffer_vector = &hood[neighbor_here].output_buffer.at(id);
		    //std::vector<double> buffer_vector = &hood[neighbor_here].output_buffer.at(id);
		    std::vector<double> buffer_vector = hood[neighbor_here].output_buffer.at(id);
		    //std::vector<double> buffer_vector = output_buffer.at(neighbor_here);
		    //for (int i=0; i<buffer_vector->size(); i++) {
		    for (int i=0; i<buffer_vector.size(); i++) {
			buffer[i] = buffer_vector[i];
			//std::cout << "buffer[" << i << "] = " << buffer[i] << " "; 
		    }
		    //std::cout << std::endl;

		    // Put elements into our own subdomain
		    FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
					      &neighbor_here,
					      &volume,
					      buffer);
		    
		    //std::cout << "domain " << id << " is exchanging with " << neighbor_here << " at timestep " << timestep << std::endl;

		    // Get outgoing boundarys from the neighbors	
		    //		    std::cout << "CPP: about to call hpx_get_elems_fort" << std::endl;		   
		    /*
		    FNAME(hpx_get_elems_fort)(&hood[neighbor_here].domainWrapper->dg,
					      &id,
					      &volume,
					      buffer);
		    */
		    
		    // Put those arrays inside current domain
		    //		    std::cout << "CPP: about to call hpx_put_elems_fort" << std::endl;
		}// end loop over neighbors

		exchange_step = false;

		if (rkstep == 2) {
		    advance_step = true;
		    update_step = false;
		    rkstep = 1;
		} else {
		    ++rkstep;
		    update_step = true;
		    advance_step = false;
		}
		
		    
	    } else if (advance_step) {	
	      //std::cout << "advancing domain " << id << " at timestep " << timestep <<std::endl;
	      //	      std::cout << "CPP: about to call dg_timestep_advance_fort" << std::endl;
		FNAME(dg_timestep_advance_fort)(&domainWrapper->size,
						&domainWrapper->dg,
						&domainWrapper->global,
						&domainWrapper->nodalattr,
						&timestep
						);
		++timestep;
		
		update_step = true;
		exchange_step = false;
		advance_step = false;
	    } else {
		std::cout << "Something went wrong!" << std::endl;
	    }
	    	    
	} else {
	    ++timestep;
	}

      } // else busywork

    }

    // Serialization function
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	//throw std::runtime_error("no serialization yet!");
	ar & neighbors_here & id & timestep & rkstep & update_step & exchange_step & advance_step & output_buffer & busywork;
    }

    /*
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	ar & domainWrapper & neighbors_here & id & timestep & rkstep & update_step & exchange_step & advance_step;
    }
    */


private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    std::vector<int> neighbors_here;
    std::map<int,std::vector<double> > output_buffer;
    int id;
    int timestep;
    int rkstep;
    bool update_step;
    bool exchange_step;
    bool advance_step;
    bool busywork;
};

//typedef LibGeoDecomp::ContainerCell<DomainReference, MAX_CELL_SIZE, int> FortranCell;

class FortranInitializer : public LibGeoDecomp::SimpleInitializer<DomainReference>
{
public:
    //using SimpleInitializer<DomainReference>::dimensions;

    typedef LibGeoDecomp::UnstructuredGrid<DomainReference, MATRICES, double, C, SIGMA> Grid;

  FortranInitializer(std::size_t numDomains, std::size_t numSteps, bool busywork) :
        SimpleInitializer<DomainReference>(LibGeoDecomp::Coord<1>(numDomains), numSteps),
        numDomains(numDomains),
	busywork(busywork)
    {
	// Empty
    }

    LibGeoDecomp::Adjacency getAdjacency() const
    {
	LibGeoDecomp::Adjacency adjacency;
        for(int id = 0; id < numDomains ; id++) {
	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    int domain_number = id;
	    
	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &domain_number);

	    // Get list of neighbors for this domain
	    std::vector<int> neighbors_here = neighboringDomainIDs(size, dg, global);
	    // Destroy the domain
	    FNAME(term_fort)(&size,&global,&dg,&nodalattr);
	    adjacency[id] = neighbors_here;
	    }
	return adjacency;
    }

    void grid(LibGeoDecomp::GridBase<DomainReference, 1> *grid)
    {
	

	// The "double" below is just the weighting, 
	// we can ignore this for DGSWEM
        std::map<LibGeoDecomp::Coord<2>, double> adjacency;	
	
	// Loop only over the domains we need to
	int origin = grid->boundingBox().origin[0];
	int end = origin + grid->boundingBox().dimensions[0];
	for (int id = origin; id < end; ++id) {
	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    int domain_number = id;
	    
	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &domain_number);

	    // Get list of neighbors for this domain
	    std::vector<int> neighbors_here = neighboringDomainIDs(size, dg, global);
	    
	    // Loop over all the neighbors
	    for (int i=0; i<neighbors_here.size(); i++) {
		// Adding adjacency information - pairs of neighbors
		// 1.0 is the weighting, which is unused for DGSWEM
		adjacency[LibGeoDecomp::Coord<2>(id, neighbors_here[i])] = 1.0;
	    }
	    
	    DomainReference cell(id, size, global, dg, nodalattr, busywork);
	    grid->set(LibGeoDecomp::Coord<1>(id), cell);
        }
	Grid *unstructuredgrid = dynamic_cast<Grid *>(grid);
	unstructuredgrid->setAdjacency(0, adjacency);
    }

    // Serialization function
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	throw std::runtime_error("no serialization yet!");
    }

private:
    std::size_t numDomains;
  bool busywork;
};

//#define LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(CARGO)			\
//    typedef LibGeoDecomp::HPXReceiver<CARGO>::receiveAction DummyReceiver_ ## CARGO ## _ReceiveAction; \
//    HPX_REGISTER_ACTION_DECLARATION(DummyReceiver_ ## CARGO ## _ReceiveAction); \
//    HPX_REGISTER_ACTION(DummyReceiver_ ## CARGO ## _ReceiveAction);	\

//    typedef hpx::components::simple_component<LibGeoDecomp::HPXReceiver<CARGO> > receiver_type_ ## CARGO; \
//    HPX_REGISTER_COMPONENT(receiver_type_ ## CARGO , DummyReceiver_ ## CARGO);

LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(DomainReference)

//typedef LibGeoDecomp::CoordBox<2> CoordBoxType;
//LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(CoordBoxType)

//typedef std::vector<LibGeoDecomp::ContainerCell<DomainReference, 20ul, int> > ContainerCellType;
//LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(std::vector<DomainReference>)

//This is a hack
typedef LibGeoDecomp::CoordBox<1> CoordBox1;
LIBGEODECOMP_REGISTER_HPX_COMM_TYPE(CoordBox1)

//Change me:
typedef LibGeoDecomp::HpxSimulator<DomainReference, LibGeoDecomp::UnstructuredStripingPartition> SimulatorType;
//typedef LibGeoDecomp::HpxSimulator<DomainReference, LibGeoDecomp::RecursiveBisectionPartition<2> > SimulatorType;

int hpx_main(variables_map & vm)
{

  
  
  // number of domains read from config file
    int n_domains;
    //    std::cout << "CPP: about to call hpx_read_n_domains" << std::endl;
    FNAME(hpx_read_n_domains)(&n_domains);

    //    int n_timesteps = 86401;
    //int n_timesteps = 2000;
    //int n_timesteps = 8000;
    //int n_timesteps = 100;
    //int n_timesteps = 2;

    int n_timesteps = vm["n_timesteps"].as<std::size_t>();
    bool busywork = vm["busywork"].as<bool>();

    std::cout << "n_timesteps from cmd line = " << n_timesteps << std::endl;
    std::cout << "busywork from cmd line = " << busywork << std::endl;

    int n_rksteps = 2;


    // Needed by HPX simulator
    std::vector<double> updateGroupSpeeds(1, 1.0);
    int ghostZoneWidth = 1;
    int total_rksteps = n_timesteps*(n_rksteps*2+1)+1;
    std::cout << "total_rksteps = " << total_rksteps << std::endl;
    //std::exit(1);
    //total_rksteps = 12;

    FortranInitializer *init = new FortranInitializer(n_domains, total_rksteps, busywork);

    SimulatorType sim(
		      init,
		      updateGroupSpeeds,
		      new LibGeoDecomp::TracingBalancer(new LibGeoDecomp::OozeBalancer()),
		      1, // Just changed this from 10 to 1. 
		      ghostZoneWidth,
		      "dgswem-hpx");
    sim.run();

    return hpx::finalize();
}

int main(int argc, char **argv)
{

  // Parse command line options
  options_description
    desc_commandline("usage: " HPX_APPLICATION_STRING " [options]");

  desc_commandline.add_options()
    (
            "n_timesteps"
	    , value<std::size_t>()->default_value(100)
	    , "Number of timesteps"
     )
    (
            "busywork"
	    , value<bool>()->default_value(false)
	    , "Option to use busywork instead of calling fortran subroutines"
     )
    ;

    // We want HPX to run hpx_main() on all localities to avoid the
    // initial overhead caused by broadcasting the work from one to
    // all other localities:
    std::vector<std::string> config(1, "hpx.run_hpx_main!=1");


    // fixme: I'm not sure that config is being used by hpx_main
    return hpx::init(desc_commandline,argc,argv,config);

}
