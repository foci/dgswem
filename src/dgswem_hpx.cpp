#include <iostream>
#include <limits>
#include <vector>
#include <libgeodecomp.h>
#include <libgeodecomp/geometry/unstructuredgridmesher.h>
#include <libgeodecomp/io/logger.h>

// Needed by HPX
#include <libgeodecomp/parallelization/hpxsimulator.h>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <libgeodecomp/loadbalancer/tracingbalancer.h>
#include <libgeodecomp/geometry/partitions/recursivebisectionpartition.h>


#include "fname.h"
#include "fortran_declarations.hpp"

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
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighboringDomainIDs(size,dg,global)),
        id(id),
        timestep(0),
	rkstep(1),
	update_step(true),
	exchange_step(false),
	advance_step(false)
    {}

    template<typename HOOD>
    void update(const HOOD& hood, int nanoStep)
    {
	if (timestep != 0) {

	    if (update_step) {
		/*
		std::cout << "updating (domain_id = " << id
			  << ", timestep = " << timestep
			  << ", rkstep = " << rkstep
			  << ")...\n";
		*/
		FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
					      &domainWrapper->dg,
					      &domainWrapper->global,
					      &domainWrapper->nodalattr,
					      &timestep,
					  &rkstep
					      );
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
		    
		    
		    /*
		      std::cout << "domain " << id << " is exchanging with " << neighbor_here
			      << std::endl;
		    */
		    // Get outgoing boundarys from the neighbors
		    FNAME(hpx_get_elems_fort)(&hood[neighbor_here].domainWrapper->dg,
					      &id,
					      &volume,
					      buffer);
		    
		    
		    // Put those arrays inside current domain
		    FNAME(hpx_put_elems_fort)(&domainWrapper->dg,
					      &neighbor_here,
					      &volume,
					      buffer);
		    
		    
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
		//std::cout << "advancing domain " << id << std::endl;
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
    }


    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	ar & domainWrapper & neighbors_here & id & timestep & rkstep & update_step & exchange_step & advance_step;
    }

private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    std::vector<int> neighbors_here;
    int id;
    int timestep;
    int rkstep;
    bool update_step;
    bool exchange_step;
    bool advance_step;

};

typedef LibGeoDecomp::ContainerCell<DomainReference, MAX_CELL_SIZE, int> FortranCell;

class FortranInitializer : public LibGeoDecomp::SimpleInitializer<FortranCell>
{
public:
    using SimpleInitializer<FortranCell>::dimensions;

    FortranInitializer(std::size_t numDomains, std::size_t numSteps) :
        SimpleInitializer<FortranCell>(LibGeoDecomp::Coord<2>(), numSteps),
        numDomains(numDomains)
    {
        // Initialize domain decomposition
        std::vector<std::vector<int> > neighbors;

        for(int i = 0; i < numDomains ; i++) {
	    // Initialize the domains temporarily to get grid information

	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    int n_domains_fort;
	    int n_rksteps_fort;
	    
	    LibGeoDecomp::FloatCoord<2> coord;

	    int domain_number = i;
	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &domain_number);

            FNAME(lgd_yield_subdomain_coord)(&global, &coord[0]);
            domainCoords << coord;
	    neighbors.push_back(neighboringDomainIDs(size, dg, global));

	    //destroy these domains
	    FNAME(term_fort)(&size,&global,&dg,&nodalattr);
        }

        mesher = LibGeoDecomp::UnstructuredGridMesher<2>(domainCoords, neighbors);
        dimensions = mesher.logicalGridDimension();

        LOG(INFO, "FortranInitializer(\n"
            << "  dimensions: " << dimensions << "\n"
            << "  cellDimensions: " << mesher.cellDimension() << ")\n");
    }

    void grid(LibGeoDecomp::GridBase<FortranCell, 2> *grid)
    {
        LibGeoDecomp::CoordBox<2> box = grid->boundingBox();

        for(int id = 0; id < numDomains ; id++) {
            // Create vectors of domain pointers and ids
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
	    
            LibGeoDecomp::FloatCoord<2> coord;
            FNAME(lgd_yield_subdomain_coord)(&global, &coord[0]);
            LibGeoDecomp::Coord<2> logicalCoord = mesher.positionToLogicalCoord(coord);

            if (box.inBounds(logicalCoord)) {
                LOG(INFO, "adding domain " << id
                    <<  " to grid cell " << logicalCoord
                    << " (is at position " << coord << ")\n");
                FortranCell cell = grid->get(logicalCoord);
                cell.insert(id, DomainReference(id, size, global, dg, nodalattr));
                grid->set(logicalCoord, cell);
            } else {
		FNAME(term_fort)(&size,&global,&dg,&nodalattr);
            }
        }
    }

    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	ar & boost::serialization::base_object<SimpleInitializer<FortranCell> >(*this);
    }

private:
    char base_path_buf[1024];
    LibGeoDecomp::UnstructuredGridMesher<2> mesher;
    std::size_t numDomains;
    std::vector<LibGeoDecomp::FloatCoord<2> > domainCoords;

};

typedef LibGeoDecomp::HpxSimulator::HpxSimulator<FortranCell, RecursiveBisectionPartition<s> > SimulatorType;

int main(int argc, char* argv[])
{
    // todo: these should be read in via config files
    int n_domains;
    FNAME(hpx_read_n_domains)(&n_domains);

    int n_timesteps = 86401;
    //int n_timesteps = 4001;
    //int n_timesteps = 2;

    int n_rksteps = 2;


    // Needed by HPX simulator
    std::vector<double> updateGroupSpeeds(1, 1.0);
    int ghostZoneWidth = 1;
    FortranInitializer *init = new FortranInitializer(n_domains, total_rksteps);

    int total_rksteps = n_timesteps*(n_rksteps*2+1)+1;
    SimulatorType sim(init,
		      updateGroupSpeeds, 
		      new TracingBalancer(new OozeBalancer()),
		      int loadBalancingPeriod = 10,
		      ghostZoneWidth,
		      "dgswem-hpx");
    sim.run();

    return 0;
}
