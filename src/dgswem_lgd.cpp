#include <iostream>
#include <limits>
#include <vector>
#include <libgeodecomp.h>
#include <libgeodecomp/geometry/unstructuredgridmesher.h>
#include <libgeodecomp/io/logger.h>

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
	    FNAME(term_fort)(&size,&global,&dg,&nodalattr);
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
        timestep(0)
    {}

    template<typename HOOD>
    void update(const HOOD& hood, int nanoStep)
    {
        std::cout << "updating...\n";
	FNAME(dg_hydro_timestep_fort)(&domainWrapper->size,
				&domainWrapper->dg,
				&domainWrapper->global,
				&domainWrapper->nodalattr,
				&timestep,
				&rkstep
				);

        std::cout << "  done\n";

	// TODO: add boundary exchange
        // Boundary exchange
	/*
        //Loop over neighbors
        for (int neighbor=0; neighbor<neighbors_here.size(); neighbor++) {
            std::cout << "  comm " << id << "<->" << neighbors_here[neighbor] << "\n";
            int neighbor_here = neighbors_here[neighbor];
            int num_nodes;
            int alive[MAX_BOUNDARY_SIZE];

            std::cout << "  pull\n";
            // Get outgoing boundarys from the neighbors
            FNAME(get_outgoing_nodes_fort)(&hood[neighbor_here].domainWrapper->domain,&id,&num_nodes,alive);

            std::cout << "  push\n";
            // Put those arrays inside current domain
            FNAME(put_incoming_nodes_fort)(&domainWrapper->domain,&neighbor_here,&num_nodes,alive);
            std::cout << "  done\n";
        }// end loop over neighbors
	*/
	

	// TODO add some logic for the RK_step and calling dg_timestep_advance_fort
        ++timestep;
    }

private:
    boost::shared_ptr<FortranPointerWrapper> domainWrapper;
    std::vector<int> neighbors_here;
    int id;
    int timestep;
    int rkstep;
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
        //std::fill(&base_path_buf[0], &base_path_buf[1024], 0);
        //std::copy(base_path.begin(), base_path.end(), base_path_buf);

        // Initialize domain decomposition
        std::vector<std::vector<int> > neighbors;

        for(int i = 0; i < numDomains ; i++) {

	    // Initialize the domains temporarily to get grid information

	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    int n_domains;
	    int n_rksteps;
	    
	    LibGeoDecomp::FloatCoord<2> coord;

	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &n_domains,
				    &i,
				    &n_rksteps);

            FNAME(lgd_yield_subdomain_coord)(&global, &coord[0]);
            domainCoords << coord;
	    
	    neighbors.push_back(neighboringDomainIDs(&size, &dg, &global));

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

	    int n_domains;
	    int n_rksteps;

	    FNAME(dgswem_init_fort)(&size,
				    &dg,
				    &global,
				    &nodalattr,
				    &n_domains,
				    &id,
				    &n_rksteps);
	    
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

private:
    char base_path_buf[1024];
    LibGeoDecomp::UnstructuredGridMesher<2> mesher;
    std::size_t numDomains;
    std::vector<LibGeoDecomp::FloatCoord<2> > domainCoords;
};

int main(int argc, char* argv[])
{
    // fixme: these should be read in via config files
    int n_domains = 32;
    int n_timesteps = 100;

    LibGeoDecomp::SerialSimulator<FortranCell> sim(
        new FortranInitializer(n_domains, n_timesteps));
    sim.run();

    return 0;
}
