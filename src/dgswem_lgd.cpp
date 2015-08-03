#include <iostream>
#include <limits>
#include <vector>
#include <libgeodecomp.h>
#include <libgeodecomp/geometry/unstructuredgridmesher.h>
#include <libgeodecomp/io/logger.h>
#include "fname.h"

#include "fortran_includes.hpp"


std::vector<int> neighboringDomainIDs(void *size, void *dg, void *global)
{

    //std::cout << "twerk1\n";
    int numneighbors_fort;
    int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    //std::cout << "twerk2 " << domain << "\n";
    if (!domain) {
        return std::vector<int>();
    }

    FNAME(get_neighbors_fort)(&size,
			      &dg,
			      &global,
			      neighbors_fort,
			      &numneighbors_fort);
    //std::cout << "twerk3\n";
    std::vector<int> ret(numneighbors_fort);
    //std::cout << "twerk4\n";
    std::copy(neighbors_fort, neighbors_fort + numneighbors_fort, ret.begin());
    //std::cout << "twerk5\n";

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
      
      /* TODO: write a destructor
	 ~FortranPointerWrapper()
	 {
	 if (domain) {
	 FNAME(term_fort)(&domain);
	 }
	 }
      */
      
      
	void *size;
	void *global;
	void *dg;
	void *nodalattr;

    };

    explicit
    DomainReference(int id = 0, void *size = 0, void *global = 0, void *dg = 0, void *nodalattr = 0) :
        domainWrapper(new FortranPointerWrapper(size, global, dg, nodalattr)),
        neighbors_here(neighboringDomainIDs(domain)),
        id(id),
        timestep(0)
    {}

    template<typename HOOD>
    void update(const HOOD& hood, int nanoStep)
    {
        std::cout << "updating...\n";
        // Update values using Game of Life rules
        //FNAME(update_fort)(&timestep, &domainWrapper->domain);
	FNAME(dg_timestep_fort)(&domainWrapper->size,
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

    FortranInitializer(const std::string& base_path, std::size_t numDomains, std::size_t numSteps) :
        SimpleInitializer<FortranCell>(LibGeoDecomp::Coord<2>(), numSteps),
        numDomains(numDomains)
    {
        std::fill(&base_path_buf[0], &base_path_buf[1024], 0);
        std::copy(base_path.begin(), base_path.end(), base_path_buf);

        // Initialize domain decomposition
        std::vector<std::vector<int> > neighbors;

        for(int i = 0; i < numDomains ; i++) {
	/*
            void *domain = NULL;
            FNAME(init_fort)(&i, &domain, base_path_buf);

            FNAME(lgd_yield_subdomain_coord)(&domain, &coord[0]);
            domainCoords << coord;

            neighbors.push_back(neighboringDomainIDs(domain));

            FNAME(term_fort)(&domain);
	*/

            LibGeoDecomp::FloatCoord<2> coord;
	    coord[0] = 0.0;
	    domainCoords << coord;

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
            void *domain = NULL;
            FNAME(init_fort)(&id, &domain, base_path_buf);

            LibGeoDecomp::FloatCoord<2> coord;
            FNAME(lgd_yield_subdomain_coord)(&domain, &coord[0]);
            LibGeoDecomp::Coord<2> logicalCoord = mesher.positionToLogicalCoord(coord);

            if (box.inBounds(logicalCoord)) {
                LOG(INFO, "adding domain " << id
                    <<  " to grid cell " << logicalCoord
                    << " (is at position " << coord << ")\n");
                FortranCell cell = grid->get(logicalCoord);
                cell.insert(id, DomainReference(id, domain));
                grid->set(logicalCoord, cell);
            } else {
                FNAME(term_fort)(&domain);
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
    // std::vector<int> ids;

    // std::vector<int> numneighbors;
    // std::vector<std::vector<int> > neighbors;

    // fixme: make these configurable via command line
    int n_domains = 32;
    int n_timesteps = 100;

    // parse command line
    if (argc != 2) {
        std::cerr << "USAGE: " << argv[0] << " PATH_TO_DECOMPOSE_MESH" << std::endl;
        return 1;
    }

    LibGeoDecomp::SerialSimulator<FortranCell> sim(
        new FortranInitializer(argv[1], n_domains, n_timesteps));
    sim.run();

    // // Get list of neighbors from FORTRAN side of each domain
    // for(int i=0; i<domains.size(); i++) {
    //     int numneighbors_fort;
    //     int neighbors_fort[MAX_DOMAIN_NEIGHBORS];
    //     std::vector<int> neighbors_here;

    //     FNAME(get_neighbors_fort)(&numneighbors_fort,neighbors_fort,&domains[i]);
    //     numneighbors.push_back(numneighbors_fort);

    //     for (int j=0; j<numneighbors_fort; j++) {
    //         std::cout << neighbors_fort[j] << " ";
    //         neighbors_here.push_back(neighbors_fort[j]);
    //     }
    //     std::cout << std::endl;

    //     neighbors.push_back(neighbors_here);

    // }

    // // Print out some information about the domains and their neighbors
    // std::cout << "*** Grid Information ***" << std::endl;
    // std::cout << "domains.size() = " << domains.size() << std::endl;
    // std::cout << "numneighbors.size() = " << numneighbors.size() << std::endl;
    // for (int domain=0; domain<numneighbors.size(); domain++) {
    //     std::cout << "numneighbors[" << domain << "] = " << numneighbors[domain] << std::endl;
    //     std::cout << "neighbors: ";
    //     std::vector<int> neighbors_here = neighbors[domain];
    //     for (int neighbor=0; neighbor<numneighbors[domain]; neighbor++) {
    //         std::cout << neighbors_here[neighbor] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "*** End Grid Information ***" << std::endl;

    // // Print initial state as timestep 0
    // int zero = 0;
    // for(int i=0; i<domains.size(); i++) {
    //     FNAME(print_fort)(&domains[i], &zero);
    // }


    // // Start timestepping loop
    // for(int timestep=1; timestep<=n_timesteps; timestep++) {

    //     // Update values using Game of Life rules
    //     // Loop over all domains
    //     for (int j=0; j<domains.size(); j++) {
    //         FNAME(update_fort)(&timestep,&domains[j]);
    //     } // End loop over domains

    //     // Boundary exchange
    //     // Loop over domains
    //     for (int domain=0; domain<domains.size(); domain++) {
    //         std::vector<int> neighbors_here = neighbors[domain];

    //         //Loop over neighbors
    //         for (int neighbor=0; neighbor<numneighbors[domain]; neighbor++) {
    //             int neighbor_here = neighbors_here[neighbor];
    //             int num_nodes;
    //             int alive[MAX_BOUNDARY_SIZE];

    //             // Get outgoing boundarys from the neighbors
    //             FNAME(get_outgoing_nodes_fort)(&domains[neighbor_here],&domain,&num_nodes,alive);

    //             // Put those arrays inside current domain
    //             FNAME(put_incoming_nodes_fort)(&domains[domain],&neighbor_here,&num_nodes,alive);

    //         }// end loop over neighbors

    //     }// end loop over domains

    //     // Print domains
    //     for(int i=0; i<domains.size(); i++) {
    //         FNAME(print_fort)(&domains[i], &timestep);
    //     }

    // } // end timestep loop

    // // Deallocate domains
    // for(int i=0; i<domains.size(); i++) {
    //     FNAME(term_fort)(&domains[i]);
    // }

    return 0;
}
