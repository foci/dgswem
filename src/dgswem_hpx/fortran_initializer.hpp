#ifndef FORTRAN_INITIALIZER
#define FORTRAN_INITIALIZER

#include "../fname.h"
#include "../fortran_declarations.hpp"

class FortranInitializer : public LibGeoDecomp::SimpleInitializer<DomainReference>
{
public:
    typedef LibGeoDecomp::UnstructuredGrid<DomainReference, MATRICES, double, C, SIGMA> Grid;

  FortranInitializer(std::size_t numDomains, std::size_t numSteps) :
        SimpleInitializer<DomainReference>(LibGeoDecomp::Coord<1>(numDomains), numSteps),
        numDomains(numDomains)
    {
	// Empty
    }

    boost::shared_ptr<LibGeoDecomp::Adjacency>  getAdjacency(const LibGeoDecomp::Region<1>& /* unused */ ) const
    {
	boost::shared_ptr<LibGeoDecomp::Adjacency> adjacency(new LibGeoDecomp::RegionBasedAdjacency());

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
	    //adjacency[id] = neighbors_here;
	    for (int i=0; i < neighbors_here.size(); i++) {
		adjacency->insert(id, neighbors_here[i]);
	    }
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

	    std::cout << "initializing (domain_id = " << id
		      << ", locality = " << hpx::find_here()
		      << ")...\n";
	    
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
	    
	    DomainReference cell(id, size, global, dg, nodalattr, neighbors_here);
	    grid->set(LibGeoDecomp::Coord<1>(id), cell);
        }
	Grid *unstructuredgrid = dynamic_cast<Grid *>(grid);
	unstructuredgrid->setWeights(0, adjacency);
    }

    // Serialization function
    template <class ARCHIVE>
    void serialize(ARCHIVE& ar, unsigned)
    {
	throw std::runtime_error("no serialization yet!");
    }

private:
    std::size_t numDomains;
};

#endif
