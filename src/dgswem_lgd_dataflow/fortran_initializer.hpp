#ifndef FORTRAN_INITIALIZER
#define FORTRAN_INITIALIZER

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "../fname.h"
#include "../fortran_declarations.hpp"

class FortranInitializer : public LibGeoDecomp::SimpleInitializer<DomainReference>
{
public:

  FortranInitializer(std::size_t numDomains, std::size_t numSteps) :
        SimpleInitializer<DomainReference>(LibGeoDecomp::Coord<1>(numDomains), numSteps),
        numDomains(numDomains)
    {
	// Empty
    }

    boost::shared_ptr<LibGeoDecomp::Adjacency>  getAdjacency(const LibGeoDecomp::Region<1>& /* unused */ ) const
    {
	boost::shared_ptr<LibGeoDecomp::Adjacency> adjacency(new LibGeoDecomp::RegionBasedAdjacency());

	hpx::cout << "calling Adjacency" << std::endl;

        for(int id = 0; id < numDomains ; id++) {
	    /*
	    void *size = NULL;
	    void *dg = NULL;
	    void *global = NULL;
	    void *nodalattr = NULL;

	    int domain_number = id;
	    	    
	    // FIXME: Replace initialization with lighter weight function to
	    // retrieve neighbors
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

	    */

	    //Read the DG.18 file directly from here
	    std::ifstream dg18file;
	    std::ostringstream filenamess;
	    filenamess << std::setw(4) << std::setfill('0') << id;
	    std::string filename = "PE"+filenamess.str()+"/DG.18";
	    //std::cout << "filename = " << filename << std::endl;
	    dg18file.open(filename);
	    std::string line;
	    int line_num = 1;
	    int num_recv_neighbors;
	    int num_send_neighbors;
	    std::vector<int> recv_neighbors;
	    std::vector<int> send_neighbors;
	    while (std::getline(dg18file, line)) 
		{
		    std::istringstream iss(line);
		    std::vector<std::string> tokens{
			std::istream_iterator<std::string>(iss),
			    std::istream_iterator<std::string>()
			    };
		    if (line_num==1) {
			if(tokens[0]!="RES") std::cout << "Bad first line in DG.18!\n";
			if(stoi(tokens[2])!=id) std::cout << "id in DG.18 doesn't match!\n";
		    }
		    ++line_num;
		    if (tokens[0]=="COMM") {
			num_recv_neighbors = stoi(tokens[2]);
			num_send_neighbors = stoi(tokens[3]);
		    }
		    if (tokens[0]=="RECV") recv_neighbors.push_back(stoi(tokens[2]));
		    if (tokens[0]=="SEND") send_neighbors.push_back(stoi(tokens[2]));
		}		    
	    /*
	    std::cout << "num_recv_neighbors = " << num_recv_neighbors << std::endl;
	    std::cout << "num_send_neighbors = " << num_send_neighbors << std::endl;
	    */

	    int num_neighbors = std::max(num_recv_neighbors,num_send_neighbors);

	    for (auto&& send_neighbor: send_neighbors) {
		if(std::find(recv_neighbors.begin(), 
			     recv_neighbors.end(), 
			     send_neighbor) != recv_neighbors.end()) {
		    /* already there */
		} else {
		    recv_neighbors.push_back(send_neighbor);
		}
	    }

	    /*

	    // Debugging output
	    std::cout << "From FORTRAN:\n";
	    std::cout << "number of neighbors = " << neighbors_here.size() << std::endl;
	    std::cout << "neighbors:";
	    for (auto&& neighbor: neighbors_here) {
		std::cout << neighbor << " ";
	    }
	    std::cout << std::endl;
	    std::cout << std::endl;
	    std::cout << "From c++ reader:\n";
	    std::cout << "num_neighbors = " << num_neighbors << std::endl;
	    std::cout << "recv_neighbors.size() = " << recv_neighbors.size() << std::endl;
	    std::cout << "neighbors:";
	    for (auto&& neighbor: recv_neighbors) {
		std::cout << neighbor << " ";
	    }
	    std::cout << std::endl;

	    for (int i=0; i < neighbors_here.size(); i++) {
		adjacency->insert(id, neighbors_here[i]);
	    }
	    */
	    for (auto&& neighbor: recv_neighbors) {
		adjacency->insert(id, neighbor);
	    }

	}
	return adjacency;
    }

    void grid(LibGeoDecomp::GridBase<DomainReference, 1> *grid)
    {
	
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
	    
	    DomainReference cell(id, size, global, dg, nodalattr, neighbors_here);
	    grid->set(LibGeoDecomp::Coord<1>(id), cell);
        }
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
