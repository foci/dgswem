#include <iostream>
#include <limits>
#include <vector>
#include <exception>
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
#include "dgswem_hpx/domain_reference.hpp"
#include "dgswem_hpx/fortran_initializer.hpp"

//Change me:
typedef LibGeoDecomp::HpxSimulator<DomainReference, LibGeoDecomp::UnstructuredStripingPartition> SimulatorType;

int hpx_main(variables_map & vm)
{
  
  // number of domains read from config file
    int n_domains;
    FNAME(hpx_read_n_domains)(&n_domains);

    int n_timesteps = vm["n_timesteps"].as<std::size_t>();

    std::cout << "n_timesteps from cmd line = " << n_timesteps << std::endl;

    int n_rksteps = 2;

    // Needed by HPX simulator
    std::vector<double> updateGroupSpeeds(1, 1.0);
    int ghostZoneWidth = 1;
    int total_rksteps = n_timesteps*(n_rksteps*2+1)+1;
    std::cout << "total_rksteps = " << total_rksteps << std::endl;

    FortranInitializer *init = new FortranInitializer(n_domains, total_rksteps);

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
    ;

    // We want HPX to run hpx_main() on all localities to avoid the
    // initial overhead caused by broadcasting the work from one to
    // all other localities, and we want to increase the size of the
    // "small" stack.
  std::vector<std::string> config;
  config.push_back("hpx.run_hpx_main!=1");
  config.push_back("hpx.stacks.small_size!=0x20000");

  return hpx::init(desc_commandline,argc,argv,config);
}
