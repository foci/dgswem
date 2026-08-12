#!/usr/bin/env python3

TEMPLTE_FORT15_FILE_START = """\
Simple Mesh                              ! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION
                                         ! 24 CHARACTER ALPHANUMERIC RUN IDENTIFICATION
1 20.0 1 50 1000.0                       ! NFOVER - NONFATAL ERROR OVERRIDE OPTION
0                                        ! NABOUT - ABREVIATED OUTPUT OPTION PARAMETER
600                                      ! NSCREEN - OUTPUT TO UNIT 6 PARAMETER
0                                        ! IHOT - HOT START OPTION PARAMETER
2                                        ! ICS - COORDINATE SYSTEM OPTION PARAMETER
0                                        ! IM - MODEL RUN TYPE: 0,10,20,30 = 2DDI, 1,11,21,31 = 3D(VS), 2 = 3D(DSS)
1                                        ! NOLIBF - NONLINEAR BOTTOM FRICTION OPTION
2                                        ! NOLIFA - OPTION TO INCLUDE FINITE AMPLITUDE TERMS
1                                        ! NOLICA - OPTION TO INCLUDE CONVECTIVE ACCELERATION TERMS
1                                        ! NOLICAT - OPTION TO CONSIDER TIME DERIVATIVE OF CONV ACC TERMS
0                                        ! NWP - Number of nodal attributes.
0                                        ! NCOR - VARIABLE CORIOLIS IN SPACE OPTION PARAMETER
0                                        ! NTIP - TIDAL POTENTIAL OPTION PARAMETER
0                                        ! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER
1                                        ! NRAMP - RAMP FUNCTION OPTION
9.81                                     ! G - ACCELERATION DUE TO GRAVITY - DETERMINES UNITS
0.05000                                  ! TAU0 - WEIGHTING FACTOR IN GWCE
{:f}                                     ! DT - TIME STEP (IN SECONDS)
0.00000                                  ! STATIM - STARTING SIMULATION TIME IN DAYS
0.00000                                  ! REFTIME - REFERENCE TIME (IN DAYS) FOR NODAL FACTORS AND EQUILIBRIUM ARGS
5.00000                                  ! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)
1.00000                                  ! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)
0.800000 0.200000 0.000000               ! TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
0.100000 2 10 0.010000                   ! H0, NODEDRYMIN, NODEWETMIN, VELMIN - MINIMUM WATER DEPTH AND DRYING/WETTING OPTIONS
-80.000000 30.000000                     ! SLAM0, SFEA0 - LONGITUDE AND LATITUDE ON WHICH THE CPP COORDINATE PROJECTION IS CENTERED
0.002500                                 ! FFACTOR - 2DDI BOTTOM FRICTION COEFFICIENT
-0.20                                    ! ESLM - SPATIALLY CONSTANT HORIZONTAL EDDY VISCOSITY FOR THE MOMENTUM EQUATIONS
0.000000                                 ! CORI - CONSTANT CORIOLIS COEFFICIENT
0                                        ! NTIF - NUMBER OF TIDAL POTENTIAL CONSTITUENTS
1                                        ! NBFR - NUMBER OF PERIODIC FORCING FREQUENCIES ON ELEVATION SPECIFIED BOUNDARIES
M2                                       ! BOUNTAG - FORCING CONSTITUENT NAME
0.000140520000000 0.96723 8.55
M2                                       ! EALPHA - FORCING CONSTITUENT NAME AGAIN
"""

TEMPLATE_FORT15_FILE_END = """\
110                                      ! ANGINN - MINIMUM ANGLE FOR TANGENTIAL FLOW
0 0.000000 0.000000 0                    ! NOUTE, TOUTSE, TOUTFE, NSPOOLE - FORT 61 OPTIONS
0                                        ! NSTAE - NUMBER OF ELEVATION RECORDING STATIONS, FOLLOWED BY LOCATIONS ON PROCEEDING LINES
0 0.000000 0.000000 0                    ! NOUTV, TOUTSV, TOUTFV, NSPOOLV - FORT 62 OPTIONS
0                                        ! NSTAV - NUMBER OF VELOCITY RECORDING STATIONS, FOLLOWED BY LOCATIONS ON PROCEEDING LINES
1 0.000000 10.000000 {:d}                ! NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE - GLOBAL ELEVATION OUTPUT INFO (UNIT 63)
0 0.000000 0.000000 0                    ! NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV - GLOBAL VELOCITY OUTPUT INFO (UNIT 64)
0                                        ! NHARF - NUMBER OF FREQENCIES IN HARMONIC ANALYSIS
0.000000 0.000000 0 0.000000             ! THAS,THAF,NHAINC,FMV - HARMONIC ANALYSIS PARAMETERS
0 0 0 0                                  ! NHASE,NHASV,NHAGE,NHAGV - CONTROL HARMONIC ANALYSIS AND OUTPUT TO UNITS 51,52,53,54
0 0                                      ! NHSTAR,NHSINC - HOT START FILE GENERATION PARAMETERS
1 0 1e-07 35 0                           ! ITITER, ISLDIA, CONVCR, ITMAX - ALGEBRAIC SOLUTION PARAMETERS
simple_mesh_generator.py
The Water Institute
padcirc
netCDF
one
no_comments
simple_mesh_generator
CF3
zcobell@thewaterinstitute.org
2023-01-01 00:00:00
&wetDryControl slim=0.000400 windlim=True directvelWD=True /
&MetControl DragLawString=garratt WindDragLimit=0.00250 invertedBarometerOnElevationBoundary=true /
"""


def main():
    import numpy as np
    import logging
    import time
    import argparse

    # ...Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s :: %(levelname)s :: %(filename)s :: %(funcName)s :: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%Z",
    )

    log = logging.getLogger(__name__)

    # ...Parse the command line arguments
    parser = argparse.ArgumentParser(description="Create a simple mesh for testing")
    parser.add_argument("target_elements", type=float, help="Target number of elements")
    parser.add_argument("output_file", type=str, help="Output file name")
    args = parser.parse_args()

    tick = time.time()

    target_number_nodes = args.target_elements / 2
    target_res = 40.0 / np.sqrt(target_number_nodes)

    log.info("Target resolution: {:f}".format(target_res))

    x_pos = np.linspace(-100, -60, int((100 - 60) / target_res))
    y_pos = np.linspace(10, 50, int((50 - 10) / target_res))

    # ...Preallocate the nodes array
    nodes = np.zeros((len(x_pos) * len(y_pos), 3), dtype=float)

    v_shore = 10.0
    v_offshore = -200.0

    log.info(
        "Attempting to create a mesh with dimensions: {} x {}".format(
            len(x_pos), len(y_pos)
        )
    )
    log.info("Total estimated nodes: {:d}".format(len(x_pos) * len(y_pos)))
    log.info(
        "Total estimated elements: {:d}".format((len(x_pos) - 1) * (len(y_pos) - 1) * 2)
    )

    # ...Create the nodes. Generate a sloping beach between -200m (south) and 10m (north)
    log.info("Creating the nodes")
    for i, x in enumerate(x_pos):
        for j, y in enumerate(y_pos):
            z = v_offshore + (v_shore - v_offshore) * (y_pos[0] - y) / (
                y_pos[0] - y_pos[-1]
            )
            nodes[i * len(y_pos) + j, :] = [x, y, z]

    # ...Generate boundary conditions. Make a list of the node index.
    # The top row and both sides of nodes are type 20 boundaries
    # ...The bottom row is a type 20 boundary
    log.info("Creating the boundary conditions")
    boundaries = {
        "left": np.arange(0, len(y_pos), 1),
        "bottom": np.arange(0, len(x_pos) * len(y_pos), len(y_pos)),
        "right": np.arange(
            len(x_pos) * len(y_pos) - len(y_pos), len(x_pos) * len(y_pos), 1
        ),
        "top": np.arange(len(y_pos) - 1, len(x_pos) * len(y_pos), len(y_pos)),
    }

    # ...Triangulate the nodes into elements
    log.info("Triangulating the nodes into elements")
    elements = np.zeros(((len(x_pos) - 1) * (len(y_pos) - 1) * 2, 3), dtype=int)
    for i in range(len(x_pos[:-1])):
        for j in range(len(y_pos[:-1])):

            elements[i * (len(y_pos) - 1) * 2 + j * 2, :] = [
                i * len(y_pos) + j,
                (i + 1) * len(y_pos) + j,
                i * len(y_pos) + j + 1,
            ]

            elements[i * (len(y_pos) - 1) * 2 + j * 2 + 1, :] = [
                i * len(y_pos) + j + 1,
                (i + 1) * len(y_pos) + j,
                (i + 1) * len(y_pos) + j + 1,
            ]

    tock = time.time()
    log.info("Time to create mesh: {:f} seconds".format(tock - tick))

    # ...Write out the nodes and elements in the ADCIRC format
    log.info("Writing out the nodes and elements in the ADCIRC format")
    with open(args.output_file, "w") as f:
        f.write("Simple Mesh\n")
        f.write("{:d} {:d}".format(len(elements), len(nodes)))
        for i in range(len(nodes)):
            f.write(
                "\n{: >12d} {:12.8f} {:12.8f} {:12.8f}".format(
                    i + 1, nodes[i, 0], nodes[i, 1], -nodes[i, 2]
                )
            )

        for j in range(len(elements)):
            f.write(
                "\n{: >12d} 3 {:d} {:d} {:d}".format(
                    j + 1, elements[j, 0] + 1, elements[j, 1] + 1, elements[j, 2] + 1
                )
            )

        # ...Write out the boundary conditions
        # ...First, open boundary indexes
        f.write("\n")
        f.write("1\n")
        f.write("{:d}\n".format(len(boundaries["bottom"])))
        f.write("{:d}\n".format(len(boundaries["bottom"])))
        for i in boundaries["bottom"]:
            f.write("{:d}\n".format(i + 1))

        # ...Second, land boundary indexes
        f.write("3\n")
        total_land_boundaries = (
            len(boundaries["top"]) + len(boundaries["left"]) + len(boundaries["right"])
        )
        f.write("{:d}\n".format(total_land_boundaries))

        f.write("{:d} {:d}\n".format(len(boundaries["top"]), 20))
        for i in boundaries["top" "" ""]:
            f.write("{:d}\n".format(i + 1))

        f.write("{:d} {:d}\n".format(len(boundaries["left"]), 20))
        for i in boundaries["left"]:
            f.write("{:d}\n".format(i + 1))

        f.write("{:d} {:d}\n".format(len(boundaries["right"]), 20))
        for i in boundaries["right"]:
            f.write("{:d}\n".format(i + 1))

    # ...Write out the fort.15 file
    time_step = float(int(0.5 * (target_res * 1000.0) / np.sqrt(9.80665)))
    hourly_output = int(3600.0 / time_step)

    log.info("Writing out the fort.15 file")
    with open(args.output_file + ".15", "w") as f:
        f.write(TEMPLTE_FORT15_FILE_START.format(time_step))
        for i in range(len(boundaries["bottom"])):
            f.write("1.000000 0.000                           ! EMO, EFA\n")
        f.write(TEMPLATE_FORT15_FILE_END.format(hourly_output))

    # ...Write some information about the mesh
    with open(args.output_file + ".info", "w") as f:
        f.write("Simple Mesh\n")
        f.write("Number of nodes: {:d}\n".format(len(nodes)))
        f.write("Number of elements: {:d}\n".format(len(elements)))
        f.write("Number of land boundaries: {:d}\n".format(total_land_boundaries))
        f.write("Number of open boundaries: {:d}\n".format(len(boundaries["bottom"])))
        f.write("Target resolution: {:f}\n".format(target_res))
        f.write("Actual number of nodes: {:d}\n".format(len(x_pos) * len(y_pos)))
        f.write(
            "Actual number of elements: {:d}\n".format(
                (len(x_pos) - 1) * (len(y_pos) - 1) * 2
            )
        )
        f.write(
            "Actual number of land boundaries: {:d}\n".format(total_land_boundaries)
        )
        f.write(
            "Actual number of open boundaries: {:d}\n".format(len(boundaries["bottom"]))
        )
        f.write(
            "Actual number of boundary nodes: {:d}\n".format(
                len(boundaries["bottom"])
                + len(boundaries["top"])
                + len(boundaries["left"])
                + len(boundaries["right"])
            )
        )
        f.write(
            "Actual number of boundary elements: {:d}\n".format(
                len(boundaries["bottom"])
                + len(boundaries["top"])
                + len(boundaries["left"])
                + len(boundaries["right"])
            )
        )

    tock = time.time()
    log.info("Total time: {:f} seconds".format(tock - tick))


if __name__ == "__main__":
    main()