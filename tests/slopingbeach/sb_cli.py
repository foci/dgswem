#!/usr/bin/env python3
"""
CLI wrapper for generating Sloping Beach meshes with WGS84 coordinates.
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

from slopingbeach import SlopingBeach, SlopingBeachMesh, u, write_to_ADCIRC


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s :: %(levelname)s :: %(filename)s :: %(funcName)s :: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%Z",
    )
    log = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Create a simple mesh for testing")
    parser.add_argument("target_elements", type=float, help="Target number of elements")
    parser.add_argument("output_file", type=str, help="Output file name")
    parser.add_argument(
        "--v-shore",
        type=float,
        default=10.0,
        help="Shore bed elevation in meters (default: 10.0)",
    )
    parser.add_argument(
        "--v-offshore",
        type=float,
        default=-200.0,
        help="Offshore bed elevation in meters (default: -200.0)",
    )
    parser.add_argument(
        "--bbox",
        type=float,
        nargs=4,
        default=[-100.0, -60.0, 10.0, 50.0],
        metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"),
        help="WGS84 bounding box coordinates in degrees: lon_min lon_max lat_min lat_max (default: -100 -60 10 50)",
    )
    args = parser.parse_args()

    tick = time.time()

    bbox = (
        (args.bbox[0] * u.deg, args.bbox[1] * u.deg),
        (args.bbox[2] * u.deg, args.bbox[3] * u.deg),
    )
    specs = SlopingBeach(
        v_shore=args.v_shore * u.m,
        v_offshore=args.v_offshore * u.m,
        bbox=bbox,
    )

    mesh = SlopingBeachMesh.from_target_elements(int(args.target_elements), specs=specs)

    log.info(
        "Attempting to create a mesh with dimensions: {} x {}".format(
            mesh.nx, mesh.ny
        )
    )
    nodes, elements, boundaries = mesh.build()
    log.info("Total estimated nodes: {:d}".format(len(nodes)))
    log.info("Total estimated elements: {:d}".format(len(elements)))

    log.info("Writing out the nodes and elements in the ADCIRC format")
    mesh.write_to_ADCIRC(args.output_file)

    tock = time.time()
    log.info("Time to create mesh: {:f} seconds".format(tock - tick))


if __name__ == "__main__":
    main()