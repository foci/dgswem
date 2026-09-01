#!/usr/bin/env python3
"""
Sloping Beach domain and mesh generation module with WGS84 coordinates and Pint unit enforcement.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pint

# Initialize UnitRegistry
u = pint.UnitRegistry()

# WGS84 authalic / mean Earth radius and physics constants
EARTH_RADIUS = 6371008.8 * u.m
G_CONST = 9.81 * u.meter / u.second**2
CFL_CONST = 0.005

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
zcobell@thewaterinstitute.org, wukenton@gmail.com
2023-01-01 00:00:00
&wetDryControl slim=0.000400 windlim=True directvelWD=True /
&MetControl DragLawString=garratt WindDragLimit=0.00250 invertedBarometerOnElevationBoundary=true /
"""


def enforce_length(val: Any, name: str = "value") -> pint.Quantity:
    """Enforce that a value is a Pint Quantity with length dimensionality."""
    if not isinstance(val, pint.Quantity):
        raise TypeError(
            f"{name} must be a pint.Quantity with length units (e.g. 10 * u.m), got {type(val).__name__}: {val}"
        )
    if not val.check("[length]"):
        raise ValueError(
            f"{name} must have dimensions of [length], got {val.dimensionality} ({val.units})"
        )
    return val


def enforce_angle(val: Any, name: str = "value") -> pint.Quantity:
    """Enforce that a value is a Pint Quantity with degree units (e.g. u.deg)."""
    if not isinstance(val, pint.Quantity):
        raise TypeError(
            f"{name} must be a pint.Quantity with angular units (e.g. 10 * u.deg), got {type(val).__name__}: {val}"
        )
    if not val.units == u.deg:
        raise ValueError(
            f"{name} must have units of u.deg, got {val.units}"
        )
    return val


def parse_wgs84_bbox(
    bbox: tuple[tuple[pint.Quantity, pint.Quantity], tuple[pint.Quantity, pint.Quantity]]
) -> tuple[pint.Quantity, pint.Quantity, pint.Quantity, pint.Quantity]:
    """
    Parse and validate a WGS84 bounding box [[lon_min, lon_max], [lat_min, lat_max]].
    Enforces that each component has angular units, lon in [-180, 180] deg,
    lat in [-90, 90] deg, and min < max.
    """

    lon_range = bbox[0]
    lat_range = bbox[1]

    lon_min = enforce_angle(lon_range[0], "bbox lon_min")
    lon_max = enforce_angle(lon_range[1], "bbox lon_max")
    lat_min = enforce_angle(lat_range[0], "bbox lat_min")
    lat_max = enforce_angle(lat_range[1], "bbox lat_max")
    
    if not (-180.0 * u.deg <= lon_min <= 180.0 * u.deg and -180.0 * u.deg <= lon_max <= 180.0 * u.deg):
        raise ValueError(
            f"Longitude values must be within [-180, 180] degrees, got [{lon_min}, {lon_max}]"
        )
    if not (-90.0 * u.deg <= lat_min <= 90.0 * u.deg and -90.0 * u.deg <= lat_max <= 90.0 * u.deg):
        raise ValueError(
            f"Latitude values must be within [-90, 90] degrees, got [{lat_min_deg}, {lat_max_deg}]"
        )

    if lon_min >= lon_max:
        raise ValueError(
            f"bbox lon_min ({lon_min_deg}) must be strictly less than lon_max ({lon_max_deg})"
        )
    if lat_min >= lat_max:
        raise ValueError(
            f"bbox lat_min ({lat_min_deg}) must be strictly less than lat_max ({lat_max_deg})"
        )

    return lon_min, lon_max, lat_min, lat_max


class SlopingBeach:
    """
    Specification for a sloping beach domain geometry in WGS84 coordinates and bathymetry.

    Parameters:
        v_shore: Bed elevation at the northern/shore boundary (pint Quantity with length units, e.g. u.m).
        v_offshore: Bed elevation at the southern/offshore boundary (pint Quantity with length units, e.g. u.m).
        bbox: WGS84 bounding box [[lon_min, lon_max], [lat_min, lat_max]] with angular units (e.g. u.deg).
    """

    def __init__(
        self,
        v_shore: pint.Quantity = 10.0 * u.m,
        v_offshore: pint.Quantity = -200.0 * u.m,
        bbox: Any = [[-100.0 * u.deg, -60.0 * u.deg], [10.0 * u.deg, 50.0 * u.deg]],
    ):
        self.v_shore = enforce_length(v_shore, "v_shore")
        self.v_offshore = enforce_length(v_offshore, "v_offshore")
        self.lon_min, self.lon_max, self.lat_min, self.lat_max = parse_wgs84_bbox(
            bbox
        )
        self.bbox = [[self.lon_min, self.lon_max], [self.lat_min, self.lat_max]]

    @property
    def dlon(self) -> pint.Quantity:
        """Angular span in longitude (lon_max - lon_min)."""
        return self.lon_max - self.lon_min

    @property
    def dlat(self) -> pint.Quantity:
        """Angular span in latitude (lat_max - lat_min)."""
        return self.lat_max - self.lat_min

    @property
    def mid_lat(self) -> pint.Quantity:
        """Mid-latitude of the domain."""
        return (self.lat_min + self.lat_max) / 2.0

    @property
    def xlen(self) -> pint.Quantity:
        """
        Approximate zonal physical width (length) in meters using equirectangular projection:
        L_x = R * cos(phi_mid) * dlon.
        """
        dlon_rad = self.dlon.to(u.radian).magnitude
        mid_lat_rad = self.mid_lat.to(u.radian).magnitude
        return EARTH_RADIUS * np.cos(mid_lat_rad) * dlon_rad

    @property
    def ylen(self) -> pint.Quantity:
        """
        Approximate meridional physical height (length) in meters using spherical arc length:
        L_y = R * dlat.
        """
        dlat_rad = self.dlat.to(u.radian).magnitude
        return EARTH_RADIUS * dlat_rad

    @property
    def area(self) -> pint.Quantity:
        """Approximate surface area of the domain in square meters."""
        return self.xlen * self.ylen

    @classmethod
    def square(
        cls,
        side: pint.Quantity = 40.0 * u.deg,
        v_shore: pint.Quantity = 10.0 * u.m,
        v_offshore: pint.Quantity = -200.0 * u.m,
        origin: tuple[pint.Quantity, pint.Quantity] = (0.0 * u.deg, 0.0 * u.deg),
    ) -> SlopingBeach:
        """
        Create a square domain defined by an angular span in longitude and latitude.

        Parameters:
            side: Angular side length (pint Quantity with angular units, e.g. 40.0 * u.deg).
            v_shore: Bed elevation at the shore boundary.
            v_offshore: Bed elevation at the offshore boundary.
            origin: Lower-left (lon_min, lat_min) origin with angular units.
        """
        side_angle = enforce_angle(side, "side")
        lon0 = enforce_angle(origin[0], "origin lon0")
        lat0 = enforce_angle(origin[1], "origin lat0")
        bbox = [[lon0, lon0 + side_angle], [lat0, lat0 + side_angle]]
        return cls(v_shore=v_shore, v_offshore=v_offshore, bbox=bbox)


class SlopingBeachMesh:
    """
    Structured triangular mesh generator for a SlopingBeach domain in WGS84 coordinates.
    """

    def __init__(self, specs: SlopingBeach, nx: int, ny: int):
        if not isinstance(specs, SlopingBeach):
            raise TypeError(
                f"specs must be an instance of SlopingBeach, got {type(specs).__name__}"
            )
        if nx < 2 or ny < 2:
            raise ValueError(
                f"nx and ny must both be at least 2, got nx={nx}, ny={ny}"
            )
        self.specs = specs
        self.nx = int(nx)
        self.ny = int(ny)
        self.nodes: np.ndarray | None = None
        self.elements: np.ndarray | None = None
        self.boundaries: dict[str, np.ndarray] | None = None

    @classmethod
    def from_target_elements(
        cls,
        ne: int,
        specs: SlopingBeach | None = None,
    ) -> SlopingBeachMesh:
        """
        Create a domain mesh from a target number of elements.
        If specs is not specified, defaults to SlopingBeach().
        """
        if specs is None:
            specs = SlopingBeach()

        if ne < 2:
            raise ValueError(f"Target elements ne must be at least 2, got {ne}")

        target_number_nodes = ne / 2.0
        target_mesh_size = max(2, int(round(np.sqrt(target_number_nodes))))

        return cls(specs=specs, nx=target_mesh_size, ny=target_mesh_size)

    def create_nodes(self) -> np.ndarray:
        """
        Create node coordinates array [NP x 3] (lon_deg, lat_deg, z_m).
        z slopes linearly from v_offshore at lat_min to v_shore at lat_max.
        """
        lon_min_deg = self.specs.lon_min.to(u.deg).magnitude
        lon_max_deg = self.specs.lon_max.to(u.deg).magnitude
        lat_min_deg = self.specs.lat_min.to(u.deg).magnitude
        lat_max_deg = self.specs.lat_max.to(u.deg).magnitude
        v_offshore_m = self.specs.v_offshore.to(u.m).magnitude
        v_shore_m = self.specs.v_shore.to(u.m).magnitude

        lon_pos = np.linspace(lon_min_deg, lon_max_deg, self.nx)
        lat_pos = np.linspace(lat_min_deg, lat_max_deg, self.ny)

        nodes = np.zeros((self.nx * self.ny, 3), dtype=float)

        for i, lon in enumerate(lon_pos):
            for j, lat in enumerate(lat_pos):
                if self.ny > 1:
                    z = v_offshore_m + (v_shore_m - v_offshore_m) * (
                        lat - lat_min_deg
                    ) / (lat_max_deg - lat_min_deg)
                else:
                    z = v_offshore_m
                nodes[i * self.ny + j, :] = [lon, lat, z]

        self.nodes = nodes
        return nodes

    def create_boundary_conditions(self) -> dict[str, np.ndarray]:
        """
        Generate boundary condition node index arrays (0-based indexing).
        Returns dict with keys: 'left', 'bottom' (open), 'right', 'top'.
        """
        boundaries = {
            "left": np.arange(0, self.ny, 1),
            "bottom": np.arange(0, self.nx * self.ny, self.ny),
            "right": np.arange(
                self.nx * self.ny - self.ny, self.nx * self.ny, 1
            ),
            "top": np.arange(self.ny - 1, self.nx * self.ny, self.ny),
        }
        self.boundaries = boundaries
        return boundaries

    def triangulate_elems(self) -> np.ndarray:
        """
        Triangulate the structured grid into triangular elements [NE x 3] (0-based indexing).
        """
        elements = np.zeros(((self.nx - 1) * (self.ny - 1) * 2, 3), dtype=int)
        for i in range(self.nx - 1):
            for j in range(self.ny - 1):
                idx = i * (self.ny - 1) * 2 + j * 2
                elements[idx, :] = [
                    i * self.ny + j,
                    (i + 1) * self.ny + j,
                    i * self.ny + j + 1,
                ]
                elements[idx + 1, :] = [
                    i * self.ny + j + 1,
                    (i + 1) * self.ny + j,
                    (i + 1) * self.ny + j + 1,
                ]
        self.elements = elements
        return elements

    def build(self) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
        """
        Build and return (nodes, elements, boundaries).
        """
        nodes = self.create_nodes()
        boundaries = self.create_boundary_conditions()
        elements = self.triangulate_elems()
        return nodes, elements, boundaries

    def write_to_ADCIRC(
        self,
        output_file: str | Path,
        title: str = "Simple Mesh",
        write_fort15: bool = True,
        write_info: bool = True,
    ) -> None:
        """
        Write mesh out in ADCIRC fort.14 (and optionally fort.15 and .info) formats.
        """
        write_to_ADCIRC(
            mesh=self,
            output_file=output_file,
            title=title,
            write_fort15=write_fort15,
            write_info=write_info,
        )


def write_to_ADCIRC(
    mesh: SlopingBeachMesh,
    output_file: str | Path,
    title: str = "Simple Mesh",
    write_fort15: bool = True,
    write_info: bool = True,
) -> None:
    """
    Write the SlopingBeachMesh out in ADCIRC fort.14 (and optionally fort.15 and .info) format.
    Coordinates are written in WGS84 decimal degrees (lon, lat, -depth).
    """
    output_file = Path(output_file)

    if mesh.nodes is None:
        mesh.create_nodes()
    if mesh.boundaries is None:
        mesh.create_boundary_conditions()
    if mesh.elements is None:
        mesh.triangulate_elems()

    nodes = mesh.nodes
    elements = mesh.elements
    boundaries = mesh.boundaries

    # 1. Write fort.14 grid file
    with open(output_file, "w") as f:
        f.write(f"{title}\n")
        f.write(f"{len(elements):d} {len(nodes):d}")
        for i in range(len(nodes)):
            f.write(
                f"\n{i + 1:>12d} {nodes[i, 0]:12.8f} {nodes[i, 1]:12.8f} {-nodes[i, 2]:12.8f}"
            )

        for j in range(len(elements)):
            f.write(
                f"\n{j + 1:>12d} 3 {elements[j, 0] + 1:d} {elements[j, 1] + 1:d} {elements[j, 2] + 1:d}"
            )

        # Open boundary (bottom)
        f.write("\n")
        f.write("1\n")
        f.write(f"{len(boundaries['bottom']):d}\n")
        f.write(f"{len(boundaries['bottom']):d}\n")
        for i in boundaries["bottom"]:
            f.write(f"{i + 1:d}\n")

        # Land boundaries (top, left, right)
        f.write("3\n")
        total_land_boundaries = (
            len(boundaries["top"])
            + len(boundaries["left"])
            + len(boundaries["right"])
        )
        f.write(f"{total_land_boundaries:d}\n")

        f.write(f"{len(boundaries['top']):d} 20\n")
        for i in boundaries["top"]:
            f.write(f"{i + 1:d}\n")

        f.write(f"{len(boundaries['left']):d} 20\n")
        for i in boundaries["left"]:
            f.write(f"{i + 1:d}\n")

        f.write(f"{len(boundaries['right']):d} 20\n")
        for i in boundaries["right"]:
            f.write(f"{i + 1:d}\n")

    # 2. Write fort.15 file
    if write_fort15:
        dx = mesh.specs.xlen / (mesh.nx - 1)
        dy = mesh.specs.ylen / (mesh.ny - 1)
        res = min(dx, dy)
        res_angular = min(mesh.specs.dlon / mesh.nx, mesh.specs.dlat / mesh.ny)
        # Implied water depth of 1 m for the CFL condition
        time_step = np.floor(CFL_CONST * res / np.sqrt(G_CONST * u.m))
        hourly_output = int(((1 * u.hour) / time_step).to(u.dimensionless).magnitude)

        fort15_file = str(output_file) + ".15"
        with open(fort15_file, "w") as f:
            f.write(TEMPLTE_FORT15_FILE_START.format(time_step.to(u.s).magnitude))
            for _ in range(len(boundaries["bottom"])):
                f.write("1.000000 0.000                           ! EMO, EFA\n")
            f.write(TEMPLATE_FORT15_FILE_END.format(hourly_output))

    # 3. Write .info file
    if write_info:
        info_file = str(output_file) + ".info"
        dx = mesh.specs.xlen / (mesh.nx - 1)
        dy = mesh.specs.ylen / (mesh.ny - 1)
        res = min(dx, dy)
        res_angular = min(mesh.specs.dlon / mesh.nx, mesh.specs.dlat / mesh.ny)
        total_land_boundaries = (
            len(boundaries["top"])
            + len(boundaries["left"])
            + len(boundaries["right"])
        )
        with open(info_file, "w") as f:
            f.write(f"{title}\n")
            f.write(f"Number of nodes: {len(nodes):d}\n")
            f.write(f"Number of elements: {len(elements):d}\n")
            f.write(f"Number of land boundaries: {total_land_boundaries:d}\n")
            f.write(f"Number of open boundaries: {len(boundaries['bottom']):d}\n")
            f.write(f"Target metric resolution: {res.to(u.m).magnitude:f} m\n")
            f.write(f"Target angular resolution: {res_angular.to(u.deg).magnitude:f} deg\n")
            f.write(f"Actual number of nodes: {len(nodes):d}\n")
            f.write(f"Actual number of elements: {len(elements):d}\n")
            f.write(f"Actual number of land boundaries: {total_land_boundaries:d}\n")
            f.write(f"Actual number of open boundaries: {len(boundaries['bottom']):d}\n")
            f.write(
                f"Actual number of boundary nodes: {len(boundaries['bottom']) + total_land_boundaries:d}\n"
            )
            f.write(
                f"Actual number of boundary elements: {len(boundaries['bottom']) + total_land_boundaries:d}\n"
            )