#!/usr/bin/env python
"""
Build the Star Destroyer fleet demo: a detailed Imperial Star Destroyer, a squadron
of TIE fighters and a flight of X-wings, where every ship is ONE PIMMS chain.

Each ship is first drawn as a set of lattice sites (voxels), each tagged with a
bead type so that the parts can be coloured separately. The sites are then
threaded into a single self-avoiding chain that visits every one of them exactly
once, with each bead bonded to the next (a bond joins two sites that differ by at
most 1 on every axis, as in PIMMS). Because the whole ship is one chain:

* the Star Destroyer is frozen with a single ``C 1`` line in ``freeze.txt``, and
* each fighter keeps its exact shape when the run uses only rigid-body moves
  (whole-chain rotations and VMMC translations), so the fighters fly around as
  rigid bodies.

Run from this directory:   ``python build_fleet.py``
Outputs:                   ``fleet.restart``, ``freeze.txt`` and ``params.prm``

The geometry is deterministic: the same ``BUILD_SEED`` always gives the same files.
"""

from __future__ import annotations

import itertools
import math
import pickle
import random
from typing import Callable, Dict, Iterable, List, Sequence, Set, Tuple

Site = Tuple[int, int, int]
Model = Dict[Site, str]
Rotation = Tuple[Tuple[int, int, int], Tuple[int, int, int], Tuple[int, int, int]]

# ---------------------------------------------------------------------------
# output files and the box
# ---------------------------------------------------------------------------
RESTART_FILENAME = "fleet.restart"
FREEZE_FILENAME = "freeze.txt"
PARAMETER_FILENAME = "params.prm"

BOX_DIMENSIONS: Site = (210, 120, 96)      # DIMENSIONS in KEYFILE.kf must be at least this
BUILD_SEED = 1977                          # seeds the greebles, the threading and the fighter placement

# ---------------------------------------------------------------------------
# bead types (single letters; PIMMS writes each as a PDB residue name, so the
# parts can be coloured by residue name in VMD / PyMOL - see readme.md)
# ---------------------------------------------------------------------------
SD_HULL = "H"             # HIS   Star Destroyer hull
SD_PANEL = "G"            # GLY   darker hull panels and surface greebles
SD_SUPERSTRUCTURE = "T"   # THR   command superstructure, tower and bridge
SD_WINDOW = "W"           # TRP   lit windows (bridge and superstructure)
SD_DOME = "D"             # ASP   shield-generator domes
SD_ENGINE = "K"           # LYS   engine housings
SD_ENGINE_GLOW = "E"      # GLU   engine glow
SD_HANGAR_LIGHT = "L"     # LEU   ventral hangar-bay lights
SD_TURRET = "A"           # ALA   turbolaser turrets

TIE_COCKPIT = "C"         # CYS   TIE cockpit ball
TIE_VIEWPORT = "V"        # VAL   TIE forward viewport
TIE_FRAME = "F"           # PHE   TIE wing frame, spars and hub
TIE_PANEL = "P"           # PRO   TIE solar panels
TIE_PYLON = "N"           # ASN   TIE wing pylons

XWING_HULL = "S"          # SER   X-wing fuselage and wings
XWING_STRIPE = "R"        # ARG   X-wing red squadron stripes
XWING_CANOPY = "Q"        # GLN   X-wing cockpit canopy
XWING_DROID = "X"         # XXX   X-wing astromech droid
XWING_CANNON = "I"        # ILE   X-wing wingtip laser cannons
XWING_ENGINE = "M"        # MET   X-wing engines
XWING_GLOW = "Y"          # TYR   X-wing engine glow

ALL_BEAD_TYPES: Tuple[str, ...] = (
    SD_HULL, SD_PANEL, SD_SUPERSTRUCTURE, SD_WINDOW, SD_DOME, SD_ENGINE, SD_ENGINE_GLOW,
    SD_HANGAR_LIGHT, SD_TURRET,
    TIE_COCKPIT, TIE_VIEWPORT, TIE_FRAME, TIE_PANEL, TIE_PYLON,
    XWING_HULL, XWING_STRIPE, XWING_CANOPY, XWING_DROID, XWING_CANNON, XWING_ENGINE, XWING_GLOW,
)

# ---------------------------------------------------------------------------
# Star Destroyer geometry, in lattice sites. Local frame: x runs from the bow
# (x = 0) to the stern, y across the ship (0 on the centreline), z up (0 is the
# plane of the equatorial trench that separates the upper and lower hulls).
# The real ship is about 1600 m long and 1000 m wide; the beam ratio is kept.
# ---------------------------------------------------------------------------
SD_LENGTH = 128                 # bow to stern face
SD_HALF_BEAM = 39               # half-width of the hull at the stern
SD_UPPER_RISE = 10              # how far the dorsal ridge rises above the hull edge, at the stern
SD_LOWER_DROP = 7               # how far the ventral keel drops below the hull edge, at the stern
SD_EDGE_UP = 2                  # thickness of the hull edge above the trench, at the stern
SD_EDGE_DOWN = 2                # thickness of the hull edge below the trench, at the stern
SD_TRENCH_DEPTH = 2             # how far the equatorial trench cuts in from the hull edge
SD_TRENCH_START = 0.12          # the trench starts this fraction of the length behind the bow
SD_SPINE_HALF_WIDTH = 2         # half-width of the raised dorsal spine
SD_SPINE_START = 0.25           # the spine starts this fraction of the length behind the bow

# the stepped command superstructure: (front, back) as fractions of the length,
# half-width and the height of its flat top above the trench plane
SD_TIERS: Tuple[Tuple[float, float, int, int], ...] = (
    (0.56, 0.96, 14, 14),
    (0.64, 0.95, 11, 17),
    (0.71, 0.94, 8, 19),
    (0.77, 0.93, 5, 21),
)
SD_TIER_WINDOW_SPACING = 2      # lit windows every this many sites along the superstructure flanks

# the tower neck, the bridge on top of it, and the two shield-generator domes
SD_NECK_X: Tuple[float, float] = (0.835, 0.865)    # front, back as fractions of the length
SD_NECK_HALF_WIDTH = 2
SD_NECK_TOP = 29
SD_BRIDGE_X: Tuple[float, float] = (0.825, 0.875)
SD_BRIDGE_HALF_WIDTH = 11
SD_BRIDGE_Z: Tuple[int, int] = (30, 32)
SD_DOME_Y = 8                   # the domes sit at y = +/- this
SD_DOME_RADIUS = 2.6
SD_DOME_CENTRE_Z = 36
SD_DOME_STALK_Z: Tuple[int, int] = (33, 34)

# ventral hangar bay (a recess in the keel, lit inside)
SD_HANGAR_X: Tuple[float, float] = (0.42, 0.60)
SD_HANGAR_HALF_WIDTH = 5
SD_HANGAR_DEPTH = 3

# engines on the stern face: (y, z, radius, length behind the stern)
SD_MAIN_ENGINES: Tuple[Tuple[int, int, float, int], ...] = (
    (0, 2, 4.6, 5), (-13, 2, 4.6, 5), (13, 2, 4.6, 5),
)
SD_SMALL_ENGINES: Tuple[Tuple[int, int, float, int], ...] = (
    (-6, -4, 2.1, 4), (6, -4, 2.1, 4), (-21, -1, 2.1, 4), (21, -1, 2.1, 4),
)
SD_NOZZLE_RING_FRACTION = 0.62  # glow fills the nozzle inside this fraction of the engine radius

# turbolaser turrets: on the superstructure flanks and along the dorsal hull edges
SD_FLANK_TURRET_X: Tuple[float, float] = (0.58, 0.92)
SD_FLANK_TURRET_SPACING = 6
SD_EDGE_TURRET_X: Tuple[float, float] = (0.30, 0.86)
SD_EDGE_TURRET_SPACING = 8
SD_EDGE_TURRET_INSET = 5        # turrets sit this many sites in from the hull edge

# surface detail
SD_N_GREEBLES = 260             # small raised blocks scattered over the dorsal hull
SD_GREEBLE_MAX_SIZE = 4         # largest greeble footprint edge
SD_N_PANELS = 220               # darker rectangular hull panels (paint only)
SD_PANEL_MAX_SIZE = 7
SD_DETAIL_START = 0.10          # no surface detail ahead of this fraction of the length

# The hull is hollow: only a skin this many sites deep is kept, which takes the
# chain from about 66,000 beads (solid) to about 44,000. Every site inside is at
# least SD_SHELL_THICKNESS + 1 from open space, so a fighter can never reach the
# hollow by a VMMC step of up to VMMC_MAX_DISPLACEMENT = 2 sites per axis, and a
# bonded chain can't straddle the skin. (A whole-box teleport, MOVE_CHAIN_TRANSLATE
# or MOVE_CLUSTER_TRANSLATE, could drop a fighter inside - so the demo leaves them off.)
SD_SHELL_THICKNESS = 3

# ---------------------------------------------------------------------------
# TIE fighter geometry. Local frame: x forward, y along the wing axis, z up.
# (Fighters are drawn far larger than true scale, or they'd be under a site.)
# ---------------------------------------------------------------------------
TIE_BALL_RADIUS = 3.5
TIE_VIEWPORT_DEPTH = 2.5        # viewport voxels are those with x >= this
TIE_PYLON_HALF_WIDTH = 1
TIE_WING_Y = 9                  # the wings are the planes y = +/- this
TIE_WING_RADIUS = 8.0           # centre-to-vertex radius of the hexagonal wings
TIE_HUB_HALF_WIDTH = 1
TIE_SPAR_HALF_WIDTH = 0.5

# ---------------------------------------------------------------------------
# X-wing geometry. Local frame: x forward (nose at +x), y along the span, z up.
# ---------------------------------------------------------------------------
XW_NOSE_X: Tuple[int, int] = (2, 14)          # nose root .. tip
XW_NOSE_ROOT_HALF_WIDTH = 2
XW_COCKPIT_X: Tuple[int, int] = (-4, 1)
XW_COCKPIT_HALF_WIDTH = 2
XW_FUSELAGE_Z: Tuple[int, int] = (-2, 2)
XW_CANOPY_X: Tuple[int, int] = (-3, 1)
XW_DROID_X: Tuple[int, int] = (-7, -6)
XW_REAR_X: Tuple[int, int] = (-13, -5)
XW_REAR_HALF_WIDTH = 3
XW_WING_Y: Tuple[int, int] = (4, 13)           # wing root .. tip, |y|
XW_WING_SLOPE = 0.5                           # rise of each S-foil per site of span
XW_WING_ROOT_X: Tuple[int, int] = (-13, -6)
XW_WING_TIP_X: Tuple[int, int] = (-12, -9)
XW_STRIPE_Y: Tuple[int, int] = (8, 10)         # the red stripe band, |y|
XW_ENGINE_Y = 5
XW_ENGINE_Z = 3
XW_ENGINE_X: Tuple[int, int] = (-15, -5)
XW_CANNON_X: Tuple[int, int] = (-12, 5)

# ---------------------------------------------------------------------------
# the fleet and where it starts (box coordinates)
# ---------------------------------------------------------------------------
SD_BOX_ORIGIN: Site = (60, 60, 36)            # where the Star Destroyer's local (0, 0, 0) goes
N_TIE_FIGHTERS = 12
N_X_WINGS = 8
# regions (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi) that each fighter's centre is drawn from
TIE_REGION: Tuple[int, int, int, int, int, int] = (66, 160, 12, 108, 54, 84)
XWING_REGION: Tuple[int, int, int, int, int, int] = (18, 40, 16, 104, 12, 84)
MAX_PLACEMENT_TRIES = 20000                   # per fighter, before giving up loudly

# the threader: how many Posa rotations it may try in all, and in a row without
# finding a way forward, before it stops and splices in what's left
MAX_ROTATIONS = 200000
MAX_FRUITLESS_ROTATIONS = 2000

# the 26 bond vectors, face steps first so ties favour tidy, axis-aligned bonds
BOND_STEPS: Tuple[Site, ...] = tuple(sorted(
    ((dx, dy, dz) for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3) if (dx, dy, dz) != (0, 0, 0)),
    key=lambda step: (abs(step[0]) + abs(step[1]) + abs(step[2]), step)))


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def round_half_up(value: float) -> int:
    """
    Round to the nearest integer, halves away from zero.

    Python's ``round`` sends halves to the nearest even integer, which makes
    sloped surfaces uneven; this keeps every step of a slope the same.

    Parameters
    ----------
    value : float
        The number to round.

    Returns
    -------
    int
        The nearest integer, with x.5 rounded away from zero.
    """
    return int(math.floor(abs(value) + 0.5)) * (1 if value >= 0 else -1)


def add_sites(site: Site, step: Site) -> Site:
    """
    Add two lattice vectors.

    Parameters
    ----------
    site : tuple of int
        A lattice site.
    step : tuple of int
        A displacement.

    Returns
    -------
    tuple of int
        ``site + step``.
    """
    return (site[0] + step[0], site[1] + step[1], site[2] + step[2])


def are_bonded(site_a: Site, site_b: Site) -> bool:
    """
    Whether two distinct sites are close enough to be bonded in PIMMS.

    Parameters
    ----------
    site_a, site_b : tuple of int
        Two lattice sites.

    Returns
    -------
    bool
        True if they differ by at most 1 on every axis and are not the same site.
    """
    gaps = [abs(a - b) for a, b in zip(site_a, site_b)]
    return max(gaps) == 1


def fill_box(model: Model, x_range: Tuple[int, int], y_range: Tuple[int, int],
             z_range: Tuple[int, int], bead: str) -> None:
    """
    Set every site of an inclusive box to one bead type.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    x_range, y_range, z_range : tuple of int
        Inclusive (low, high) bounds on each axis.
    bead : str
        The bead type to write.

    Returns
    -------
    None
    """
    for x in range(x_range[0], x_range[1] + 1):
        for y in range(y_range[0], y_range[1] + 1):
            for z in range(z_range[0], z_range[1] + 1):
                model[(x, y, z)] = bead


def fill_sphere(model: Model, centre: Tuple[float, float, float], radius: float, bead: str) -> None:
    """
    Set every site within ``radius`` of ``centre`` to one bead type.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    centre : tuple of float
        Centre of the sphere, in sites.
    radius : float
        Radius of the sphere, in sites.
    bead : str
        The bead type to write.

    Returns
    -------
    None
    """
    reach = int(math.ceil(radius))
    cx, cy, cz = centre
    for x in range(int(math.floor(cx)) - reach, int(math.ceil(cx)) + reach + 1):
        for y in range(int(math.floor(cy)) - reach, int(math.ceil(cy)) + reach + 1):
            for z in range(int(math.floor(cz)) - reach, int(math.ceil(cz)) + reach + 1):
                if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= radius ** 2:
                    model[(x, y, z)] = bead


# ---------------------------------------------------------------------------
# the Star Destroyer
# ---------------------------------------------------------------------------
def sd_half_width(x: int) -> int:
    """
    Half-width of the Star Destroyer's hull at length position ``x``.

    The hull is the classic arrowhead seen from above: zero at the bow, growing
    linearly to ``SD_HALF_BEAM`` at the stern.

    Parameters
    ----------
    x : int
        Position along the ship, 0 at the bow.

    Returns
    -------
    int
        The half-width in sites.
    """
    return round_half_up(SD_HALF_BEAM * x / (SD_LENGTH - 1))


def sd_top(x: int, y: int) -> int:
    """
    Height of the dorsal (upper) hull surface at ``(x, y)``.

    The upper hull is a low pyramid: a thin edge wall that rises in a straight
    slope to a ridge on the centreline, all of it growing toward the stern,
    with a raised spine along the centreline behind ``SD_SPINE_START``.

    Parameters
    ----------
    x : int
        Position along the ship, 0 at the bow.
    y : int
        Position across the ship, 0 on the centreline; must be within the hull.

    Returns
    -------
    int
        The z of the highest hull site in this column.
    """
    fraction = x / (SD_LENGTH - 1)
    half_width = max(sd_half_width(x), 1)
    edge = max(1, round_half_up(SD_EDGE_UP * fraction))
    ridge = SD_UPPER_RISE * fraction * (1 - abs(y) / half_width)
    spine = 1 if (fraction >= SD_SPINE_START and abs(y) <= SD_SPINE_HALF_WIDTH) else 0
    return edge + round_half_up(ridge) + spine


def sd_bottom(x: int, y: int) -> int:
    """
    Depth of the ventral (lower) hull surface at ``(x, y)``.

    The mirror of :func:`sd_top`: a thin edge wall below the trench that drops
    in a straight slope to a keel on the centreline.

    Parameters
    ----------
    x : int
        Position along the ship, 0 at the bow.
    y : int
        Position across the ship, 0 on the centreline; must be within the hull.

    Returns
    -------
    int
        The z of the lowest hull site in this column (negative).
    """
    fraction = x / (SD_LENGTH - 1)
    half_width = max(sd_half_width(x), 1)
    edge = max(1, round_half_up(SD_EDGE_DOWN * fraction))
    keel = SD_LOWER_DROP * fraction * (1 - abs(y) / half_width)
    return -(edge + round_half_up(keel))


def build_sd_hull(model: Model) -> None:
    """
    Add the upper and lower hulls, with the equatorial trench between them.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    trench_start = SD_TRENCH_START * SD_LENGTH
    for x in range(SD_LENGTH):
        half_width = sd_half_width(x)
        for y in range(-half_width, half_width + 1):
            in_trench_band = x >= trench_start and abs(y) > half_width - SD_TRENCH_DEPTH
            for z in range(sd_bottom(x, y), sd_top(x, y) + 1):
                if z == 0 and in_trench_band:
                    continue
                model[(x, y, z)] = SD_HULL


def build_sd_hangar(model: Model) -> None:
    """
    Cut the ventral hangar bay into the keel and light its ceiling.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    x_lo = int(SD_HANGAR_X[0] * SD_LENGTH)
    x_hi = int(SD_HANGAR_X[1] * SD_LENGTH)
    for x in range(x_lo, x_hi + 1):
        for y in range(-SD_HANGAR_HALF_WIDTH, SD_HANGAR_HALF_WIDTH + 1):
            floor = sd_bottom(x, y)
            if floor + SD_HANGAR_DEPTH > -1:
                continue                          # too shallow here: the bay would cut through the hull
            for z in range(floor, floor + SD_HANGAR_DEPTH):
                model.pop((x, y, z), None)
            model[(x, y, floor + SD_HANGAR_DEPTH)] = SD_HANGAR_LIGHT


def build_sd_superstructure(model: Model) -> None:
    """
    Add the stepped command superstructure, with rows of lit windows.

    Each tier is a flat-topped block that fills everything from the hull up to
    its top; later tiers are narrower, shorter and taller, giving the stepped
    ziggurat behind the bridge tower.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    for front, back, half_width, top in SD_TIERS:
        x_lo = int(front * SD_LENGTH)
        x_hi = min(int(back * SD_LENGTH), SD_LENGTH - 1)
        for x in range(x_lo, x_hi + 1):
            for y in range(-half_width, half_width + 1):
                for z in range(1, top + 1):
                    model[(x, y, z)] = SD_SUPERSTRUCTURE
            # a row of lit windows just under the top edge of each flank
            if (x - x_lo) % SD_TIER_WINDOW_SPACING == 0:
                model[(x, -half_width, top - 1)] = SD_WINDOW
                model[(x, half_width, top - 1)] = SD_WINDOW


def build_sd_bridge(model: Model) -> None:
    """
    Add the tower neck, the bridge on top of it and the two shield domes.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    top_tier = SD_TIERS[-1][3]
    neck_x = (int(SD_NECK_X[0] * SD_LENGTH), int(SD_NECK_X[1] * SD_LENGTH))
    bridge_x = (int(SD_BRIDGE_X[0] * SD_LENGTH), int(SD_BRIDGE_X[1] * SD_LENGTH))
    fill_box(model, neck_x, (-SD_NECK_HALF_WIDTH, SD_NECK_HALF_WIDTH), (top_tier + 1, SD_NECK_TOP),
             SD_SUPERSTRUCTURE)
    fill_box(model, bridge_x, (-SD_BRIDGE_HALF_WIDTH, SD_BRIDGE_HALF_WIDTH), SD_BRIDGE_Z, SD_SUPERSTRUCTURE)
    # the bridge windows run along its forward face
    window_z = (SD_BRIDGE_Z[0] + SD_BRIDGE_Z[1]) // 2
    for y in range(-SD_BRIDGE_HALF_WIDTH + 1, SD_BRIDGE_HALF_WIDTH):
        model[(bridge_x[0], y, window_z)] = SD_WINDOW
    # each dome sits on a 2x2 stalk on the bridge roof
    dome_x = (bridge_x[0] + bridge_x[1]) / 2
    for side in (-1, 1):
        dome_y = side * SD_DOME_Y
        stalk_x = (int(math.floor(dome_x)), int(math.floor(dome_x)) + 1)
        stalk_y = (dome_y, dome_y + side) if side > 0 else (dome_y + side, dome_y)
        fill_box(model, stalk_x, (min(stalk_y), max(stalk_y)), SD_DOME_STALK_Z, SD_SUPERSTRUCTURE)
        fill_sphere(model, (dome_x, dome_y + side * 0.5, SD_DOME_CENTRE_Z), SD_DOME_RADIUS, SD_DOME)


def add_engine(model: Model, y_centre: int, z_centre: int, radius: float, length: int,
               stern_x: int) -> None:
    """
    Add one engine: a cylinder out of the stern, with a recessed glowing nozzle.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    y_centre, z_centre : int
        The engine's axis, across and up.
    radius : float
        Engine radius in sites.
    length : int
        How far the engine sticks out behind the stern face.
    stern_x : int
        x of the stern face; the engine starts at ``stern_x + 1``.

    Returns
    -------
    None
    """
    reach = int(math.ceil(radius))
    glow_radius = radius * SD_NOZZLE_RING_FRACTION
    for x in range(stern_x + 1, stern_x + length + 1):
        is_glow_layer = x == stern_x + length - 1
        is_rim_layer = x == stern_x + length
        for dy in range(-reach, reach + 1):
            for dz in range(-reach, reach + 1):
                distance = math.hypot(dy, dz)
                if distance > radius:
                    continue
                site = (x, y_centre + dy, z_centre + dz)
                if is_rim_layer and distance < glow_radius:
                    continue                      # the open mouth of the nozzle
                if is_glow_layer and distance < glow_radius:
                    model[site] = SD_ENGINE_GLOW
                else:
                    model[site] = SD_ENGINE


def build_sd_engines(model: Model) -> None:
    """
    Add the three main engines and four smaller ones on the stern face.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    for y_centre, z_centre, radius, length in SD_MAIN_ENGINES + SD_SMALL_ENGINES:
        add_engine(model, y_centre, z_centre, radius, length, SD_LENGTH - 1)


def add_turret(model: Model, x: int, y: int) -> None:
    """
    Add one turbolaser turret on the hull surface at ``(x, y)``.

    A turret is a level 3x3 base (built up from the hull, which slopes) with a
    raised gun housing on top. A turret that would land on another one is
    skipped.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    x, y : int
        Centre of the turret's footprint.

    Returns
    -------
    None
    """
    footprint = [(x + dx, y + dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]
    column_tops = {column: top_of_column(model, column[0], column[1]) for column in footprint}
    if any(model[(cx, cy, top)] == SD_TURRET for (cx, cy), top in column_tops.items()):
        return
    base = max(column_tops.values()) + 1
    for (cx, cy), top in column_tops.items():
        for z in range(top + 1, base + 1):
            model[(cx, cy, z)] = SD_TURRET
    model[(x, y, base + 1)] = SD_TURRET
    model[(x + 1, y, base + 1)] = SD_TURRET


def top_of_column(model: Model, x: int, y: int) -> int:
    """
    The highest occupied z in the column at ``(x, y)`` above the trench plane.

    Parameters
    ----------
    model : dict
        Site -> bead type.
    x, y : int
        The column.

    Returns
    -------
    int
        The top occupied z; :func:`sd_top` if nothing above the hull has been added.

    Raises
    ------
    ValueError
        If the hull surface site of the column is missing (the caller asked for
        a column off the hull).
    """
    z = sd_top(x, y)
    if (x, y, z) not in model:
        raise ValueError("column (%d, %d) is not on the hull" % (x, y))
    while (x, y, z + 1) in model:
        z += 1
    return z


def build_sd_turrets(model: Model) -> None:
    """
    Add turret rows along the superstructure flanks and the dorsal hull edges.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.

    Returns
    -------
    None
    """
    flank_half_width = SD_TIERS[0][2] + 3
    for x in range(int(SD_FLANK_TURRET_X[0] * SD_LENGTH), int(SD_FLANK_TURRET_X[1] * SD_LENGTH),
                   SD_FLANK_TURRET_SPACING):
        for side in (-1, 1):
            add_turret(model, x, side * flank_half_width)
    for x in range(int(SD_EDGE_TURRET_X[0] * SD_LENGTH), int(SD_EDGE_TURRET_X[1] * SD_LENGTH),
                   SD_EDGE_TURRET_SPACING):
        for side in (-1, 1):
            add_turret(model, x, side * (sd_half_width(x) - SD_EDGE_TURRET_INSET))


def open_dorsal_cells(model: Model, rng: random.Random) -> Tuple[int, int]:
    """
    Draw a random column on the open dorsal hull (clear of the superstructure).

    Parameters
    ----------
    model : dict
        Site -> bead type.
    rng : random.Random
        Source of randomness.

    Returns
    -------
    tuple of int
        An ``(x, y)`` column whose top site is plain hull.
    """
    while True:
        x = rng.randrange(int(SD_DETAIL_START * SD_LENGTH), SD_LENGTH - 1)
        half_width = sd_half_width(x) - SD_TRENCH_DEPTH - 1
        if half_width < 1:
            continue
        y = rng.randrange(-half_width, half_width + 1)
        if model.get((x, y, top_of_column(model, x, y))) == SD_HULL:
            return (x, y)


def build_sd_surface_detail(model: Model, rng: random.Random) -> None:
    """
    Scatter darker hull panels and small raised greebles over the dorsal hull.

    Panels only repaint surface sites; greebles add a one-site-high block that
    follows the hull surface, so each of their sites sits on the hull.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    rng : random.Random
        Source of randomness (seeded, so the ship is the same every build).

    Returns
    -------
    None
    """
    for _ in range(SD_N_PANELS):
        x0, y0 = open_dorsal_cells(model, rng)
        length = rng.randint(2, SD_PANEL_MAX_SIZE)
        width = rng.randint(1, SD_PANEL_MAX_SIZE // 2)
        for x in range(x0, min(x0 + length, SD_LENGTH - 1)):
            for y in range(y0, y0 + width):
                if abs(y) > sd_half_width(x) - SD_TRENCH_DEPTH - 1:
                    continue
                top = top_of_column(model, x, y)
                if model[(x, y, top)] == SD_HULL:
                    model[(x, y, top)] = SD_PANEL
    for _ in range(SD_N_GREEBLES):
        x0, y0 = open_dorsal_cells(model, rng)
        length = rng.randint(1, SD_GREEBLE_MAX_SIZE)
        width = rng.randint(1, SD_GREEBLE_MAX_SIZE)
        for x in range(x0, min(x0 + length, SD_LENGTH - 1)):
            for y in range(y0, y0 + width):
                if abs(y) > sd_half_width(x) - SD_TRENCH_DEPTH - 1:
                    continue
                top = top_of_column(model, x, y)
                if model[(x, y, top)] in (SD_HULL, SD_PANEL):
                    model[(x, y, top + 1)] = SD_PANEL


def star_destroyer_model(rng: random.Random) -> Model:
    """
    Draw the Imperial Star Destroyer as tagged lattice sites.

    Parameters
    ----------
    rng : random.Random
        Source of randomness for the surface detail.

    Returns
    -------
    dict
        Site (local frame, bow at x = 0) -> bead type.
    """
    model: Model = {}
    build_sd_hull(model)
    build_sd_hangar(model)
    build_sd_superstructure(model)
    build_sd_bridge(model)
    build_sd_engines(model)
    build_sd_turrets(model)
    build_sd_surface_detail(model, rng)
    hollow_out(model, SD_SHELL_THICKNESS)
    return model


def hollow_out(model: Model, thickness: int) -> None:
    """
    Remove every site deeper than ``thickness`` below the surface.

    A site's depth is its Chebyshev distance to the nearest empty site, so the
    sites kept form a skin ``thickness`` sites deep everywhere, and every site
    removed is at least ``thickness + 1`` from open space.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    thickness : int
        Depth of the skin to keep, in sites.

    Returns
    -------
    None
    """
    sites = set(model)
    depth = {site: 1 for site in sites
             if any(add_sites(site, step) not in sites for step in BOND_STEPS)}
    layer = list(depth)
    for current in range(2, thickness + 1):
        next_layer = []
        for site in layer:
            for step in BOND_STEPS:
                neighbour = add_sites(site, step)
                if neighbour in sites and neighbour not in depth:
                    depth[neighbour] = current
                    next_layer.append(neighbour)
        layer = next_layer
    for site in sites:
        if site not in depth:
            del model[site]


# ---------------------------------------------------------------------------
# the TIE fighter
# ---------------------------------------------------------------------------
def in_tie_wing(x: int, z: int) -> bool:
    """
    Whether ``(x, z)`` lies inside a TIE wing's hexagon.

    The hexagon has a vertex at the top and the bottom, as on a TIE/ln.

    Parameters
    ----------
    x, z : int
        Position in the plane of the wing, relative to its centre.

    Returns
    -------
    bool
        True inside the hexagon.
    """
    half_span = math.sqrt(3) / 2 * TIE_WING_RADIUS
    if abs(x) > half_span + 0.25:
        return False
    return abs(z) <= TIE_WING_RADIUS - abs(x) / math.sqrt(3) + 0.5


def tie_wing_bead(x: int, z: int) -> str:
    """
    The bead type of a wing site: frame (rim, spars, hub) or solar panel.

    Parameters
    ----------
    x, z : int
        Position in the plane of the wing, relative to its centre.

    Returns
    -------
    str
        ``TIE_FRAME`` or ``TIE_PANEL``.
    """
    on_rim = any(not in_tie_wing(x + dx, z + dz) for dx, dz in ((1, 0), (-1, 0), (0, 1), (0, -1)))
    on_hub = abs(x) <= TIE_HUB_HALF_WIDTH and abs(z) <= TIE_HUB_HALF_WIDTH
    on_vertical_spar = abs(x) <= TIE_SPAR_HALF_WIDTH
    on_diagonal_spar = any(abs(z - sign * x / math.sqrt(3)) <= TIE_SPAR_HALF_WIDTH for sign in (-1, 1))
    if on_rim or on_hub or on_vertical_spar or on_diagonal_spar:
        return TIE_FRAME
    return TIE_PANEL


def tie_fighter_model() -> Model:
    """
    Draw a TIE fighter as tagged lattice sites.

    A ball cockpit with a forward viewport, two square pylons, and two flat
    hexagonal wings with a frame of rim, hub and six spars.

    Returns
    -------
    dict
        Site (local frame, cockpit centred on the origin, facing +x) -> bead type.
    """
    model: Model = {}
    fill_sphere(model, (0.0, 0.0, 0.0), TIE_BALL_RADIUS, TIE_COCKPIT)
    for site in list(model):
        if site[0] >= TIE_VIEWPORT_DEPTH and abs(site[1]) <= 1 and abs(site[2]) <= 1:
            model[site] = TIE_VIEWPORT
    ball_reach = int(math.floor(TIE_BALL_RADIUS))
    for side in (-1, 1):
        fill_box(model, (-TIE_PYLON_HALF_WIDTH, TIE_PYLON_HALF_WIDTH),
                 (min(side * (ball_reach + 1), side * (TIE_WING_Y - 1)),
                  max(side * (ball_reach + 1), side * (TIE_WING_Y - 1))),
                 (-TIE_PYLON_HALF_WIDTH, TIE_PYLON_HALF_WIDTH), TIE_PYLON)
        reach = int(math.ceil(TIE_WING_RADIUS))
        for x in range(-reach, reach + 1):
            for z in range(-reach, reach + 1):
                if in_tie_wing(x, z):
                    model[(x, side * TIE_WING_Y, z)] = tie_wing_bead(x, z)
    return model


# ---------------------------------------------------------------------------
# the X-wing
# ---------------------------------------------------------------------------
def x_wing_model() -> Model:
    """
    Draw an X-wing as tagged lattice sites.

    A tapering nose, a cockpit with a raised canopy and an astromech behind it,
    a wider rear fuselage, four S-foils in the open (X) position with red
    stripes, an engine at the root of each foil and a laser cannon at each tip.

    Returns
    -------
    dict
        Site (local frame, nose toward +x) -> bead type.
    """
    model: Model = {}
    z_lo, z_hi = XW_FUSELAGE_Z
    # the nose tapers from the cockpit to a 1-wide, 2-tall tip
    nose_length = XW_NOSE_X[1] - XW_NOSE_X[0]
    for x in range(XW_NOSE_X[0], XW_NOSE_X[1] + 1):
        fraction = (x - XW_NOSE_X[0]) / nose_length
        half_width = round_half_up(XW_NOSE_ROOT_HALF_WIDTH * (1 - fraction))
        bottom = z_lo + round_half_up(1 * fraction)
        top = z_hi - 1 - round_half_up(1 * fraction)
        fill_box(model, (x, x), (-half_width, half_width), (bottom, top), XWING_HULL)
    fill_box(model, XW_COCKPIT_X, (-XW_COCKPIT_HALF_WIDTH, XW_COCKPIT_HALF_WIDTH), XW_FUSELAGE_Z, XWING_HULL)
    fill_box(model, XW_CANOPY_X, (-1, 1), (z_hi + 1, z_hi + 1), XWING_CANOPY)
    fill_box(model, (XW_DROID_X[0], XW_CANOPY_X[0] - 1), (-1, 1), (z_hi + 1, z_hi + 1), XWING_HULL)
    fill_box(model, XW_DROID_X, (-1, 1), (z_hi + 1, z_hi + 1), XWING_DROID)
    fill_box(model, XW_DROID_X, (0, 0), (z_hi + 2, z_hi + 2), XWING_DROID)
    fill_box(model, XW_REAR_X, (-XW_REAR_HALF_WIDTH, XW_REAR_HALF_WIDTH), XW_FUSELAGE_Z, XWING_HULL)

    span = XW_WING_Y[1] - XW_WING_Y[0]
    for y_sign, z_sign in itertools.product((-1, 1), repeat=2):
        # the S-foil: a sloped plate, one site thick, tapering toward the tip
        for y_abs in range(XW_WING_Y[0], XW_WING_Y[1] + 1):
            fraction = (y_abs - XW_WING_Y[0]) / span
            z = z_sign * round_half_up((y_abs - XW_WING_Y[0] + 1) * XW_WING_SLOPE)
            lead = round_half_up(XW_WING_ROOT_X[1] + (XW_WING_TIP_X[1] - XW_WING_ROOT_X[1]) * fraction)
            trail = round_half_up(XW_WING_ROOT_X[0] + (XW_WING_TIP_X[0] - XW_WING_ROOT_X[0]) * fraction)
            in_stripe = XW_STRIPE_Y[0] <= y_abs <= XW_STRIPE_Y[1]
            for x in range(trail, lead + 1):
                model[(x, y_sign * y_abs, z)] = XWING_STRIPE if in_stripe else XWING_HULL
        # the engine at the root of the foil, glowing at the back
        y_centre = y_sign * XW_ENGINE_Y
        z_centre = z_sign * XW_ENGINE_Z
        for x in range(XW_ENGINE_X[0], XW_ENGINE_X[1] + 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    glowing = x == XW_ENGINE_X[0] and dy == 0 and dz == 0
                    model[(x, y_centre + dy, z_centre + dz)] = XWING_GLOW if glowing else XWING_ENGINE
        # the cannon along the tip: two sites wide, so a chain can run out and back
        tip_z = z_sign * round_half_up((span + 1) * XW_WING_SLOPE)
        for x in range(XW_CANNON_X[0], XW_CANNON_X[1] + 1):
            model[(x, y_sign * XW_WING_Y[1], tip_z)] = XWING_CANNON
            model[(x, y_sign * (XW_WING_Y[1] + 1), tip_z)] = XWING_CANNON
    return model


# ---------------------------------------------------------------------------
# threading a ship's sites into one chain
# ---------------------------------------------------------------------------
def thread_single_chain(sites: Iterable[Site], rng: random.Random) -> List[Site]:
    """
    Order a set of lattice sites into one chain that visits each site exactly once.

    This is a Hamiltonian path on the graph whose edges join sites that are
    bonded in PIMMS (Chebyshev distance 1). It's found greedily with
    Warnsdorff's rule (always step to the unvisited neighbour with the fewest
    unvisited neighbours of its own), with two repairs:

    * when the walk gets stuck, a Posa rotation reverses the tail of the path
      after a site bonded to the current end, giving a new end that can
      continue (Posa, Discrete Math. 14, 359-364, 1976);
    * a site left behind on its own is spliced in between two consecutive
      chain sites that are both bonded to it.

    Parameters
    ----------
    sites : iterable of tuple of int
        The sites to thread. They must form one bonded piece.
    rng : random.Random
        Breaks ties, so the result is reproducible for a given seed.

    Returns
    -------
    list of tuple of int
        Every site once, each bonded to the next.

    Raises
    ------
    RuntimeError
        If no single chain was found (some sites could not be threaded).
    """
    site_set: Set[Site] = set(sites)
    neighbours: Dict[Site, List[Site]] = {}
    for site in site_set:
        neighbours[site] = [add_sites(site, step) for step in BOND_STEPS if add_sites(site, step) in site_set]
    free_degree: Dict[Site, int] = {site: len(bonded) for site, bonded in neighbours.items()}
    visited: Set[Site] = set()
    path: List[Site] = []
    index_of: Dict[Site, int] = {}

    def visit(site: Site) -> None:
        visited.add(site)
        index_of[site] = len(path)
        path.append(site)
        for neighbour in neighbours[site]:
            free_degree[neighbour] -= 1

    def reverse_tail(start: int) -> None:
        path[start:] = path[start:][::-1]
        for position in range(start, len(path)):
            index_of[path[position]] = position

    def rotate_to_open_end() -> bool:
        # a Posa rotation whose new end still has unvisited neighbours
        end = path[-1]
        last = len(path) - 1
        options = [index_of[n] for n in neighbours[end] if n in visited and index_of[n] < last - 1]
        good = [i for i in options if free_degree[path[i + 1]] > 0]
        if good:
            reverse_tail(max(good) + 1)
            return True
        if free_degree[path[0]] > 0:
            path.reverse()
            for position, site in enumerate(path):
                index_of[site] = position
            return True
        if not options:
            return False
        # no rotation opens the end at once: take a short random one and look again
        options.sort(reverse=True)
        reverse_tail(rng.choice(options[:3]) + 1)
        return True

    def stranded_by(site: Site) -> int:
        # unvisited sites whose only unvisited neighbour is `site`: stepping
        # onto `site` leaves each of them reachable only as the very next step
        return sum(1 for n in neighbours[site] if n not in visited and free_degree[n] == 1)

    start_site = min(site_set, key=lambda site: (free_degree[site], site))
    visit(start_site)
    rotations = 0
    fruitless_rotations = 0
    while len(path) < len(site_set) and rotations < MAX_ROTATIONS:
        end = path[-1]
        options = [n for n in neighbours[end] if n not in visited]
        if options:
            next_site = min(options, key=lambda n: (stranded_by(n), free_degree[n], rng.random()))
            visit(next_site)
            fruitless_rotations = 0
            continue
        rotations += 1
        fruitless_rotations += 1
        if fruitless_rotations > MAX_FRUITLESS_ROTATIONS or not rotate_to_open_end():
            break

    # splice in any sites left behind, one at a time, until none can be placed
    unvisited = [site for site in sorted(site_set) if site not in visited]
    progress = True
    while unvisited and progress:
        progress = False
        for orphan in list(unvisited):
            if not splice_in(orphan, path, index_of, neighbours):
                continue
            visited.add(orphan)
            unvisited.remove(orphan)
            progress = True
    if unvisited:
        raise RuntimeError("could not thread %d of %d sites into one chain (first few: %s)"
                           % (len(unvisited), len(site_set), unvisited[:5]))
    return path


def splice_in(orphan: Site, path: List[Site], index_of: Dict[Site, int],
              neighbours: Dict[Site, List[Site]]) -> bool:
    """
    Splice an unvisited site into the chain, if there's a way to.

    Three ways are tried, simplest first:

    * at either end of the chain, if the orphan is bonded to that end;
    * between two consecutive chain sites that are both bonded to it;
    * a 2-opt splice: with ``a = path[i]`` and ``b = path[j]`` (``i < j``) both
      bonded to the orphan, the chain becomes ``... a, orphan, b, path[j-1], ...,
      path[i+1], path[j+1] ...``, which works when ``path[i+1]`` is bonded to
      ``path[j+1]`` (or ``b`` is the last site).

    Parameters
    ----------
    orphan : tuple of int
        The site to splice in.
    path : list of tuple of int
        The chain so far; changed in place on success.
    index_of : dict
        Site -> its index in ``path``; kept up to date.
    neighbours : dict
        Site -> the sites bonded to it.

    Returns
    -------
    bool
        True if the orphan was spliced in.
    """
    bonded = set(neighbours[orphan])
    if path[-1] in bonded:
        new_path = path + [orphan]
    elif path[0] in bonded:
        new_path = [orphan] + path
    else:
        positions = sorted(index_of[n] for n in bonded if n in index_of)
        new_path = []
        for position in positions:
            if position + 1 < len(path) and path[position + 1] in bonded:
                new_path = path[:position + 1] + [orphan] + path[position + 1:]
                break
        if not new_path:
            for i, j in itertools.combinations(positions, 2):
                if j + 1 == len(path) or are_bonded(path[i + 1], path[j + 1]):
                    new_path = path[:i + 1] + [orphan] + path[i + 1:j + 1][::-1] + path[j + 1:]
                    break
        if not new_path:
            return False
    path[:] = new_path
    for position, site in enumerate(path):
        index_of[site] = position
    return True


def check_chain(path: Sequence[Site], model: Model) -> None:
    """
    Check that a threaded chain covers its ship exactly and is bonded throughout.

    Parameters
    ----------
    path : sequence of tuple of int
        The chain.
    model : dict
        The ship it was threaded from.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If a site is repeated or missing, or two consecutive sites aren't bonded.
    """
    if len(set(path)) != len(path) or set(path) != set(model):
        raise RuntimeError("the chain does not visit every site of the ship exactly once")
    for site_a, site_b in zip(path, path[1:]):
        if not are_bonded(site_a, site_b):
            raise RuntimeError("chain sites %s and %s are not bonded" % (site_a, site_b))


# ---------------------------------------------------------------------------
# lattice rotations and placing the fleet
# ---------------------------------------------------------------------------
def proper_rotations() -> List[Rotation]:
    """
    The 24 rotations that map the cubic lattice onto itself.

    Returns
    -------
    list of 3x3 tuple of int
        Signed permutation matrices with determinant +1.
    """
    rotations: List[Rotation] = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            rows = []
            for row in range(3):
                entries = [0, 0, 0]
                entries[permutation[row]] = signs[row]
                rows.append((entries[0], entries[1], entries[2]))
            matrix: Rotation = (rows[0], rows[1], rows[2])
            if determinant(matrix) == 1:
                rotations.append(matrix)
    return rotations


def determinant(matrix: Rotation) -> int:
    """
    Determinant of a 3x3 integer matrix.

    Parameters
    ----------
    matrix : 3x3 tuple of int
        The matrix.

    Returns
    -------
    int
        Its determinant.
    """
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def rotate(site: Site, matrix: Rotation) -> Site:
    """
    Apply a lattice rotation to a site.

    Parameters
    ----------
    site : tuple of int
        The site.
    matrix : 3x3 tuple of int
        The rotation.

    Returns
    -------
    tuple of int
        ``matrix @ site``.
    """
    return (sum(matrix[0][k] * site[k] for k in range(3)),
            sum(matrix[1][k] * site[k] for k in range(3)),
            sum(matrix[2][k] * site[k] for k in range(3)))


IDENTITY: Rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
HALF_TURN_ABOUT_Z: Rotation = ((-1, 0, 0), (0, -1, 0), (0, 0, 1))


def place_fighter(path: Sequence[Site], matrix: Rotation, region: Tuple[int, int, int, int, int, int],
                  occupied: Set[Site], rng: random.Random) -> List[Site]:
    """
    Find a spot for one fighter that is inside the box and touches nothing.

    The fighter is rotated, then shifted to random centres drawn from
    ``region`` until none of its sites is on, or bonded to, an occupied site,
    so every fighter starts free of every other ship.

    Parameters
    ----------
    path : sequence of tuple of int
        The fighter's chain, in its local frame.
    matrix : 3x3 tuple of int
        Rotation applied before shifting.
    region : tuple of int
        (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi), inclusive, for the shift.
    occupied : set of tuple of int
        Sites already taken (box coordinates). Not changed here.
    rng : random.Random
        Source of randomness.

    Returns
    -------
    list of tuple of int
        The fighter's chain in box coordinates.

    Raises
    ------
    RuntimeError
        If no free spot is found in ``MAX_PLACEMENT_TRIES`` draws.
    """
    rotated = [rotate(site, matrix) for site in path]
    x_lo, x_hi, y_lo, y_hi, z_lo, z_hi = region
    for _ in range(MAX_PLACEMENT_TRIES):
        shift = (rng.randint(x_lo, x_hi), rng.randint(y_lo, y_hi), rng.randint(z_lo, z_hi))
        placed = [add_sites(site, shift) for site in rotated]
        if not all(0 <= site[k] < BOX_DIMENSIONS[k] for site in placed for k in range(3)):
            continue
        touching = any(add_sites(site, step) in occupied for site in placed for step in BOND_STEPS + ((0, 0, 0),))
        if not touching:
            return placed
    raise RuntimeError("no room found for a fighter in region %s" % (region,))


# ---------------------------------------------------------------------------
# writing the files
# ---------------------------------------------------------------------------
def write_parameter_file(filename: str) -> None:
    """
    Write a parameter file in which every interaction is zero.

    The fleet is pure excluded volume: ships bounce off each other and the
    frozen Star Destroyer, but nothing attracts or repels.

    Parameters
    ----------
    filename : str
        Where to write it.

    Returns
    -------
    None
    """
    lines = ["##",
             "## Star Destroyer fleet demo parameter file (written by build_fleet.py).",
             "##",
             "## Every interaction is ZERO: the system is pure excluded volume, so the",
             "## fighters fly about freely and bounce off each other and off the frozen",
             "## Star Destroyer. The many bead types are only there to colour the parts.",
             "##",
             "",
             "## ANGLE PENALTIES (every bead type needs one, even if zero)"]
    lines += ["ANGLE_PENALTY\t%s\t0\t0\t0" % bead for bead in ALL_BEAD_TYPES]
    lines += ["", "## PAIRWISE INTERACTIONS (all zero)"]
    for index, bead_a in enumerate(ALL_BEAD_TYPES):
        for bead_b in ALL_BEAD_TYPES[index:]:
            lines.append("%s\t%s\t0" % (bead_a, bead_b))
    lines += ["", "## BEAD-SOLVENT INTERACTIONS (all zero)"]
    lines += ["%s\t0\t0" % bead for bead in ALL_BEAD_TYPES]
    with open(filename, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def write_restart_and_freeze(chains: Dict[int, Tuple[List[Site], str, int]]) -> None:
    """
    Write the restart file holding the fleet and the freeze file for the Star Destroyer.

    Parameters
    ----------
    chains : dict
        chainID -> (box positions in chain order, sequence, chainType).

    Returns
    -------
    None
    """
    restart = {
        "DIMENSIONS": list(BOX_DIMENSIONS),
        "ENERGY": 0,
        "HARDWALL": True,
        "CHAINS": {chain_id: [[list(site) for site in positions], sequence, chain_type]
                   for chain_id, (positions, sequence, chain_type) in chains.items()},
    }
    with open(RESTART_FILENAME, "wb") as handle:
        pickle.dump(restart, handle)
    with open(FREEZE_FILENAME, "w") as handle:
        handle.write("# freeze the Star Destroyer (chain 1, the whole ship); the fighters fly\n")
        handle.write("C 1\n")


def build_ship(name: str, make_model: Callable[[], Model], rng: random.Random) -> Tuple[List[Site], str]:
    """
    Draw one ship, thread it into a single chain and check the result.

    Parameters
    ----------
    name : str
        For the progress message.
    make_model : callable
        Returns the ship's tagged sites.
    rng : random.Random
        Breaks ties in the threading.

    Returns
    -------
    tuple
        (chain sites in the local frame, sequence of bead types along the chain).
    """
    model = make_model()
    path = thread_single_chain(model, rng)
    check_chain(path, model)
    sequence = "".join(model[site] for site in path)
    print("%-15s %6d beads, one chain" % (name, len(path)))
    return (path, sequence)


def main() -> None:
    """
    Build the fleet and write ``fleet.restart``, ``freeze.txt`` and ``params.prm``.

    Returns
    -------
    None
    """
    rng = random.Random(BUILD_SEED)
    sd_path, sd_sequence = build_ship("Star Destroyer", lambda: star_destroyer_model(rng), rng)
    tie_path, tie_sequence = build_ship("TIE fighter", tie_fighter_model, rng)
    xwing_path, xwing_sequence = build_ship("X-wing", x_wing_model, rng)

    chains: Dict[int, Tuple[List[Site], str, int]] = {}
    sd_box = [add_sites(site, SD_BOX_ORIGIN) for site in sd_path]
    if not all(0 <= site[k] < BOX_DIMENSIONS[k] for site in sd_box for k in range(3)):
        raise RuntimeError("the Star Destroyer does not fit in the box")
    chains[1] = (sd_box, sd_sequence, 0)
    occupied: Set[Site] = set(sd_box)

    next_id = 2
    # TIEs face back along -x (toward the oncoming X-wings); X-wings face +x
    for chain_type, path, sequence, matrix, region, count in (
            (1, tie_path, tie_sequence, HALF_TURN_ABOUT_Z, TIE_REGION, N_TIE_FIGHTERS),
            (2, xwing_path, xwing_sequence, IDENTITY, XWING_REGION, N_X_WINGS)):
        for _ in range(count):
            placed = place_fighter(path, matrix, region, occupied, rng)
            occupied.update(placed)
            chains[next_id] = (placed, sequence, chain_type)
            next_id += 1

    write_restart_and_freeze(chains)
    write_parameter_file(PARAMETER_FILENAME)
    n_beads = sum(len(positions) for positions, _, _ in chains.values())
    print("fleet: %d chains (1 Star Destroyer, %d TIE fighters, %d X-wings), %d beads in a %dx%dx%d box"
          % (len(chains), N_TIE_FIGHTERS, N_X_WINGS, n_beads, *BOX_DIMENSIONS))
    print("wrote %s, %s and %s" % (RESTART_FILENAME, FREEZE_FILENAME, PARAMETER_FILENAME))


if __name__ == "__main__":
    main()
