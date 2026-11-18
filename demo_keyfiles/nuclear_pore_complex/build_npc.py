#!/usr/bin/env python
"""
Build the nuclear pore complex demo: a schematic human nuclear pore complex (NPC)
in a nuclear envelope, with mobile FG-nucleoporins, transport receptors and inert
cargo.

Everything structural is frozen; only disordered and soluble things move:

* the nuclear envelope (outer membrane, pore membrane and inner membrane) and
  the NPC scaffold (inner ring, cytoplasmic and nuclear rings of Y-complexes,
  nuclear basket, luminal Gp210 ring) are frozen chains, one per component;
* FG-nucleoporins are chains anchored to the scaffold by their first bead, which
  is frozen (``B`` lines in ``freeze.txt``); the rest of each FG domain moves;
* transport receptors (sticky for FG motifs) and inert cargo of the same size
  start in the cytoplasm and are free.

The geometry is a schematic built from published dimensions (sources at each
constant), at 2 nm per lattice site. It is not a reconstruction of any one
structure: the scaffold is drawn from simple solids with the right sizes and
8-fold symmetry, and the FG domains and interactions are illustrative.

Run from this directory:   ``python build_npc.py``
Outputs:                   ``npc.restart``, ``freeze.txt`` and ``params.prm``

The same ``BUILD_SEED`` always gives the same files.
"""

from __future__ import annotations

import itertools
import math
import pickle
import random
from typing import Callable, Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np
from scipy import ndimage

Site = Tuple[int, int, int]
Model = Dict[Site, str]
ChainRecord = Tuple[List[Site], str, int]

# ---------------------------------------------------------------------------
# output files, the box and the lattice scale
# ---------------------------------------------------------------------------
RESTART_FILENAME = "npc.restart"
FREEZE_FILENAME = "freeze.txt"
PARAMETER_FILENAME = "params.prm"
BUILD_SEED = 2022

NM_PER_SITE = 2.0               # one lattice site is 2 nm (LATTICE_TO_ANGSTROMS : 20 in KEYFILE.kf)
BOX_HALF_WIDTH = 46             # the box runs from -46 to +46 across (93 sites = 186 nm)
BOX_Z: Tuple[int, int] = (-58, 52)   # local z range: nucleus below, cytoplasm above
BOX_DIMENSIONS: Site = (2 * BOX_HALF_WIDTH + 1, 2 * BOX_HALF_WIDTH + 1, BOX_Z[1] - BOX_Z[0] + 1)

# ---------------------------------------------------------------------------
# literature dimensions (nm), converted to sites below. Local frame: the pore
# axis is the z axis through (0, 0); z = 0 is the middle of the nuclear
# envelope; the cytoplasm is at +z and the nucleus at -z.
# ---------------------------------------------------------------------------
# Mosalaganti et al., Science 376, eabm9506 (2022); preprint bioRxiv
# 10.1101/2021.10.26.465776: "a 120 nm-wide protein scaffold of three stacked
# rings", each of "eight spokes"; the dilated in-cell inner ring is 54 nm across
# (42 nm constricted), and the membrane pore 92 nm (82 nm constricted, Fig. 1).
NPC_DIAMETER_NM = 120.0
CHANNEL_DIAMETER_NM = 54.0
MEMBRANE_PORE_DIAMETER_NM = 92.0
N_SPOKES = 8
# Hetzer, Cold Spring Harb. Perspect. Biol. 2, a000539 (2010): the inner and
# outer nuclear membranes are kept "at an even distance of ~50 nm".
MEMBRANE_SPACING_NM = 50.0
# Mitra et al., PNAS 101, 4083-4088 (2004): the rough ER bilayer (continuous
# with the nuclear envelope) is 37.5 A phosphate to phosphate - about 4 nm.
BILAYER_THICKNESS_NM = 4.0

CHANNEL_RADIUS = CHANNEL_DIAMETER_NM / 2 / NM_PER_SITE              # 13.5 sites
PORE_WAIST_RADIUS = MEMBRANE_PORE_DIAMETER_NM / 2 / NM_PER_SITE     # 23 sites, mid-membrane
MEMBRANE_HALF_SPACING = MEMBRANE_SPACING_NM / 2 / NM_PER_SITE       # 12.5 sites, mid-plane to mid-membrane
MEMBRANE_HALF_THICKNESS = BILAYER_THICKNESS_NM / 2 / NM_PER_SITE    # 1 site
# the pore membrane is the inner half of a torus joining the two membranes
TORUS_MAJOR_RADIUS = PORE_WAIST_RADIUS + MEMBRANE_HALF_SPACING      # 35.5 sites
NPC_RADIUS = NPC_DIAMETER_NM / 2 / NM_PER_SITE                      # 30 sites

# ---------------------------------------------------------------------------
# scaffold geometry (sites). Illustrative values chosen to fit the dimensions
# above; the outer rings' Y-complex layers end up 70 nm apart, top to bottom
# (82 nm including the domains on top of them).
# ---------------------------------------------------------------------------
IR_HALF_HEIGHT = 9              # the inner ring spans z = -9..9 (38 nm)
IR_SHELL_THICKNESS = 2          # the inner ring is a hollow shell this deep (its inside is sealed off)
IR_GROOVE_HALF_WIDTH = 0.8      # shallow grooves on the channel face mark the 8 spokes
OUTER_RING_Z: Tuple[int, int] = (15, 17)     # cytoplasmic ring; the nuclear ring is its mirror image
Y_RING_RADII: Tuple[float, float] = (17.5, 24.5)   # the two concentric, head-to-tail rings of Y-complexes
Y_STEM_HALF_WIDTH = 1.0
ROD_RADIUS = 1.5                # arms, bridges, legs and connectors: at least 2 sites wide everywhere,
                                # so a chain can run out along each and back (a 1-wide tip forces a chain end)
Y_ARM_LENGTH = 5.0
Y_ARM_ANGLE_DEG = 35.0          # arms splay this far from the ring tangent, one in, one out
Y_HEAD_ANGLE_DEG = 40.0         # each Y's head sits this far round its spoke
BRIDGE_ANGLE_DEG = 20.0         # Nup93 bridges join the inner and outer Y rings here
LEG_ANGLES_DEG: Tuple[float, ...] = (10.0, 30.0)    # membrane-binding legs under the outer ring
CONNECTOR_ANGLE_DEG = 25.0      # Nup155 connectors from the outer rings down to the inner ring
CONNECTOR_RADIUS = 16.0
DOMAIN_RADIUS = 1.8             # small globular domains (Nup358, Nup214, ELYS)
NUP358_ANGLE_DEG = 15.0
NUP214_ANGLE_DEG = 35.0
ELYS_ANGLE_DEG = 30.0
BASKET_TOP: Tuple[float, float] = (17.5, -17.0)     # (radius, z) where each basket filament leaves the nuclear ring
BASKET_BOTTOM: Tuple[float, float] = (9.0, -46.0)   # (radius, z) of the distal ring: a basket ~60 nm long (illustrative)
BASKET_ANGLE_DEG = 22.5
BASKET_FILAMENT_RADIUS = 1.5
LUMINAL_RING_RADIUS = 28.0      # the Gp210 ring, in the lumen just outside the pore membrane
LUMINAL_RING_TUBE = 1.5
LUMINAL_ARCH_ANGLE_DEG = 22.5   # where each spoke's Gp210 reaches through to the pore membrane

# ---------------------------------------------------------------------------
# bead types: one letter each, written as a PDB residue name (see readme.md)
# ---------------------------------------------------------------------------
MEMBRANE = "M"           # MET  nuclear envelope (both membranes and the pore membrane)
IR_EVEN = "I"            # ILE  inner ring, even spokes
IR_ODD = "V"             # VAL  inner ring, odd spokes
CR = "C"                 # CYS  cytoplasmic ring (Y-complexes, bridges, legs, connectors)
NUP358_DOMAIN = "R"      # ARG  Nup358 (RanBP2) domains on the cytoplasmic ring
NUP214_DOMAIN = "W"      # TRP  Nup214 complex domains on the cytoplasmic ring
NR = "N"                 # ASN  nuclear ring
ELYS_DOMAIN = "E"        # GLU  ELYS on the nuclear ring
BASKET = "T"             # THR  nuclear basket (Tpr) and its distal ring
LUMINAL = "L"            # LEU  luminal ring (Gp210)
FG_ANCHOR = "A"          # ALA  first bead of each FG-Nup: frozen onto the scaffold
FXFG = "F"               # PHE  an FxFG/FG motif
GLFG = "Q"               # GLN  a GLFG motif (Nup98)
FG_SPACER = "S"          # SER  spacer between FG motifs
RECEPTOR = "K"           # LYS  transport receptor (binds FG motifs)
CARGO = "D"              # ASP  inert cargo (binds nothing)
ALL_BEAD_TYPES: Tuple[str, ...] = (MEMBRANE, IR_EVEN, IR_ODD, CR, NUP358_DOMAIN, NUP214_DOMAIN, NR, ELYS_DOMAIN,
                                   BASKET, LUMINAL, FG_ANCHOR, FXFG, GLFG, FG_SPACER, RECEPTOR, CARGO)

# ---------------------------------------------------------------------------
# interactions (energies in units where kT = TEMPERATURE in KEYFILE.kf, so
# -20 at TEMPERATURE 100 is -0.2 kT per contact). Illustrative: FG motifs are
# cohesive (GLFG more than FxFG), receptors bind FG motifs, cargo binds nothing,
# and the scaffold and membrane are inert walls. Every pair not listed is zero.
# Receptor binding was tuned by hand: at -0.6 kT per contact receptors pile into
# the FG meshwork and stay there; at -0.4 kT they pass through.
# ---------------------------------------------------------------------------
INTERACTIONS: Dict[Tuple[str, str], int] = {
    (FXFG, FXFG): -20,
    (FXFG, GLFG): -25,
    (GLFG, GLFG): -30,
    (RECEPTOR, FXFG): -40,
    (RECEPTOR, GLFG): -40,
}

# ---------------------------------------------------------------------------
# FG-nucleoporins. One bead stands for ~5 residues, and every third bead is an
# FG motif (one per ~15 residues). Lengths are illustrative, of the order of
# the disordered FG domains (a few hundred residues). Per-pore copy numbers of
# Nups come in multiples of 8 (Ori et al., Mol. Syst. Biol. 9, 648, 2013); the
# counts per spoke here are illustrative.
# (name, motif bead, beads per chain, where it's anchored, copies per spoke)
# ---------------------------------------------------------------------------
FG_MOTIF_SPACING = 3
FG_NUPS: Tuple[Tuple[str, str, int, str, int], ...] = (
    ("Nup62", FXFG, 40, "inner ring", 4),
    ("Nup98", GLFG, 80, "inner ring", 4),
    ("Nup214", FXFG, 60, "cytoplasmic ring", 1),
    ("Nup358", FXFG, 64, "cytoplasmic ring top", 1),
    ("Nup153", FXFG, 100, "nuclear ring", 1),
)

# ---------------------------------------------------------------------------
# transport receptors and inert cargo: same length, so the only difference is
# whether they bind FG motifs. Both start in the cytoplasm.
# ---------------------------------------------------------------------------
N_RECEPTORS = 24
N_CARGO = 24
SOLUBLE_LENGTH = 6
SOLUBLE_START_Z: Tuple[int, int] = (28, 48)
MAX_GROWTH_TRIES = 2000

# the threader: how many Posa rotations it may try in all, and in a row
MAX_ROTATIONS = 200000
MAX_FRUITLESS_ROTATIONS = 2000

# the 26 bond vectors, face steps first so ties favour tidy, axis-aligned bonds
BOND_STEPS: Tuple[Site, ...] = tuple(sorted(
    ((dx, dy, dz) for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3) if (dx, dy, dz) != (0, 0, 0)),
    key=lambda step: (abs(step[0]) + abs(step[1]) + abs(step[2]), step)))


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
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
    return max(abs(a - b) for a, b in zip(site_a, site_b)) == 1


def in_box(site: Site) -> bool:
    """
    Whether a local-frame site lies inside the box.

    Parameters
    ----------
    site : tuple of int
        A site in the local frame (pore axis through x = y = 0).

    Returns
    -------
    bool
        True inside the box.
    """
    x, y, z = site
    return abs(x) <= BOX_HALF_WIDTH and abs(y) <= BOX_HALF_WIDTH and BOX_Z[0] <= z <= BOX_Z[1]


def cylindrical(site: Site) -> Tuple[float, float, int]:
    """
    Cylindrical coordinates of a site about the pore axis.

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    tuple
        (radius in sites, angle in degrees in [0, 360), z).
    """
    x, y, z = site
    return (math.hypot(x, y), math.degrees(math.atan2(y, x)) % 360.0, z)


def polar_site(radius: float, angle_deg: float, z: float) -> Tuple[float, float, float]:
    """
    The Cartesian point at a given radius, angle and height.

    Parameters
    ----------
    radius : float
        Distance from the pore axis, in sites.
    angle_deg : float
        Angle round the axis, in degrees.
    z : float
        Height.

    Returns
    -------
    tuple of float
        (x, y, z).
    """
    angle = math.radians(angle_deg)
    return (radius * math.cos(angle), radius * math.sin(angle), z)


def spoke_angles() -> List[float]:
    """
    The angle of each spoke's centre, in degrees.

    Returns
    -------
    list of float
        ``N_SPOKES`` angles, 0, 45, 90, ...
    """
    return [index * 360.0 / N_SPOKES for index in range(N_SPOKES)]


def add_segment(model: Model, start: Tuple[float, float, float], end: Tuple[float, float, float],
                radius: float, bead: str) -> None:
    """
    Add every site within ``radius`` of a line segment (a rod with round ends).

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    start, end : tuple of float
        The segment's ends.
    radius : float
        Rod radius, in sites.
    bead : str
        Bead type.

    Returns
    -------
    None
    """
    sx, sy, sz = start
    ex, ey, ez = end
    dx, dy, dz = ex - sx, ey - sy, ez - sz
    length_squared = dx * dx + dy * dy + dz * dz
    reach = int(math.ceil(radius)) + 1
    for x in range(int(math.floor(min(sx, ex))) - reach, int(math.ceil(max(sx, ex))) + reach + 1):
        for y in range(int(math.floor(min(sy, ey))) - reach, int(math.ceil(max(sy, ey))) + reach + 1):
            for z in range(int(math.floor(min(sz, ez))) - reach, int(math.ceil(max(sz, ez))) + reach + 1):
                t = 0.0 if length_squared == 0 else ((x - sx) * dx + (y - sy) * dy + (z - sz) * dz) / length_squared
                t = min(1.0, max(0.0, t))
                px, py, pz = sx + t * dx, sy + t * dy, sz + t * dz
                if (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2 <= radius ** 2:
                    model[(x, y, z)] = bead


def add_sphere(model: Model, centre: Tuple[float, float, float], radius: float, bead: str) -> None:
    """
    Add every site within ``radius`` of a point.

    Parameters
    ----------
    model : dict
        Site -> bead type; changed in place.
    centre : tuple of float
        The centre.
    radius : float
        Radius, in sites.
    bead : str
        Bead type.

    Returns
    -------
    None
    """
    add_segment(model, centre, centre, radius, bead)


def mirror_z(model: Model, bead_map: Dict[str, str]) -> Model:
    """
    Reflect a model through z = 0, renaming its bead types.

    Parameters
    ----------
    model : dict
        Site -> bead type.
    bead_map : dict
        Old bead type -> new bead type.

    Returns
    -------
    dict
        The reflected model.
    """
    return {(x, y, -z): bead_map[bead] for (x, y, z), bead in model.items()}


# ---------------------------------------------------------------------------
# the nuclear envelope
# ---------------------------------------------------------------------------
def is_membrane(site: Site) -> bool:
    """
    Whether a site belongs to the nuclear envelope.

    Away from the pore, the outer and inner membranes are flat sheets
    ``2 * MEMBRANE_HALF_THICKNESS`` thick, centred at z = +/- MEMBRANE_HALF_SPACING.
    Nearer the axis than ``TORUS_MAJOR_RADIUS`` they curve into the pore
    membrane: the inner half of a torus that joins them, whose waist is
    ``PORE_WAIST_RADIUS`` from the axis.

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    bool
        True if the site is membrane.
    """
    radius, _, z = cylindrical(site)
    distance: float
    if radius >= TORUS_MAJOR_RADIUS:
        distance = abs(z)
    else:
        distance = math.hypot(radius - TORUS_MAJOR_RADIUS, z)
    return abs(distance - MEMBRANE_HALF_SPACING) <= MEMBRANE_HALF_THICKNESS


def in_pore(site: Site) -> bool:
    """
    Whether a site is inside the membrane pore (on the pore side of the pore membrane).

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    bool
        True if the site is nearer the axis than the pore membrane and between the
        two membranes' heights.
    """
    radius, _, z = cylindrical(site)
    if radius >= TORUS_MAJOR_RADIUS or abs(z) > MEMBRANE_HALF_SPACING + MEMBRANE_HALF_THICKNESS:
        return False
    return math.hypot(radius - TORUS_MAJOR_RADIUS, z) > MEMBRANE_HALF_SPACING + MEMBRANE_HALF_THICKNESS


def nuclear_envelope_model() -> Model:
    """
    Draw the nuclear envelope across the whole box, pierced by the pore.

    Returns
    -------
    dict
        Site -> ``MEMBRANE``.
    """
    model: Model = {}
    reach = int(MEMBRANE_HALF_SPACING + MEMBRANE_HALF_THICKNESS) + 1
    for x in range(-BOX_HALF_WIDTH, BOX_HALF_WIDTH + 1):
        for y in range(-BOX_HALF_WIDTH, BOX_HALF_WIDTH + 1):
            for z in range(-reach, reach + 1):
                if is_membrane((x, y, z)):
                    model[(x, y, z)] = MEMBRANE
    return model


# ---------------------------------------------------------------------------
# the scaffold
# ---------------------------------------------------------------------------
def inner_ring_model() -> Model:
    """
    Draw the inner ring: it fills the pore from the channel wall to the pore membrane.

    The ring is solid enough to seal the pore - the only way through is the
    central channel - and hollow inside (a shell ``IR_SHELL_THICKNESS`` deep).
    The eight spokes alternate bead type, and shallow grooves on the channel
    face mark where they meet.

    Returns
    -------
    dict
        Site -> ``IR_EVEN`` or ``IR_ODD``.
    """
    model: Model = {}
    reach = int(TORUS_MAJOR_RADIUS) + 1
    spoke_width = 360.0 / N_SPOKES
    for x in range(-reach, reach + 1):
        for y in range(-reach, reach + 1):
            for z in range(-IR_HALF_HEIGHT, IR_HALF_HEIGHT + 1):
                site = (x, y, z)
                radius, angle, _ = cylindrical(site)
                if radius < CHANNEL_RADIUS or not in_pore(site):
                    continue
                spoke = int(((angle + spoke_width / 2) % 360.0) // spoke_width)
                # spokes meet at angles (k + 1/2) * spoke_width; how far round
                # the ring (in sites) this site is from the nearest such boundary
                offset = (angle - spoke_width / 2) % spoke_width
                along = math.radians(min(offset, spoke_width - offset)) * radius
                if radius < CHANNEL_RADIUS + 1.5 and along <= IR_GROOVE_HALF_WIDTH:
                    continue
                model[site] = IR_EVEN if spoke % 2 == 0 else IR_ODD
    hollow_out(model, IR_SHELL_THICKNESS)
    return model


def outer_ring_model(z_range: Tuple[int, int], ring_bead: str, cytoplasmic: bool) -> Model:
    """
    Draw a cytoplasmic-ring-like outer ring at heights ``z_range``.

    Two concentric rings of Y-complexes, arranged head to tail (each Y's stem
    runs along its spoke and its two arms splay out at the head), joined by
    Nup93 bridges, with legs down to the outer membrane and connectors down to
    the inner ring. The cytoplasmic ring also carries Nup358 domains on top and
    Nup214 domains facing the channel; the nuclear ring (drawn here and then
    mirrored by the caller) carries ELYS instead.

    Parameters
    ----------
    z_range : tuple of int
        Bottom and top z of the Y-complex layer (for the cytoplasmic side).
    ring_bead : str
        Bead type of the ring.
    cytoplasmic : bool
        True for the cytoplasmic ring, False for the nuclear ring.

    Returns
    -------
    dict
        Site -> bead type (cytoplasmic orientation: +z).
    """
    model: Model = {}
    z_low, z_high = z_range
    z_mid = (z_low + z_high) / 2
    reach = int(max(Y_RING_RADII) + Y_ARM_LENGTH) + 2
    # the stems: a continuous band round each Y ring
    for x in range(-reach, reach + 1):
        for y in range(-reach, reach + 1):
            radius = math.hypot(x, y)
            for ring_radius in Y_RING_RADII:
                if abs(radius - ring_radius) <= Y_STEM_HALF_WIDTH:
                    for z in range(z_low, z_high):
                        model[(x, y, z)] = ring_bead
    for spoke in spoke_angles():
        for ring_radius in Y_RING_RADII:
            # the two arms at each Y's head
            head_angle = spoke + Y_HEAD_ANGLE_DEG
            head = polar_site(ring_radius, head_angle, z_mid)
            tangent = math.radians(head_angle + 90.0)
            for side in (-1, 1):
                direction = tangent + side * math.radians(90.0 - Y_ARM_ANGLE_DEG)
                end = (head[0] + Y_ARM_LENGTH * math.cos(direction),
                       head[1] + Y_ARM_LENGTH * math.sin(direction), z_mid)
                add_segment(model, head, end, ROD_RADIUS, ring_bead)
        # Nup93 bridge between the two Y rings
        add_segment(model, polar_site(Y_RING_RADII[0], spoke + BRIDGE_ANGLE_DEG, z_mid),
                    polar_site(Y_RING_RADII[1], spoke + BRIDGE_ANGLE_DEG, z_mid), ROD_RADIUS, ring_bead)
        # legs from the outer Y ring down onto the membrane
        for leg_angle in LEG_ANGLES_DEG:
            top = polar_site(Y_RING_RADII[1], spoke + leg_angle, z_low)
            bottom = polar_site(Y_RING_RADII[1], spoke + leg_angle, 0.0)
            add_segment(model, top, bottom, ROD_RADIUS, ring_bead)
        # Nup155-like connectors from the inner Y ring down to the inner ring
        add_segment(model, polar_site(CONNECTOR_RADIUS, spoke + CONNECTOR_ANGLE_DEG, z_low),
                    polar_site(CONNECTOR_RADIUS, spoke + CONNECTOR_ANGLE_DEG, IR_HALF_HEIGHT + 1),
                    ROD_RADIUS, ring_bead)
        if cytoplasmic:
            add_sphere(model, polar_site(Y_RING_RADII[1], spoke + NUP358_ANGLE_DEG, z_high + DOMAIN_RADIUS - 0.5),
                       DOMAIN_RADIUS, NUP358_DOMAIN)
            add_sphere(model, polar_site(Y_RING_RADII[0] - 2.5, spoke + NUP214_ANGLE_DEG, z_high),
                       DOMAIN_RADIUS, NUP214_DOMAIN)
        else:
            add_sphere(model, polar_site(Y_RING_RADII[0], spoke + ELYS_ANGLE_DEG, z_high + DOMAIN_RADIUS - 0.5),
                       DOMAIN_RADIUS, ELYS_DOMAIN)
    # keep only what's on the cytoplasmic side of the envelope and out of the
    # inner ring's space: the legs stop at the membrane, the connectors at the ring
    return {site: bead for site, bead in model.items() if on_cytoplasmic_side(site) and not in_pore_body(site)}


def on_cytoplasmic_side(site: Site) -> bool:
    """
    Whether a site is on the cytoplasmic side of the nuclear envelope.

    That's above the outer membrane, or inside the pore above the envelope's
    mid-plane - not in the membrane, the lumen or anywhere on the nuclear side.

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    bool
        True on the cytoplasmic side.
    """
    radius, _, z = cylindrical(site)
    if z <= 0:
        return False
    outer_surface = MEMBRANE_HALF_SPACING + MEMBRANE_HALF_THICKNESS
    if radius >= TORUS_MAJOR_RADIUS:
        return z > outer_surface
    return math.hypot(radius - TORUS_MAJOR_RADIUS, z) > outer_surface


def in_pore_body(site: Site) -> bool:
    """
    Whether a site lies in the inner ring's volume (whether or not it's hollowed).

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    bool
        True if the inner ring would occupy it.
    """
    radius, _, z = cylindrical(site)
    return abs(z) <= IR_HALF_HEIGHT and radius >= CHANNEL_RADIUS - 1.5 and in_pore(site)


def basket_model() -> Model:
    """
    Draw the nuclear basket: eight Tpr filaments from the nuclear ring to a distal ring.

    Returns
    -------
    dict
        Site -> ``BASKET``.
    """
    model: Model = {}
    for spoke in spoke_angles():
        top = polar_site(BASKET_TOP[0], spoke + BASKET_ANGLE_DEG, BASKET_TOP[1])
        bottom = polar_site(BASKET_BOTTOM[0], spoke + BASKET_ANGLE_DEG, BASKET_BOTTOM[1])
        add_segment(model, top, bottom, BASKET_FILAMENT_RADIUS, BASKET)
    steps = 96
    for index in range(steps):
        start = polar_site(BASKET_BOTTOM[0], index * 360.0 / steps, BASKET_BOTTOM[1])
        end = polar_site(BASKET_BOTTOM[0], (index + 1) * 360.0 / steps, BASKET_BOTTOM[1])
        add_segment(model, start, end, BASKET_FILAMENT_RADIUS, BASKET)
    return model


def luminal_ring_model() -> Model:
    """
    Draw the Gp210 luminal ring, with one connector per spoke through to the pore membrane.

    Returns
    -------
    dict
        Site -> ``LUMINAL``.
    """
    model: Model = {}
    reach = int(LUMINAL_RING_RADIUS + LUMINAL_RING_TUBE) + 1
    for x in range(-reach, reach + 1):
        for y in range(-reach, reach + 1):
            radius = math.hypot(x, y)
            for z in range(-2, 3):
                if math.hypot(radius - LUMINAL_RING_RADIUS, z) <= LUMINAL_RING_TUBE:
                    model[(x, y, z)] = LUMINAL
    pore_membrane_radius = TORUS_MAJOR_RADIUS - MEMBRANE_HALF_SPACING
    for spoke in spoke_angles():
        add_segment(model, polar_site(pore_membrane_radius, spoke + LUMINAL_ARCH_ANGLE_DEG, 0.0),
                    polar_site(LUMINAL_RING_RADIUS, spoke + LUMINAL_ARCH_ANGLE_DEG, 0.0), ROD_RADIUS, LUMINAL)
    return {site: bead for site, bead in model.items() if not is_membrane(site)}


def hollow_out(model: Model, thickness: int) -> None:
    """
    Remove every site deeper than ``thickness`` below the surface.

    A site's depth is its Chebyshev distance to the nearest empty site, so the
    sites kept form a skin ``thickness`` sites deep, and every site removed is
    at least ``thickness + 1`` from open space - out of reach of any bead
    outside, since a bead can only move to a site bonded to a neighbour.

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
    depth = {site: 1 for site in sites if any(add_sites(site, step) not in sites for step in BOND_STEPS)}
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
# threading a component's sites into one chain (the same algorithm as the
# star_destroyer_fleet demo's build_fleet.py)
# ---------------------------------------------------------------------------
def thread_single_chain(sites: Iterable[Site], rng: random.Random) -> List[Site]:
    """
    Order a set of lattice sites into one chain that visits each site exactly once.

    This is a Hamiltonian path on the graph whose edges join sites that are
    bonded in PIMMS (Chebyshev distance 1), found greedily with Warnsdorff's
    rule (step to the unvisited neighbour with the fewest unvisited neighbours,
    avoiding steps that strand a site). When the walk gets stuck, a Posa
    rotation reverses the tail of the path after a site bonded to the end
    (Posa, Discrete Math. 14, 359-364, 1976); sites left over are spliced in.

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
        If no single chain was found.
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
        options.sort(reverse=True)
        reverse_tail(rng.choice(options[:3]) + 1)
        return True

    def stranded_by(site: Site) -> int:
        return sum(1 for n in neighbours[site] if n not in visited and free_degree[n] == 1)

    start_site = min(site_set, key=lambda site: (free_degree[site], site))
    visit(start_site)
    rotations = 0
    fruitless_rotations = 0
    while len(path) < len(site_set) and rotations < MAX_ROTATIONS:
        options = [n for n in neighbours[path[-1]] if n not in visited]
        if options:
            visit(min(options, key=lambda n: (stranded_by(n), free_degree[n], rng.random())))
            fruitless_rotations = 0
            continue
        rotations += 1
        fruitless_rotations += 1
        if fruitless_rotations > MAX_FRUITLESS_ROTATIONS or not rotate_to_open_end():
            break

    unvisited = [site for site in sorted(site_set) if site not in visited]
    progress = True
    while unvisited and progress:
        progress = False
        for orphan in list(unvisited):
            if splice_in(orphan, path, index_of, neighbours):
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

    Tried in order: at either end of the chain; between two consecutive chain
    sites both bonded to it; or a 2-opt splice ``... a, orphan, b, path[j-1],
    ..., path[i+1], path[j+1] ...`` for ``a = path[i]``, ``b = path[j]`` both
    bonded to it, when ``path[i+1]`` is bonded to ``path[j+1]``.

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


def check_chain(path: Sequence[Site], expected_sites: Iterable[Site]) -> None:
    """
    Check that a chain covers exactly the expected sites and is bonded throughout.

    Parameters
    ----------
    path : sequence of tuple of int
        The chain.
    expected_sites : iterable of tuple of int
        The sites it should visit, once each.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If a site is repeated or missing, or two consecutive sites aren't bonded.
    """
    if len(set(path)) != len(path) or set(path) != set(expected_sites):
        raise RuntimeError("the chain does not visit every site exactly once")
    for site_a, site_b in zip(path, path[1:]):
        if not are_bonded(site_a, site_b):
            raise RuntimeError("chain sites %s and %s are not bonded" % (site_a, site_b))


# ---------------------------------------------------------------------------
# growing the mobile chains
# ---------------------------------------------------------------------------
def grow_chain(start: Site, length: int, occupied: Set[Site], rng: random.Random,
               allowed: Callable[[Site], bool]) -> List[Site]:
    """
    Grow a self-avoiding chain from ``start`` into free space.

    Each new bead goes on a random free site bonded to the last one, picked
    with weight equal to the number of free sites around it (so the chain
    steers away from dead ends); if the chain does get stuck it backs up a
    few beads and tries again.

    Parameters
    ----------
    start : tuple of int
        The first bead's site (must be free).
    length : int
        Number of beads.
    occupied : set of tuple of int
        Sites already taken; the new chain's sites are added on success.
    rng : random.Random
        Source of randomness.
    allowed : callable
        Returns False for sites the chain may not use (besides occupied ones).

    Returns
    -------
    list of tuple of int
        The chain's sites.

    Raises
    ------
    RuntimeError
        If the chain can't be grown in ``MAX_GROWTH_TRIES`` back-ups.
    """
    path = [start]
    taken = set(occupied)
    taken.add(start)
    backups = 0

    def free(site: Site) -> bool:
        return site not in taken and in_box(site) and allowed(site)

    while len(path) < length:
        options = [add_sites(path[-1], step) for step in BOND_STEPS if free(add_sites(path[-1], step))]
        if not options:
            backups += 1
            if backups > MAX_GROWTH_TRIES or len(path) == 1:
                raise RuntimeError("could not grow a chain of %d from %s" % (length, start))
            for _ in range(min(3, len(path) - 1)):
                taken.discard(path.pop())
            continue
        weights = [1 + sum(1 for step in BOND_STEPS if free(add_sites(site, step))) for site in options]
        chosen = rng.choices(options, weights=weights)[0]
        path.append(chosen)
        taken.add(chosen)
    occupied.update(path)
    return path


def fg_sequence(motif: str, length: int) -> str:
    """
    The bead sequence of an FG domain: an anchor bead, then spacers with every third bead a motif.

    Parameters
    ----------
    motif : str
        The motif bead type (``FXFG`` or ``GLFG``).
    length : int
        Number of beads, anchor included.

    Returns
    -------
    str
        The sequence.
    """
    return FG_ANCHOR + "".join(motif if index % FG_MOTIF_SPACING == FG_MOTIF_SPACING - 1 else FG_SPACER
                               for index in range(length - 1))


def free_site_near(target: Tuple[float, float, float], occupied: Set[Site], touching: Set[Site]) -> Site:
    """
    The free site nearest ``target`` that is bonded to a site in ``touching``.

    Parameters
    ----------
    target : tuple of float
        Where the anchor should be, roughly.
    occupied : set of tuple of int
        Sites already taken.
    touching : set of tuple of int
        The anchor must be bonded to one of these (the scaffold it hangs from).

    Returns
    -------
    tuple of int
        The chosen site.

    Raises
    ------
    RuntimeError
        If there is no such site within 4 sites of the target.
    """
    tx, ty, tz = target
    best: List[Tuple[float, Site]] = []
    reach = 4
    for x in range(int(round(tx)) - reach, int(round(tx)) + reach + 1):
        for y in range(int(round(ty)) - reach, int(round(ty)) + reach + 1):
            for z in range(int(round(tz)) - reach, int(round(tz)) + reach + 1):
                site = (x, y, z)
                if site in occupied or not in_box(site):
                    continue
                if not any(add_sites(site, step) in touching for step in BOND_STEPS):
                    continue
                best.append(((x - tx) ** 2 + (y - ty) ** 2 + (z - tz) ** 2, site))
    if not best:
        raise RuntimeError("no free anchor site near %s" % (target,))
    return min(best)[1]


def fg_anchor_targets(where: str, spoke: float, copy: int, copies: int) -> Tuple[float, float, float]:
    """
    Where a given copy of an FG-Nup is anchored, roughly.

    Parameters
    ----------
    where : str
        "inner ring", "cytoplasmic ring", "cytoplasmic ring top" or "nuclear ring".
    spoke : float
        The spoke's angle, in degrees.
    copy : int
        Which copy on this spoke (0-based).
    copies : int
        How many copies this spoke carries.

    Returns
    -------
    tuple of float
        A target point in the local frame.
    """
    spoke_width = 360.0 / N_SPOKES
    if where == "inner ring":
        angle = spoke - spoke_width / 2 + (copy + 0.5) * spoke_width / copies
        z = -IR_HALF_HEIGHT + 2 + copy * (2 * IR_HALF_HEIGHT - 4) / max(copies - 1, 1)
        return polar_site(CHANNEL_RADIUS - 1.0, angle, z)
    if where == "cytoplasmic ring":
        return polar_site(Y_RING_RADII[0] - 4.5, spoke + NUP214_ANGLE_DEG, OUTER_RING_Z[1])
    if where == "cytoplasmic ring top":
        return polar_site(Y_RING_RADII[1], spoke + NUP358_ANGLE_DEG, OUTER_RING_Z[1] + 2 * DOMAIN_RADIUS + 0.5)
    if where == "nuclear ring":
        return polar_site(Y_RING_RADII[0] - 3.0, spoke + 5.0, -OUTER_RING_Z[1])
    raise ValueError("unknown anchor place %r" % where)


# ---------------------------------------------------------------------------
# checking the envelope seals the nucleus from the cytoplasm
# ---------------------------------------------------------------------------
def check_sealing(occupied: Set[Site]) -> None:
    """
    Check that the only way between cytoplasm and nucleus is the central channel.

    Empty sites are grouped into connected pockets (two empty sites connect if
    they're bonded, i.e. a bead could step from one to the other or a bond
    could span them). The cytoplasm (top face) must reach the nucleus (bottom
    face) through the open channel, must not reach it once the channel is
    plugged, and the lumen between the membranes must reach neither.

    Parameters
    ----------
    occupied : set of tuple of int
        Every frozen site (membrane and scaffold), local frame.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If the envelope leaks.
    """
    offset = (BOX_HALF_WIDTH, BOX_HALF_WIDTH, -BOX_Z[0])
    filled = np.zeros(BOX_DIMENSIONS, dtype=bool)
    for x, y, z in occupied:
        filled[x + offset[0], y + offset[1], z + offset[2]] = True
    connectivity = ndimage.generate_binary_structure(3, 3)

    def pockets(grid: np.ndarray) -> np.ndarray:
        labels, _ = ndimage.label(~grid, structure=connectivity)
        return np.asarray(labels)

    def label_at(labels: np.ndarray, site: Site) -> int:
        return int(labels[site[0] + offset[0], site[1] + offset[1], site[2] + offset[2]])

    corner = BOX_HALF_WIDTH - 1
    cytoplasm, nucleus, lumen = (corner, corner, BOX_Z[1]), (corner, corner, BOX_Z[0]), (corner, corner, 0)
    open_labels = pockets(filled)
    if label_at(open_labels, cytoplasm) != label_at(open_labels, nucleus):
        raise RuntimeError("the central channel does not connect the cytoplasm to the nucleus")
    if label_at(open_labels, lumen) in (label_at(open_labels, cytoplasm), 0):
        raise RuntimeError("the lumen leaks into the cytoplasm or nucleus")
    plugged = filled.copy()
    for x in range(-BOX_HALF_WIDTH, BOX_HALF_WIDTH + 1):
        for y in range(-BOX_HALF_WIDTH, BOX_HALF_WIDTH + 1):
            # the plug spans the channel and the spoke grooves cut into its wall
            if math.hypot(x, y) < CHANNEL_RADIUS + 2.5:
                plugged[x + offset[0], y + offset[1], offset[2] - 1:offset[2] + 2] = True
    plugged_labels = pockets(plugged)
    if label_at(plugged_labels, cytoplasm) == label_at(plugged_labels, nucleus):
        raise RuntimeError("the cytoplasm reaches the nucleus without passing through the central channel")


# ---------------------------------------------------------------------------
# writing the files
# ---------------------------------------------------------------------------
def write_parameter_file(filename: str) -> None:
    """
    Write the parameter file: every pair zero except those in ``INTERACTIONS``.

    Parameters
    ----------
    filename : str
        Where to write it.

    Returns
    -------
    None
    """
    lines = ["##",
             "## Nuclear pore complex demo parameter file (written by build_npc.py).",
             "##",
             "## Energies are in units where kT = TEMPERATURE (100 in KEYFILE.kf). FG motifs",
             "## (F = FxFG, Q = GLFG) are cohesive, transport receptors (K) bind them, the",
             "## inert cargo (D) binds nothing, and the membrane and scaffold are inert walls.",
             "##",
             "",
             "## ANGLE PENALTIES (every bead type needs one, even if zero)"]
    lines += ["ANGLE_PENALTY\t%s\t0\t0\t0" % bead for bead in ALL_BEAD_TYPES]
    lines += ["", "## PAIRWISE INTERACTIONS"]
    for index, bead_a in enumerate(ALL_BEAD_TYPES):
        for bead_b in ALL_BEAD_TYPES[index:]:
            energy = INTERACTIONS.get((bead_a, bead_b), INTERACTIONS.get((bead_b, bead_a), 0))
            lines.append("%s\t%s\t%d" % (bead_a, bead_b, energy))
    lines += ["", "## BEAD-SOLVENT INTERACTIONS (all zero)"]
    lines += ["%s\t0\t0" % bead for bead in ALL_BEAD_TYPES]
    with open(filename, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def to_box(site: Site) -> List[int]:
    """
    Convert a local-frame site to box coordinates.

    Parameters
    ----------
    site : tuple of int
        A site in the local frame.

    Returns
    -------
    list of int
        The same site in box coordinates (all >= 0).
    """
    return [site[0] + BOX_HALF_WIDTH, site[1] + BOX_HALF_WIDTH, site[2] - BOX_Z[0]]


def write_restart_and_freeze(chains: Dict[int, ChainRecord], frozen_ids: Sequence[int],
                             anchored_ids: Sequence[int]) -> None:
    """
    Write the restart file and the freeze file.

    Parameters
    ----------
    chains : dict
        chainID -> (local-frame sites in chain order, sequence, chainType).
    frozen_ids : sequence of int
        Chains frozen whole (``C`` line).
    anchored_ids : sequence of int
        Chains whose first bead is frozen (``B <id> 0`` lines).

    Returns
    -------
    None
    """
    restart = {
        "DIMENSIONS": list(BOX_DIMENSIONS),
        "ENERGY": 0,
        "HARDWALL": True,
        "CHAINS": {chain_id: [[to_box(site) for site in sites], sequence, chain_type]
                   for chain_id, (sites, sequence, chain_type) in chains.items()},
    }
    with open(RESTART_FILENAME, "wb") as handle:
        pickle.dump(restart, handle)
    with open(FREEZE_FILENAME, "w") as handle:
        handle.write("# the nuclear envelope and the NPC scaffold, frozen whole\n")
        handle.write("C %s\n" % " ".join(str(chain_id) for chain_id in frozen_ids))
        handle.write("# each FG-Nup is held by its first bead (bead 0), anchored on the scaffold\n")
        for chain_id in anchored_ids:
            handle.write("B %d 0\n" % chain_id)


# ---------------------------------------------------------------------------
# putting it together
# ---------------------------------------------------------------------------
def main() -> None:
    """
    Build the NPC demo and write ``npc.restart``, ``freeze.txt`` and ``params.prm``.

    Returns
    -------
    None
    """
    rng = random.Random(BUILD_SEED)
    cytoplasmic_ring = outer_ring_model(OUTER_RING_Z, CR, cytoplasmic=True)
    nuclear_ring = mirror_z(outer_ring_model(OUTER_RING_Z, NR, cytoplasmic=False),
                            {NR: NR, ELYS_DOMAIN: ELYS_DOMAIN})
    components: List[Tuple[str, Model]] = [
        ("nuclear envelope", nuclear_envelope_model()),
        ("inner ring", inner_ring_model()),
        ("cytoplasmic ring", cytoplasmic_ring),
        ("nuclear ring", nuclear_ring),
        ("nuclear basket", basket_model()),
        ("luminal ring", luminal_ring_model()),
    ]
    # later components give way to earlier ones where they overlap
    occupied: Set[Site] = set()
    chains: Dict[int, ChainRecord] = {}
    next_type = 0
    for name, model in components:
        model = {site: bead for site, bead in model.items() if site not in occupied and in_box(site)}
        path = thread_single_chain(model, rng)
        check_chain(path, model)
        chain_id = len(chains) + 1
        chains[chain_id] = (path, "".join(model[site] for site in path), next_type)
        next_type += 1
        occupied.update(path)
        print("%-18s chain %d, %6d beads (frozen)" % (name, chain_id, len(path)))
    frozen_ids = sorted(chains)
    check_sealing(occupied)
    print("the envelope seals: the central channel is the only way between cytoplasm and nucleus")

    scaffold_by_place = {
        "inner ring": {site for site in chains[2][0]},
        "cytoplasmic ring": {site for site, bead in cytoplasmic_ring.items() if bead == NUP214_DOMAIN},
        "cytoplasmic ring top": {site for site, bead in cytoplasmic_ring.items() if bead == NUP358_DOMAIN},
        "nuclear ring": set(chains[4][0]),
    }
    anchored_ids: List[int] = []
    for name, motif, length, where, copies in FG_NUPS:
        sequence = fg_sequence(motif, length)
        n_chains = 0
        for spoke in spoke_angles():
            for copy in range(copies):
                anchor = free_site_near(fg_anchor_targets(where, spoke, copy, copies), occupied,
                                        scaffold_by_place[where])
                path = grow_chain(anchor, length, occupied, rng, lambda site: True)
                chain_id = len(chains) + 1
                chains[chain_id] = (path, sequence, next_type)
                anchored_ids.append(chain_id)
                n_chains += 1
        print("%-18s %d chains of %d beads (bead 0 frozen on the %s)" % (name, n_chains, length, where))
        next_type += 1

    for name, bead, count in (("transport receptor", RECEPTOR, N_RECEPTORS), ("inert cargo", CARGO, N_CARGO)):
        for _ in range(count):
            while True:
                start = (rng.randint(-BOX_HALF_WIDTH, BOX_HALF_WIDTH), rng.randint(-BOX_HALF_WIDTH, BOX_HALF_WIDTH),
                         rng.randint(*SOLUBLE_START_Z))
                if start not in occupied:
                    break
            path = grow_chain(start, SOLUBLE_LENGTH, occupied, rng, lambda site: site[2] >= SOLUBLE_START_Z[0])
            chains[len(chains) + 1] = (path, bead * SOLUBLE_LENGTH, next_type)
        print("%-18s %d chains of %d beads (free, starting in the cytoplasm)" % (name, count, SOLUBLE_LENGTH))
        next_type += 1

    write_restart_and_freeze(chains, frozen_ids, anchored_ids)
    write_parameter_file(PARAMETER_FILENAME)
    n_beads = sum(len(sites) for sites, _, _ in chains.values())
    print("NPC: %d chains, %d beads in a %dx%dx%d box (%.0f x %.0f x %.0f nm)"
          % (len(chains), n_beads, *BOX_DIMENSIONS, *(edge * NM_PER_SITE for edge in BOX_DIMENSIONS)))
    print("wrote %s, %s and %s" % (RESTART_FILENAME, FREEZE_FILENAME, PARAMETER_FILENAME))


if __name__ == "__main__":
    main()
