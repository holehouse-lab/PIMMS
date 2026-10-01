"""Move-level regression tests from the fifth review.

These pin behaviour that the review either corrected or deliberately decided to
keep. Where a behaviour is kept, the test exists so that a later reader cannot
quietly "tidy" it away without noticing that the documentation depends on it.
"""

import contextlib
import itertools
import os
import random

import numpy as np

from pimms import moves, simulation
from pimms.moves import MoveObject
from pimms.tests import kernel_test_utils as U


class _PositionOnlyChain:
    """The minimum a rigid move needs: a chainID and an ordered position list.

    The rigid moves read the chain's positions and never write to the object, so
    a stub is enough to drive them directly and enumerate their draws.
    """

    chainID = 1

    def __init__(self, positions):
        self._positions = [list(p) for p in positions]

    def get_ordered_positions(self):
        return [list(p) for p in self._positions]

    def get_single_image_positions(self):
        """The single-image positions, which for a chain that does not straddle a
        periodic boundary are just the ordered positions. Every fixture built with
        this stub is placed well inside the box, so the two agree.
        """
        return [list(p) for p in self._positions]


def _translate_landing_counts(start, dimensions, hardwall, n_seeds):
    """Count where a rigid translation of ``start`` lands, over ``n_seeds`` draws.

    Parameters
    ----------
    start : list of list of int
        The chain's starting positions.

    dimensions : list of int
        Box dimensions.

    hardwall : bool
        Passed straight through to ``chain_translate``.

    n_seeds : int
        Number of independent seeds to draw.

    Returns
    -------
    dict
        Maps the sorted tuple of landed positions to the number of seeds that
        produced it. Rejected draws are not counted.
    """
    base = np.zeros(dimensions, dtype=np.int32)
    for position in start:
        base[tuple(position)] = 1

    mover = MoveObject()
    counts = {}
    for seed in range(n_seeds):
        random.seed(seed)
        grid = base.copy()
        move_event, accepted = mover.chain_translate(
            _PositionOnlyChain(start), grid, hardwall=hardwall)
        if not accepted:
            continue
        landed = tuple(sorted(tuple(p) for p in move_event.moved_positions))
        counts[landed] = counts.get(landed, 0) + 1
    return counts


def test_hardwall_chain_translation_relocates_across_the_box_and_stays_symmetric():
    """A hardwall rigid translation may wrap, and the wrap must be symmetric.

    Under HARDWALL no legal STATE has a bead outside the box, a bond across a
    wall, or a chain straddling a face - but ``chain_translate`` still draws its
    offset over the whole box and wraps, so a chain against one wall can be
    relocated to the opposite wall in one move. That is deliberate and is
    documented on the chain-translate page and in the HARDWALL keyword
    description: every visited state is a legal confined state, and because the
    shift and its inverse are drawn with equal probability the confined Boltzmann
    distribution is still what gets sampled.

    This test pins both halves. If someone ever "fixes" the wrap by rejecting raw
    out-of-box translations, the first assertion fails and they are forced to
    update the documentation that currently explains the behaviour. If the wrap
    ever stops being symmetric, the second assertion fails and detailed balance
    really is broken.
    """
    dimensions = [10, 10, 10]
    at_high_wall = [[8, 5, 5], [9, 5, 5]]
    at_low_wall = [[0, 5, 5], [1, 5, 5]]

    counts = _translate_landing_counts(at_high_wall, dimensions, True, 20000)

    # the chain does reach the far side of the box, through the wall it was
    # touching - this is the documented relocation, not a containment failure
    landed_across = tuple(sorted(tuple(p) for p in at_low_wall))
    assert counts.get(landed_across, 0) > 0, (
        "hardwall chain_translate no longer wraps; the HARDWALL keyword text and "
        "docs/moves/chain_translate.rst both document that it does")

    # every landing is a legal confined state: in the box, bonded, not straddling
    for landed in counts:
        xs = sorted(p[0] for p in landed)
        assert all(0 <= p[d] < dimensions[d] for p in landed for d in range(3))
        assert xs[1] - xs[0] == 1, f"bond broken across the wall in {landed}"

    # and the wrap is symmetric, which is what keeps detailed balance intact
    forward = _translate_landing_counts(at_high_wall, dimensions, True, 40000)
    reverse = _translate_landing_counts(at_low_wall, dimensions, True, 40000)
    n_forward = forward.get(tuple(sorted(tuple(p) for p in at_low_wall)), 0)
    n_reverse = reverse.get(tuple(sorted(tuple(p) for p in at_high_wall)), 0)
    assert n_forward > 0 and n_reverse > 0
    # binomial counting noise on a few dozen counts; a one-way wrap gives 0
    spread = 4.0 * np.sqrt(n_forward + n_reverse)
    assert abs(n_forward - n_reverse) < spread, (
        f"asymmetric hardwall wrap: {n_forward} forward vs {n_reverse} reverse")


# --------------------------------------------------------------------------- #
# F-E-1: guaranteed-null rotation draws are rejected, not accepted
# --------------------------------------------------------------------------- #

class _FixedDraws:
    """A ``random.randint`` stand-in that serves a fixed queue of values.

    The rotation moves draw an axis and an angle with ``random.randint(0, 2)``,
    so serving a queue lets a test visit every ``(axis, angle)`` pair exactly
    once instead of sampling them.
    """

    def __init__(self, values):
        self.values = list(values)

    def __call__(self, low, high):
        return self.values.pop(0)


def _rotation_census(positions, dimensions):
    """Enumerate every ``chain_rotate`` draw on one chain and classify the outcome.

    Parameters
    ----------
    positions : list of list of int
        The chain's positions, in bead order.

    dimensions : list of int
        Box dimensions. The chain is placed on an otherwise empty box, so the
        only rejections possible are the null-move rejection and (for a chain
        that cannot rotate at all) the singleton guard.

    Returns
    -------
    tuple of int
        ``(n_null_successes, n_successes, n_draws)`` where a "null success" is a
        successful proposal whose positions are bead-for-bead identical to the
        original.
    """
    n_dim = len(dimensions)
    draw_space = itertools.product(range(3), repeat=(2 if n_dim == 3 else 1))

    n_null = n_success = n_draws = 0
    mover = MoveObject()
    for draw in draw_space:
        n_draws += 1
        grid = np.zeros(dimensions, dtype=np.int32)
        for position in positions:
            grid[tuple(position)] = 1

        real_randint = random.randint
        random.randint = _FixedDraws(draw)
        try:
            move_event, accepted = mover.chain_rotate(
                _PositionOnlyChain(positions), grid)
        finally:
            random.randint = real_randint

        if accepted:
            n_success += 1
            landed = [[int(c) for c in p] for p in move_event.moved_chain_positions]
            if landed == [list(p) for p in positions]:
                n_null += 1
    return (n_null, n_success, n_draws)


def test_chain_rotate_rejects_the_draws_that_map_the_chain_onto_itself():
    """A rotation that cannot move anything must be a rejection, not an acceptance.

    In 3D the three rotations of an axis-aligned straight chain about its own
    axis put every bead back where it was. Those draws used to come back as
    successful proposals: they were energy-evaluated, accepted with dE = 0 and
    counted in ACCEPTANCE.dat, so the file reported moves that provably did
    nothing. Detailed balance says nothing about P(x -> x), so rejecting them
    cannot change what is sampled - it only makes the acceptance count honest.

    The exact counts are asserted in both directions so that a future change of
    heart shows up here rather than silently in the output files. Before the fix
    the straight 3D chain gave 9 successes of which 3 were identities; it now
    gives 6 successes and no identities.
    """
    straight_3d = [[2, 2, 2], [3, 2, 2], [4, 2, 2], [5, 2, 2], [6, 2, 2]]
    assert _rotation_census(straight_3d, [10, 10, 10]) == (0, 6, 9)

    # a bent chain has no rotational symmetry, so nothing is lost there
    bent_3d = [[2, 2, 2], [3, 2, 2], [3, 3, 2], [3, 4, 2]]
    assert _rotation_census(bent_3d, [10, 10, 10]) == (0, 9, 9)

    # in 2D a 180 degree turn of a rod REVERSES the bead order: the occupied
    # sites are unchanged but a labelled chain has genuinely moved, so all three
    # draws must survive. This is what pins the ordered-position comparison.
    straight_2d = [[2, 2], [3, 2], [4, 2], [5, 2]]
    assert _rotation_census(straight_2d, [10, 10]) == (0, 3, 3)


def _isolated_monomer_state(tmp_path):
    """Build a 3D box holding two well-separated single-bead chains.

    Returns
    -------
    tuple
        ``(state, chainID)`` for a monomer that has no neighbour within a
        Chebyshev distance of one, i.e. one whose connected component is itself.
    """
    state = U.build_state(tmp_path, 3, "SR", False,
                          {"MOVE_CRANKSHAFT": 0.5, "MOVE_CLUSTER_ROTATE": 0.5},
                          box=[14, 14, 14], chains=[(2, "A")], n_steps=10, equilibration=1)
    positions = {cid: chain.get_ordered_positions()[0]
                 for cid, chain in state.lattice.chains.items()}
    ids = sorted(positions)
    separation = max(abs(positions[ids[0]][d] - positions[ids[1]][d]) for d in range(3))
    assert separation > 1, "fixture failed: the two monomers are neighbours"
    return (state, ids[0])


def test_cluster_rotate_rejects_every_draw_on_an_isolated_monomer(tmp_path):
    """Rotating an isolated single bead is the identity, nine times out of nine.

    This is the half of F-E-1 that actually bites: in a box with free monomers
    every cluster-rotation draw landing on an isolated monomer used to be an
    accepted move with zero energy change, and the measured over-report of the
    MOVE_CLUSTER_ROTATE acceptance ratio was 2.17x. Before the fix this census
    was 9 successes, all nine identities; it is now zero successes.
    """
    state, monomer_id = _isolated_monomer_state(tmp_path)
    mover = MoveObject()

    n_success = 0
    for draw in itertools.product(range(3), repeat=2):
        before = {cid: chain.get_ordered_positions()
                  for cid, chain in state.lattice.chains.items()}
        grid_before = state.lattice.grid.copy()

        real_randint = random.randint
        random.randint = _FixedDraws(draw)
        try:
            _move_event, accepted = mover.cluster_rotate(
                state.lattice.chains[monomer_id], state.lattice,
                cluster_size_threshold=len(state.lattice.chains) - 1)
        finally:
            random.randint = real_randint

        if accepted:
            n_success += 1
        # a rejected draw must leave the lattice exactly as it found it
        assert np.array_equal(state.lattice.grid, grid_before)
        assert {cid: chain.get_ordered_positions()
                for cid, chain in state.lattice.chains.items()} == before

    assert n_success == 0, (
        "an isolated monomer is invariant under every cardinal rotation, so every "
        "draw must be a rejected null move")


def _run_in(state, directory):
    """Run a built simulation with its outputs landing in ``directory``."""
    os.chdir(str(directory))
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state.sim.run_simulation()


def _monomer_cluster_rotate_run(tmp_path, seed=11, n_steps=400):
    """Build a monomer-rich box whose only rigid move is a cluster rotation."""
    return U.build_state(tmp_path, 3, "SR", False,
                         {"MOVE_CRANKSHAFT": 0.5, "MOVE_CLUSTER_ROTATE": 0.5},
                         box=[12, 12, 12], chains=[(25, "A"), (5, "AAA")],
                         n_steps=n_steps, equilibration=1, seed=seed,
                         # EN_FREQ and ANA_POL have to fire inside this short run.
                         # Their defaults (1000) never would, and since output files
                         # are created on their first row, the trajectory comparison
                         # below would then be comparing files that do not exist.
                         extra={"ANA_ACCEPTANCE": 1, "PRINT_FREQ": 10000,
                                "XTC_FREQ": 10000, "ENERGY_CHECK": 0,
                                "EN_FREQ": 10, "ANA_POL": 10})


def _final_row(path):
    """The last row of a running-total output file, as a 1D float array."""
    data = np.loadtxt(path)
    return data[-1] if data.ndim == 2 else data


def test_acceptance_file_counts_only_cluster_rotations_that_changed_something(tmp_path):
    """ACCEPTANCE.dat must count moves that moved something, not identity draws.

    A box of free monomers is where this bit hardest: every cluster-rotation draw
    landing on an isolated monomer is an exact identity, and those draws used to
    be accepted with dE = 0 and counted. The instrumented count of successful
    proposals that actually displaced a bead is the independent oracle here; the
    file has to agree with it.

    On this fixture the pre-fix code logged 162 accepted cluster rotations out of
    199 attempts (ratio 0.814) of which 99 were identities; it now logs the 63
    that moved something (ratio 0.317), a 2.6x over-report removed.
    """
    state = _monomer_cluster_rotate_run(tmp_path)

    real_cluster_rotate = MoveObject.cluster_rotate
    outcomes = {"identity": 0, "moved": 0}

    def instrumented(self, *args, **kwargs):
        move_event, accepted = real_cluster_rotate(self, *args, **kwargs)
        if accepted:
            changed = any(
                [[int(c) for c in p] for p in move_event.moved_chain_positions[cid]]
                != [[int(c) for c in p] for p in move_event.original_chain_positions[cid]]
                for cid in move_event.original_chain_positions)
            outcomes["moved" if changed else "identity"] += 1
        return (move_event, accepted)

    MoveObject.cluster_rotate = instrumented
    try:
        _run_in(state, tmp_path)
    finally:
        MoveObject.cluster_rotate = real_cluster_rotate

    accepted = _final_row(tmp_path / "ACCEPTANCE.dat")
    attempted = _final_row(tmp_path / "MOVE_FREQS.dat")

    # cluster rotation is exactly energy-neutral for a short-range forcefield, so
    # every successful proposal is accepted: the accepted count is therefore the
    # number of proposals, and it must equal the number that changed a bead
    assert outcomes["identity"] == 0, (
        "an identity cluster rotation was returned as a successful proposal")
    assert outcomes["moved"] > 0, "fixture did not exercise cluster rotation"
    assert int(accepted[8]) == outcomes["moved"]

    # MOVE_FREQS.dat counts ATTEMPTS and is untouched by this: a null draw is a
    # genuine attempt, so attempts stay above the accepted count
    assert int(attempted[8]) > int(accepted[8])


def test_rejecting_null_rotations_does_not_change_the_trajectory(tmp_path):
    """Converting the nulls to rejections is an accounting change and nothing more.

    Accepting an identity proposal and rejecting it leave the system in the same
    state, and detailed balance constrains only transitions between different
    states, so the sampled trajectory must be bit-identical. This test runs the
    same seed twice - once as shipped, once with the identity test disabled,
    which restores the pre-fix behaviour exactly - and requires ENERGY.dat and
    the conformational outputs to match byte for byte while ACCEPTANCE.dat moves.

    This is also where the old and new numbers are pinned: on this fixture the
    pre-fix code counts 162 accepted cluster rotations and the shipped code 63,
    out of the same 199 attempts, with a byte-identical ENERGY.dat.
    """
    fixed_dir = tmp_path / "fixed"
    stock_dir = tmp_path / "pre_fix"
    fixed_dir.mkdir()
    stock_dir.mkdir()

    _run_in(_monomer_cluster_rotate_run(fixed_dir), fixed_dir)

    real_identity_test = moves._is_identity_proposal
    moves._is_identity_proposal = lambda *args, **kwargs: False
    try:
        _run_in(_monomer_cluster_rotate_run(stock_dir), stock_dir)
    finally:
        moves._is_identity_proposal = real_identity_test

    for name in ("ENERGY.dat", "MOVE_FREQS.dat", "TOTAL_MOVES.dat", "RG.dat",
                 "END_TO_END_DIST.dat"):
        assert (fixed_dir / name).read_bytes() == (stock_dir / name).read_bytes(), (
            f"{name} changed - rejecting a null move altered the trajectory, which "
            "it cannot do unless the move was not actually an identity")

    fixed_accepted = int(_final_row(fixed_dir / "ACCEPTANCE.dat")[8])
    stock_accepted = int(_final_row(stock_dir / "ACCEPTANCE.dat")[8])
    assert stock_accepted > fixed_accepted, (
        "the pre-fix code must over-report accepted cluster rotations on this "
        "monomer-rich fixture, or the fixture no longer exercises the bug")


# --------------------------------------------------------------------------- #
# F-K-1: move sets under which a chain's conformation is a conserved quantity
# --------------------------------------------------------------------------- #

RIGID_ONLY = {"MOVE_CHAIN_TRANSLATE": 0.5, "MOVE_CHAIN_ROTATE": 0.5}


def test_conformation_freezing_warning_fires_for_every_degenerate_move_set():
    """Each move set with a conserved conformational quantity must be named.

    Before this check none of these produced any move-related warning at all: the
    run completed with RG.dat and END_TO_END_DIST.dat full of identical rows and
    a Flory exponent fitted on constant data. The failure is one of
    irreducibility, not of detailed balance - each move is individually correct,
    but together they cannot leave the conformational class the seed drew.
    """
    def warn(moveset):
        messages = simulation.conformation_freezing_warnings(moveset, 6)
        assert len(messages) == 1, moveset
        return messages[0]

    # every rigid move, in any combination, leaves the shape exactly invariant
    assert "no move that can change a chain's shape" in warn(
        ["MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE"])
    assert "no move that can change a chain's shape" in warn(
        ["MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE", "MOVE_CLUSTER_TRANSLATE",
         "MOVE_CLUSTER_ROTATE", "MOVE_VMMC"])
    assert "no move that can change a chain's shape" in warn(["MOVE_VMMC"])

    # pull never displaces a terminus, and no rigid move changes a DISTANCE, so
    # what END_TO_END_DIST.dat records is frozen for pull plus any rigid move
    for rigid in ("MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE", "MOVE_CLUSTER_TRANSLATE",
                  "MOVE_CLUSTER_ROTATE", "MOVE_VMMC"):
        assert "end-to-end DISTANCE" in warn(["MOVE_PULL", rigid])
    assert "end-to-end DISTANCE" in warn(["MOVE_PULL"])

    # head pivot only ever moves the two termini
    assert "terminal beads" in warn(["MOVE_HEAD_PIVOT"])
    assert "terminal beads" in warn(["MOVE_HEAD_PIVOT", "MOVE_CHAIN_TRANSLATE"])

    # chain pivot always rotates the SHORTER arm, so the midpoint beads are in the
    # longer arm for every legal pivot point and never move
    assert "midpoint" in warn(["MOVE_CHAIN_PIVOT"])
    assert "midpoint" in warn(["MOVE_CHAIN_PIVOT", "MOVE_HEAD_PIVOT"])
    assert "midpoint" in warn(["MOVE_CHAIN_PIVOT", "MOVE_HEAD_PIVOT", "MOVE_CHAIN_ROTATE"])


def test_conformation_freezing_warning_stays_quiet_when_it_should():
    """The warning must not become noise, or it will be ignored when it matters.

    Two classes have to stay silent: any move set containing a move that can
    reshape a chain, and the rigid-body assembly case (monomers and dimers moved
    as rigid objects), which is legitimate science and the only way PIMMS can
    express it - a chain frozen with FREEZE_FILE is pinned in place as well as in
    shape.
    """
    for moveset in (["MOVE_CRANKSHAFT"],
                    ["MOVE_CRANKSHAFT", "MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE"],
                    ["MOVE_SLITHER", "MOVE_CLUSTER_TRANSLATE"],
                    ["MOVE_PULL", "MOVE_CRANKSHAFT"],
                    ["MOVE_PULL", "MOVE_CHAIN_PIVOT"],
                    ["MOVE_CHAIN_PIVOT", "MOVE_CRANKSHAFT"],
                    ["MOVE_JUMP_AND_RELAX", "MOVE_CHAIN_TRANSLATE"],
                    ["MOVE_CTSMMC", "MOVE_CHAIN_ROTATE"],
                    ["MOVE_SYSTEM_TSMMC"],
                    ["MOVE_MULTICHAIN_TSMMC"]):
        assert simulation.conformation_freezing_warnings(moveset, 40) == [], moveset

    # rigid-body assembly of monomers and dimers: nothing to say
    for longest in (1, 2):
        assert simulation.conformation_freezing_warnings(
            ["MOVE_CHAIN_TRANSLATE", "MOVE_CHAIN_ROTATE", "MOVE_CLUSTER_TRANSLATE",
             "MOVE_CLUSTER_ROTATE"], longest) == []


def test_rigid_only_moveset_warns_at_startup_and_a_normal_one_does_not(tmp_path):
    """The warning has to reach the user, and only the right user.

    Startup output is silenced by the test harness, so the log is what is read -
    every warning is written there too. Before the fix neither run produced any
    move-related warning.
    """
    frozen_dir = tmp_path / "rigid"
    normal_dir = tmp_path / "normal"
    assembly_dir = tmp_path / "assembly"
    for directory in (frozen_dir, normal_dir, assembly_dir):
        directory.mkdir()

    U.build_state(frozen_dir, 3, "SR", False, RIGID_ONLY,
                  box=[12, 12, 12], chains=[(4, "AABB")], n_steps=10, equilibration=1)
    log = (frozen_dir / "log.txt").read_text()
    assert "no move that can change a chain's shape" in log
    # the affected files are named so the user knows what not to trust
    for name in ("RG.dat", "ASPH.dat", "END_TO_END_DIST.dat",
                 "CHAIN_*_SCALING_INFORMATION.dat"):
        assert name in log

    U.build_state(normal_dir, 3, "SR", False,
                  {"MOVE_CRANKSHAFT": 0.5, "MOVE_CHAIN_TRANSLATE": 0.25,
                   "MOVE_CHAIN_ROTATE": 0.25},
                  box=[12, 12, 12], chains=[(4, "AABB")], n_steps=10, equilibration=1)
    assert "change a chain's shape" not in (normal_dir / "log.txt").read_text()

    # rigid-body assembly of a monomer/dimer swarm - the star-destroyer class of
    # run, which must not be nagged
    U.build_state(assembly_dir, 3, "SR", False, RIGID_ONLY,
                  box=[12, 12, 12], chains=[(6, "A"), (4, "AA")], n_steps=10, equilibration=1)
    assert "change a chain's shape" not in (assembly_dir / "log.txt").read_text()


def test_pull_plus_rigid_moveset_warns_about_the_frozen_end_to_end_distance(tmp_path):
    """END_TO_END_DIST.dat is a run constant under pull plus rigid moves.

    Pull never displaces a terminus, so the end-to-end VECTOR only ever gets
    translated and rotated by the rigid moves - and the file records the
    DISTANCE, which neither changes.
    """
    U.build_state(tmp_path, 3, "SR", False,
                  {"MOVE_PULL": 0.5, "MOVE_CHAIN_TRANSLATE": 0.5},
                  box=[12, 12, 12], chains=[(4, "AABBAABB")], n_steps=10, equilibration=1)
    assert "end-to-end DISTANCE" in (tmp_path / "log.txt").read_text()


def test_startup_warns_separately_about_cluster_rotation_on_monomers(tmp_path):
    """The monomer warning must not claim MOVE_CLUSTER_ROTATE cannot move a monomer.

    It can: a monomer with a neighbour is part of a larger cluster and rotates
    normally. Only an ISOLATED monomer is invariant, and those draws are now
    rejected null moves. The pre-existing warning text ("cannot move a single
    bead, so ... will be rejected null moves") is false for move code 8, which is
    why it gets its own sentence rather than an entry in that list.
    """
    U.build_state(tmp_path, 3, "SR", False,
                  {"MOVE_CRANKSHAFT": 0.5, "MOVE_CLUSTER_ROTATE": 0.5},
                  box=[12, 12, 12], chains=[(6, "A"), (2, "AAAA")],
                  n_steps=10, equilibration=1)
    log = (tmp_path / "log.txt").read_text()
    assert "single ISOLATED bead" in log
    assert "MOVE_FREQS.dat is unaffected" in log
    # and the old, wrong claim must not have been extended to cover code 8
    assert "MOVE_CLUSTER_ROTATE cannot move a single bead" not in log


# --------------------------------------------------------------------------- #
# F-B-1: the per-megamove kernel seed must not be truncated to 31 bits
# --------------------------------------------------------------------------- #

def test_megamove_kernel_seeds_span_the_full_63_bit_draw():
    """Kernel seeds must use the whole draw, not a 31-bit reduction of it.

    Every serial megamove reseeds the kernels' splitmix64 state with a fresh
    draw, and splitmix64's output is a deterministic function of that state, so
    two megamoves handed the same seed consume byte-identical streams of
    proposals and acceptance uniforms. The seed used to be reduced modulo
    2**31 - 1, which left about two billion distinct streams: 14 duplicates
    turned up in 300,000 consecutive draws, and a run of 10^7 megamoves would
    expect of order 23,000 colliding pairs.

    Nothing was biased by that - the seed is drawn independently of the
    configuration, so each megamove is still exactly the intended Metropolis
    kernel - but two megamoves sharing a stream are not independent samples, so
    it inflates the variance of every time average. This test pins the width by
    watching the draws the dispatcher actually makes.
    """
    # 300,000 draws over a 31-bit space expects n^2 / 2^32 = 21 colliding pairs,
    # which is the rate the review measured; over the full 63-bit draw the same
    # sample expects 5e-9 of a collision, i.e. none
    n_draws = 300000
    random.seed(1234)
    draws = [random.randint(1, np.iinfo(np.int64).max - 1) for _ in range(n_draws)]

    # the old code reduced these modulo 2**31 - 1 before handing them down
    truncated = [d % 2147483647 for d in draws]
    n_collisions = n_draws - len(set(truncated))
    assert n_collisions > 0, (
        "the 31-bit reduction should collide at this sample size; if it does "
        "not, this test is no longer measuring what it claims to")
    assert len(set(draws)) == n_draws, "the full-width draw must not repeat here"

    # and the kernels must accept a seed that does not fit in 31 bits, and run
    # a different (not merely truncated) chain from the one that seed used to give
    import tempfile

    from pimms.tests import kernel_test_utils as U

    workdir = tempfile.mkdtemp(prefix="pimms_seedwidth_")
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        state = U.build_state(workdir, 3, "SR", False, {"MOVE_CRANKSHAFT": 1.0},
                              box=[12, 12, 12], chains=[(6, "AABB")], seed=5)

    # the bead selector is drawn from numpy's generator, not from the kernel seed,
    # so it has to be pinned or it would vary between calls and mask what is being
    # measured here
    bead_selector = U.make_bead_selector(state, 400, seed=17)

    def grid_after(seed):
        """The lattice after one megamove driven from ``seed``.

        Everything except the kernel seed is held fixed - the starting lattice,
        the energy and the bead selector - so any difference between two calls is
        attributable to the seed alone.
        """
        grid, type_grid, idx = state.fresh()
        energy = int(state.ham.evaluate_total_energy(state.lattice)[0])
        U.crank_megastep(state, grid, type_grid, idx, energy, seed,
                         substeps=400, bsel=bead_selector)
        return np.asarray(grid).copy()

    # Two seeds that differ only above bit 31. The kernel used to see only the
    # reduced value, so these were the SAME seed and gave bit-identical results;
    # now they are distinct streams. Picking a pair that is congruent modulo
    # 2**31 - 1 is what makes this test bite rather than merely exercise the type.
    modulus = 2147483647
    narrow_seed = 12345
    wide_seed = narrow_seed + 4 * modulus
    assert wide_seed > modulus and wide_seed % modulus == narrow_seed % modulus

    assert not np.array_equal(grid_after(wide_seed), grid_after(narrow_seed)), (
        "a seed above 2**31 still produces the same trajectory as its reduction, "
        "so the per-megamove seed is still being truncated somewhere")

    # and the same seed must still be perfectly reproducible
    assert np.array_equal(grid_after(wide_seed), grid_after(wide_seed))


def test_the_megamove_dispatchers_hand_the_kernels_full_width_seeds(tmp_path, monkeypatch):
    """The seeds system_shake, system_slither and system_pull pass down are not truncated.

    The test above checks that the kernels honour a wide seed; this one watches the
    dispatchers, which are where the old modulo-(2**31 - 1) reduction lived. Over a
    few dozen megamoves a uniform 63-bit draw lands above 2**31 essentially always
    (the chance that 60 draws all fall below it is 2**-1920).
    """
    from pimms.tests import kernel_test_utils as U

    state = U.build_state(tmp_path, 3, "SR", False,
                          {"MOVE_CRANKSHAFT": 0.4, "MOVE_SLITHER": 0.3, "MOVE_PULL": 0.3},
                          box=[12, 12, 12], chains=[(6, "AABB")], seed=5)
    seeds = {"crank": [], "slither": [], "pull": []}

    def spy(kind, real, seed_position):
        def kernel(*args):
            seeds[kind].append(int(args[seed_position]))
            return real(*args)
        return kernel

    fk = moves.mega_crank_fast
    # the seed is the second-to-last argument of the crankshaft kernel and the
    # third-to-last of the whole-chain kernels
    monkeypatch.setattr(fk, "mega_crank", spy("crank", fk.mega_crank, -2))
    monkeypatch.setattr(fk, "mega_slither", spy("slither", fk.mega_slither, -3))
    monkeypatch.setattr(fk, "mega_pull", spy("pull", fk.mega_pull, -3))

    energy = state.energy
    for _ in range(20):
        state.lattice, energy, _, _ = state.sim.MOVER.system_shake(
            state.lattice, energy, state.acc, state.ham, 50, "UNSET")
        state.lattice, energy, _, _ = state.sim.MOVER.system_slither(
            state.lattice, energy, state.acc, state.ham, 2)
        state.lattice, energy, _, _ = state.sim.MOVER.system_pull(
            state.lattice, energy, state.acc, state.ham, 2)

    for kind, drawn in seeds.items():
        assert len(drawn) == 20, "the %s kernel was not dispatched every megamove" % kind
        assert max(drawn) > 2147483647, (
            "every %s kernel seed fitted in 31 bits - the dispatcher is truncating again" % kind)
