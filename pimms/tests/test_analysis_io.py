from types import SimpleNamespace

import numpy as np
import pytest

from pimms import analysis_IO


@pytest.fixture
def cfg_paths(monkeypatch, tmp_path):
    mapping = {
        "OUTNAME_ENERGY": "energy.out",
        "OUTNAME_CLUSTERS": "clusters.out",
        "OUTNAME_NUM_CLUSTERS": "num_clusters.out",
        "OUTNAME_LR_CLUSTERS": "lr_clusters.out",
        "OUTNAME_NUM_LR_CLUSTERS": "num_lr_clusters.out",
        "OUTNAME_CLUSTER_RG": "cluster_rg.out",
        "OUTNAME_CLUSTER_ASPH": "cluster_asph.out",
        "OUTNAME_CLUSTER_VOL": "cluster_vol.out",
        "OUTNAME_CLUSTER_AREA": "cluster_area.out",
        "OUTNAME_CLUSTER_DENSITY": "cluster_density.out",
        "OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE": "cluster_radial_density.out",
        "OUTNAME_LR_CLUSTER_RG": "lr_cluster_rg.out",
        "OUTNAME_LR_CLUSTER_ASPH": "lr_cluster_asph.out",
        "OUTNAME_LR_CLUSTER_VOL": "lr_cluster_vol.out",
        "OUTNAME_LR_CLUSTER_AREA": "lr_cluster_area.out",
        "OUTNAME_LR_CLUSTER_DENSITY": "lr_cluster_density.out",
        "OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE": "lr_cluster_radial_density.out",
        "OUTNAME_INTERNAL_SCALING": "internal_scaling.out",
        "OUTNAME_INTERNAL_SCALING_SQUARED": "internal_scaling_squared.out",
        "OUTNAME_SCALING_INFORMATION": "scaling_information.out",
        "OUTNAME_DMAP": "dmap.out",
        "OUTNAME_RG": "rg.out",
        "OUTNAME_ASPH": "asph.out",
        "OUTNAME_E2E": "e2e.out",
        "OUTNAME_R2R": "r2r.out",
        "OUTNAME_MOVES": "moves.out",
        "OUTNAME_ACCEPTANCE": "acceptance.out",
        "OUTNAME_TOTAL_MOVES": "total_moves.out",
        "OUTNAME_PERFORMANCE": "performance.out",
        "QUENCHFILE_NAME": "quench.out",
    }

    out = {}
    for attr, fname in mapping.items():
        path = tmp_path / fname
        monkeypatch.setattr(analysis_IO.CONFIG, attr, str(path))
        out[attr] = path

    return out


def _read(path):
    return path.read_text()


def test_write_energy_appends(cfg_paths):
    analysis_IO.write_energy(1, -12.5)
    analysis_IO.write_energy(2, -11.0)

    assert _read(cfg_paths["OUTNAME_ENERGY"]) == "1\t  -12.5000\n2\t  -11.0000\n"


def test_write_clusters_writes_counts_and_type_fractions(cfg_paths):
    clusters = [[1, 2], [3]]
    id_to_type = {1: 10, 2: 20, 3: 10}

    analysis_IO.write_clusters(5, clusters, id_to_type)

    assert _read(cfg_paths["OUTNAME_CLUSTERS"]) == "5, 2, 1, \n"
    assert _read(cfg_paths["OUTNAME_NUM_CLUSTERS"]) == "5\t2\n"

    chain10 = cfg_paths["OUTNAME_CLUSTERS"].parent / f"CHAIN_{10}_{cfg_paths['OUTNAME_CLUSTERS'].name}"
    chain20 = cfg_paths["OUTNAME_CLUSTERS"].parent / f"CHAIN_{20}_{cfg_paths['OUTNAME_CLUSTERS'].name}"

    assert chain10.read_text() == "0.5000, 1.0000, \n"
    assert chain20.read_text() == "0.5000, 0.0000, \n"


def test_write_clusters_single_type_skips_composition_files(cfg_paths):
    clusters = [[1], [2]]
    id_to_type = {1: 7, 2: 7}

    analysis_IO.write_clusters(1, clusters, id_to_type)

    chain_file = cfg_paths["OUTNAME_CLUSTERS"].parent / f"CHAIN_{7}_{cfg_paths['OUTNAME_CLUSTERS'].name}"
    assert not chain_file.exists()


def test_write_clusters_handles_empty_cluster_without_division_error(cfg_paths):
    clusters = [[]]
    id_to_type = {1: 1, 2: 2}

    analysis_IO.write_clusters(2, clusters, id_to_type)

    chain1 = cfg_paths["OUTNAME_CLUSTERS"].parent / f"CHAIN_{1}_{cfg_paths['OUTNAME_CLUSTERS'].name}"
    chain2 = cfg_paths["OUTNAME_CLUSTERS"].parent / f"CHAIN_{2}_{cfg_paths['OUTNAME_CLUSTERS'].name}"

    assert chain1.read_text() == "0.0000, \n"
    assert chain2.read_text() == "0.0000, \n"


def test_write_lr_clusters_writes_counts_and_type_fractions(cfg_paths):
    clusters = [[1, 2], [3]]
    id_to_type = {1: 10, 2: 20, 3: 10}

    analysis_IO.write_LR_clusters(5, clusters, id_to_type)

    assert _read(cfg_paths["OUTNAME_LR_CLUSTERS"]) == "5, 2, 1, \n"
    assert _read(cfg_paths["OUTNAME_NUM_LR_CLUSTERS"]) == "5\t2\n"

    chain10 = cfg_paths["OUTNAME_LR_CLUSTERS"].parent / f"CHAIN_{10}_{cfg_paths['OUTNAME_LR_CLUSTERS'].name}"
    chain20 = cfg_paths["OUTNAME_LR_CLUSTERS"].parent / f"CHAIN_{20}_{cfg_paths['OUTNAME_LR_CLUSTERS'].name}"

    assert chain10.read_text() == "0.5000, 1.0000, \n"
    assert chain20.read_text() == "0.5000, 0.0000, \n"


def test_write_cluster_properties_writes_all_outputs(cfg_paths):
    analysis_IO.write_cluster_properties(
        step=4,
        cluster_polymeric_properties_list=[(1.2, 0.1), (2.3, 0.2)],
        cluster_size_list=[(10.0, 8.0, 0.8), (15.0, 11.0, 0.7)],
        cluster_radial_density=[[0.1, 0.2], [0.3]],
    )

    assert _read(cfg_paths["OUTNAME_CLUSTER_RG"]) == "4, 1.2000, 2.3000, \n"
    assert _read(cfg_paths["OUTNAME_CLUSTER_ASPH"]) == "4, 0.1000, 0.2000, \n"
    assert _read(cfg_paths["OUTNAME_CLUSTER_VOL"]) == "4, 10.0000, 15.0000, \n"
    assert _read(cfg_paths["OUTNAME_CLUSTER_AREA"]) == "4, 8.0000, 11.0000, \n"
    assert _read(cfg_paths["OUTNAME_CLUSTER_DENSITY"]) == "4, 0.8000, 0.7000, \n"
    assert _read(cfg_paths["OUTNAME_CLUSTER_RADIAL_DENSITY_PROFILE"]) == "4, C1, 0.1000, 0.2000, \n4, C2, 0.3000, \n"


def test_write_lr_cluster_properties_writes_all_outputs(cfg_paths):
    analysis_IO.write_LR_cluster_properties(
        step=4,
        LR_cluster_polymeric_properties_list=[(1.2, 0.1)],
        LR_cluster_size_list=[(10.0, 8.0, 0.8)],
        LR_cluster_radial_density=[[0.1]],
    )

    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_RG"]) == "4, 1.2000, \n"
    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_ASPH"]) == "4, 0.1000, \n"
    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_VOL"]) == "4, 10.0000, \n"
    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_AREA"]) == "4, 8.0000, \n"
    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_DENSITY"]) == "4, 0.8000, \n"
    assert _read(cfg_paths["OUTNAME_LR_CLUSTER_RADIAL_DENSITY_PROFILE"]) == "4, C1, 0.1000, \n"


def test_write_internal_scaling_and_prefix(cfg_paths):
    analysis_IO.write_internal_scaling([1.1, 2.2], [3.3])

    assert _read(cfg_paths["OUTNAME_INTERNAL_SCALING"]) == "1\t1.1000\n2\t2.2000\n"
    assert _read(cfg_paths["OUTNAME_INTERNAL_SCALING_SQUARED"]) == "1\t3.3000\n"

    analysis_IO.write_internal_scaling([9.9], [8.8], prefix="pre_")
    pref_is = cfg_paths["OUTNAME_INTERNAL_SCALING"].parent / f"pre_{cfg_paths['OUTNAME_INTERNAL_SCALING'].name}"
    pref_is2 = cfg_paths["OUTNAME_INTERNAL_SCALING"].parent / f"pre_{cfg_paths['OUTNAME_INTERNAL_SCALING_SQUARED'].name}"

    assert pref_is.read_text() == "1\t9.9000\n"
    assert pref_is2.read_text() == "1\t8.8000\n"


def test_write_scaling_information_and_length_validation(cfg_paths):
    analysis_IO.write_scaling_information([0.1, 0.2], [1.0, 2.0])
    assert _read(cfg_paths["OUTNAME_SCALING_INFORMATION"]) == "0.1000\t1.0000\n0.2000\t2.0000\n"

    with pytest.raises(ValueError, match="same length"):
        analysis_IO.write_scaling_information([0.1], [1.0, 2.0])


def test_write_distance_map(cfg_paths):
    dmap = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
    analysis_IO.write_distance_map(dmap)

    assert _read(cfg_paths["OUTNAME_DMAP"]) == "1.0000\t2.0000\t\n3.0000\t4.0000\t\n"


def test_write_rg_asph_e2e(cfg_paths):
    analysis_IO.write_radius_of_gyration(1, [1.2345, 2.0])
    analysis_IO.write_asphericity(1, [0.1234])
    analysis_IO.write_end_to_end(1, [9.8765])

    assert _read(cfg_paths["OUTNAME_RG"]) == "1\t1.234\t2.000\t\n"
    assert _read(cfg_paths["OUTNAME_ASPH"]) == "1\t0.123\t\n"
    assert _read(cfg_paths["OUTNAME_E2E"]) == "1\t9.877\t\n"


def test_write_residue_residue_distance_and_length_validation(cfg_paths):
    analysis_IO.write_residue_residue_distance(7, [(1, 3), (2, 4)], [[1.0], [2.0, 3.0]])

    assert _read(cfg_paths["OUTNAME_R2R"]) == "7\t1\t3\t1.000\t\n7\t2\t4\t2.000\t3.000\t\n"

    with pytest.raises(ValueError, match="same length"):
        analysis_IO.write_residue_residue_distance(8, [(1, 3)], [[1.0], [2.0]])


def test_write_acceptance_statistics_and_total_moves(cfg_paths):
    acceptance = SimpleNamespace(
        move_count=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        accepted_count=[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        alt_Markov_chain_moves=5,
    )

    analysis_IO.write_acceptance_statistics(20, acceptance)

    assert _read(cfg_paths["OUTNAME_MOVES"]) == "20\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t\n"
    assert _read(cfg_paths["OUTNAME_ACCEPTANCE"]) == "20\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t\n"

    # total uses moves [1..8,12,13] + alt_Markov_chain_moves
    expected_total = sum([1, 2, 3, 4, 5, 6, 7, 8, 12, 13]) + 5
    assert _read(cfg_paths["OUTNAME_TOTAL_MOVES"]) == f"20\t{expected_total}\n"


def test_write_acceptance_statistics_raises_when_move_count_length_changes(cfg_paths):
    acceptance = SimpleNamespace(
        move_count=[0] * 13,
        accepted_count=[0] * 13,
        alt_Markov_chain_moves=0,
    )

    with pytest.raises(analysis_IO.AcceptanceException):
        analysis_IO.write_acceptance_statistics(1, acceptance)


def test_write_performance_and_quench_file(cfg_paths):
    analysis_IO.write_performance(30, "E", 12.345, "00:01:00", "00:59:00")
    analysis_IO.write_quench_file(30, 298.15, -42.0)

    perf = _read(cfg_paths["OUTNAME_PERFORMANCE"])
    assert perf.startswith("30\tE\t12.35")
    assert perf.endswith("\t00:59:00\n")

    assert _read(cfg_paths["QUENCHFILE_NAME"]) == "30\t298.15\t  -42.0000\n"
