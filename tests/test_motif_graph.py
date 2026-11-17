"""Contracts for the bag-of-tasks campaign DAG that replaced the serial motifs loop.

Pure scheduling/graph-construction tests: no GPU, DFT, UMA or any QM is launched.
"""
from pathlib import Path

import pytest

from redox.workflows import motifs
from redox.workflows.stages import validate_graph

CONFIG = Path("config/custom_candidates.py")
GROUPS = [{"id": "q1", "states": [("neu", 0, 1, 2), ("ox1", 1, 2, 1), ("ox2", 2, 1, 1)]},
          {"id": "q2", "states": [("neu", 0, 1, 2), ("ox1", 1, 2, 1)]}]


def _graph(groups=GROUPS, thermal_definition="native"):
    return motifs.build_graph(CONFIG, groups, Path("results/run"), "gpu", 2, 3,
                              thermal_definition, ferrocene_states=("neu", "ox"))


def _by_name(nodes):
    return {n.name: n for n in nodes}


def test_graph_is_acyclic_and_every_dependency_exists():
    nodes = _graph()
    validate_graph(nodes)  # raises on a cycle or a missing dependency
    names = {n.name for n in nodes}
    for node in nodes:
        assert set(node.dependencies) <= names


def test_intra_candidate_edges_pipeline_per_state_dft():
    nodes = _by_name(_graph())
    # uma -> dft(per state) -> sp : DFT is no longer one barrier node.
    assert nodes["uma:q1"].dependencies == ("build:q1",)
    dft_states = [f"dft:q1:{s}" for s in ("neu", "ox1", "ox2")]
    for name in dft_states:
        assert nodes[name].dependencies == ("uma:q1",)
    assert set(nodes["sp:q1"].dependencies) == set(dft_states)
    # sp -> coverage -> reorg (coverage gates reorg); reorg also needs sp.
    assert "sp:q1" in nodes["coverage:q1"].dependencies
    assert nodes["cross_points:q1"].dependencies == ("coverage:q1", "sp:q1")
    # native thermal follows sp; coverage (readiness --strict) needs the record.
    assert "sp:q1" in nodes["native_thermal:q1"].dependencies
    assert "native_thermal:q1" in nodes["coverage:q1"].dependencies


def test_ferrocene_is_computed_once_and_fc_sync_follows_every_native_thermal():
    nodes = _graph(groups=[{"id": f"c{i}", "states": [("neu", 0, 1, 2), ("ox1", 1, 2, 1)]}
                           for i in range(10)])
    fc_prereq = [n for n in nodes if n.name == "native_thermal:ferrocene"]
    assert len(fc_prereq) == 1
    assert fc_prereq[0].lock == "ferrocene" and fc_prereq[0].dependencies == ()
    # No candidate's native_thermal re-touches the shared ferrocene namespace.
    for n in nodes:
        if n.label == "native_thermal" and n.gid != "ferrocene":
            assert n.argv[n.argv.index("--only") + 1] == n.gid
            assert "native_thermal:ferrocene" in n.dependencies
    fc_sync = _by_name(nodes)["fc_sync"]
    native_thermal_names = {n.name for n in nodes if n.label == "native_thermal"}
    assert set(fc_sync.dependencies) == native_thermal_names


def test_candidate_argv_is_byte_identical_to_the_serial_command_list():
    # Scheduling-only refactor: every per-stage argv must match commands() exactly,
    # except native_thermal's --only target, which is split off ferrocene.
    serial = motifs.commands(CONFIG, ["q1"], Path("results/run"), "gpu", 2, 3, "native")
    by_label = {}
    for label, _, argv in serial:
        by_label.setdefault(label, argv)
    nodes = _by_name(_graph())
    assert list(nodes["build:q1"].argv) == by_label["build"]
    assert list(nodes["uma:q1"].argv) == by_label["uma"]
    assert list(nodes["sp:q1"].argv) == by_label["sp"]
    assert list(nodes["coverage:q1"].argv) == by_label["coverage"]
    assert list(nodes["cross_points:q1"].argv) == by_label["cross_points"]
    # dft per-state argv is the serial dft argv with the state appended to --only.
    serial_dft = by_label["dft"]
    at = serial_dft.index("--only") + 1
    assert list(nodes["dft:q1:neu"].argv) == [*serial_dft[:at], "q1:neu", *serial_dft[at + 1:]]


def test_legacy_mode_has_no_ferrocene_or_fc_sync_nodes():
    labels = {n.label for n in _graph(thermal_definition="legacy")}
    names = {n.name for n in _graph(thermal_definition="legacy")}
    assert "native_thermal" not in labels and "fc_reference_sync" not in labels
    assert "native_thermal:ferrocene" not in names and "fc_sync" not in names


def test_resume_predicate_skips_only_complete_valid_units(tmp_path, monkeypatch):
    from redox.core import common
    from redox.backends import dft as dft_backend
    dft_root = tmp_path / "dft"
    uma_root = tmp_path / "uma"
    monkeypatch.setattr(common, "DFT", dft_root)
    monkeypatch.setattr(common, "UMA", uma_root)

    # dft node: needs geometry + result and a current optimization protocol.
    dft_node = motifs.MotifNode("dft:q1:neu", "dft", (), True, gid="q1", state="neu")
    assert motifs.node_done(dft_node) is False  # nothing on disk
    sd = dft_root / "q1" / "neu"; sd.mkdir(parents=True)
    (sd / "result.json").write_text('{"charge": 0}')
    (sd / "opt.xyz").write_text("x")
    monkeypatch.setattr(dft_backend, "_opt_protocol_mismatch", lambda r, q: "stale")
    assert motifs.node_done(dft_node) is False  # stale geometry -> recompute
    monkeypatch.setattr(dft_backend, "_opt_protocol_mismatch", lambda r, q: "")
    assert motifs.node_done(dft_node) is True

    # uma node: every declared state needs a result on disk.
    uma_node = motifs.MotifNode("uma:q1", "uma", (), True, gid="q1", states=("neu", "ox1"))
    assert motifs.node_done(uma_node) is False
    for s in ("neu", "ox1"):
        d = uma_root / "q1" / s; d.mkdir(parents=True)
        (d / "result.json").write_text("{}")
    assert motifs.node_done(uma_node) is True

    # sp node: every state needs an active-protocol single point.
    sp_node = motifs.MotifNode("sp:q1", "sp", (), True, gid="q1", states=("neu", "ox1"))
    monkeypatch.setattr(common, "read_result", lambda gid, s: {"sp_status": "ok"})
    assert motifs.node_done(sp_node) is True
    monkeypatch.setattr(common, "read_result", lambda gid, s: {"sp_status": "missing"})
    assert motifs.node_done(sp_node) is False

    # native thermal node: a legacy RRHO record must NOT count as done.
    nt_node = motifs.MotifNode("native_thermal:q1", "native_thermal", (), True,
                               gid="q1", states=("neu", "ox1"))
    monkeypatch.setattr(common, "read_result",
                        lambda gid, s: {"thermal_definition": "legacy_gfn2_gas_rrho", "thermal_qc": "ok"})
    assert motifs.node_done(nt_node) is False
    monkeypatch.setattr(common, "read_result",
                        lambda gid, s: {"thermal_definition": "native_solution_dft_qrrho_spin",
                                        "thermal_qc": "ok"})
    assert motifs.node_done(nt_node) is True

    # cheap CPU / orchestration nodes are always re-run (idempotent, no GPU reserved).
    assert motifs.node_done(motifs.MotifNode("build:q1", "build", (), False, gid="q1")) is False
    assert motifs.node_done(motifs.MotifNode("fc_sync", "fc_reference_sync", (), False)) is False
