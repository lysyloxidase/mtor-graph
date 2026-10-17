"""The command line is part of the public surface: every subcommand must run."""

import pytest

from mtor_graph.cli import main

COMMANDS = [
    ["stats"],
    ["stats", "--json"],
    ["search", "sestrin"],
    ["node", "RHEB"],
    ["node", "RHEB", "--refs"],
    ["module", "circadian"],
    ["path", "LEUCINE", "RPS6", "-k", "2"],
    ["refs", "1", "2"],
    ["validate"],
    ["tree", "--root", "MTORC1", "--depth", "1"],
    ["tree", "--root", "MTORC1", "--depth", "1", "--markdown"],
    ["tree", "--root", "MTORC1", "--depth", "1", "--json"],
    ["tree", "--root", "MTORC1", "--depth", "1", "--ascii", "--refs", "--show-evidence"],
    ["--evidence", "established,emerging", "stats"],
    ["--modules", "circadian,pharmacology", "stats"],
]


@pytest.mark.parametrize("argv", COMMANDS, ids=[" ".join(c) for c in COMMANDS])
def test_subcommand_runs(argv, capsys):
    main(argv)
    assert capsys.readouterr().out.strip()


def test_global_filter_rejects_an_unknown_value(capsys):
    with pytest.raises(SystemExit):
        main(["--evidence", "nonsense", "stats"])


def test_export_draw_and_html(tmp_path, capsys):
    main(["export", str(tmp_path / "exports")])
    main(["draw", str(tmp_path / "graph.png"), "--module", "circadian"])
    main(["html", str(tmp_path / "graph.html"), "--module", "circadian"])
    assert (tmp_path / "exports" / "mtor_interactome.graphml").exists()
    assert (tmp_path / "graph.png").stat().st_size > 1000
    assert "vis-network" in (tmp_path / "graph.html").read_text(encoding="utf-8")
