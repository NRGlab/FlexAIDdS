import csv
from pathlib import Path

from flexaidds.results import load_results


def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {line}\n" for line in remarks]
    lines.extend(
        [
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n",
            "END\n",
        ]
    )
    path.write_text("".join(lines), encoding="utf-8")


def test_load_results_groups_poses_by_mode(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_1_pose_2.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 2",
            "CF = -39.0",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            "free_energy = -34.2",
            "temperature = 300.0",
        ],
    )

    result = load_results(tmp_path)

    assert result.n_modes == 2
    assert result.temperature == 300.0
    assert result.binding_modes[0].mode_id == 1
    assert result.binding_modes[0].n_poses == 2
    assert result.binding_modes[0].best_cf == -42.5
    assert result.binding_modes[0].free_energy == -41.0
    assert result.binding_modes[1].mode_id == 2
    assert result.binding_modes[1].best_cf == -35.0


def test_load_results_uses_filename_heuristics_when_remarks_missing(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "mode_7_pose_3.pdb",
        [
            "CF = -11.5",
            "temperature = 298.15",
        ],
    )

    result = load_results(tmp_path)
    mode = result.binding_modes[0]
    pose = mode.poses[0]

    assert mode.mode_id == 7
    assert pose.pose_rank == 3
    assert mode.best_cf == -11.5


def _make_two_mode_result(tmp_path: Path):
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        ["binding_mode = 1", "pose_rank = 1", "CF = -42.5", "free_energy = -41.0", "temperature = 300.0"],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        ["binding_mode = 2", "pose_rank = 1", "CF = -35.0", "free_energy = -34.2", "temperature = 300.0"],
    )
    return load_results(tmp_path)


def test_to_csv_returns_string_when_path_is_none(tmp_path: Path) -> None:
    result = _make_two_mode_result(tmp_path)
    csv_text = result.to_csv()
    assert csv_text is not None
    rows = list(csv.DictReader(csv_text.splitlines()))
    assert len(rows) == 2
    assert rows[0]["mode_id"] == "1"
    assert rows[1]["mode_id"] == "2"


def test_to_csv_writes_file(tmp_path: Path) -> None:
    result_dir = tmp_path / "results"
    result_dir.mkdir()
    result = _make_two_mode_result(result_dir)
    csv_path = tmp_path / "output.csv"
    ret = result.to_csv(csv_path)
    assert ret is None
    assert csv_path.exists()
    rows = list(csv.DictReader(csv_path.read_text(encoding="utf-8").splitlines()))
    assert len(rows) == 2
    assert float(rows[0]["free_energy"]) == -41.0


def test_to_csv_columns_match_to_records(tmp_path: Path) -> None:
    result = _make_two_mode_result(tmp_path)
    records = result.to_records()
    csv_text = result.to_csv()
    rows = list(csv.DictReader(csv_text.splitlines()))
    assert list(rows[0].keys()) == list(records[0].keys())
