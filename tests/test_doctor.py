import pytest

from healthcare import cli
from healthcare.doctor import collect_report


def test_report_contains_environment_facts():
    report = collect_report()
    assert report["python"].startswith("3.12")
    assert report["disk_free_gib"] > 0
    assert report["torch"]
    if report["cuda_available"]:
        assert report["gpu"]
        assert report["gpu_memory_gib"] > 0


def test_cli_requires_a_command():
    with pytest.raises(SystemExit):
        cli.main([])
