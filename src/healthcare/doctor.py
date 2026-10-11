import platform
import shutil
import sys

import torch

from healthcare.paths import PROJECT_ROOT

GIB = 1024**3


def collect_report() -> dict:
    cuda_available = torch.cuda.is_available()
    report = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "torch": torch.__version__,
        "cuda_available": cuda_available,
        "cuda_version": torch.version.cuda,
        "gpu": None,
        "gpu_memory_gib": None,
        "disk_free_gib": round(shutil.disk_usage(PROJECT_ROOT).free / GIB, 1),
    }
    if cuda_available:
        props = torch.cuda.get_device_properties(0)
        report["gpu"] = props.name
        report["gpu_memory_gib"] = round(props.total_memory / GIB, 1)
    return report


def run() -> int:
    report = collect_report()
    width = max(map(len, report))
    for key, value in report.items():
        print(f"{key.ljust(width)}  {value}")
    if not report["cuda_available"]:
        print("CUDA: unavailable", file=sys.stderr)
        return 1
    print("CUDA: available")
    return 0
