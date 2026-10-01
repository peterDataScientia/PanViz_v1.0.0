from pathlib import Path


def get_version() -> str:
    version_file = Path(__file__).with_name("VERSION.txt")
    try:
        value = version_file.read_text(encoding="utf-8").strip()
    except OSError:
        return "0.0.0"
    return value or "0.0.0"


PANVIZ_VERSION = get_version()
