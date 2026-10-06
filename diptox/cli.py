"""Command-line entry point for the local DiPTox web interface."""

import os
import sys
from pathlib import Path


if not __package__:
    # Direct file execution (including IDE launches) must import this checkout,
    # not a different DiPTox version installed in the active environment.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "diptox"


def run_gui() -> int:
    os.environ["DIPTOX_GUI_MODE"] = "true"
    from .logger import log_manager

    log_manager.configure(enable_file=False)
    try:
        from .gui import main
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("nicegui"):
            raise RuntimeError(
                "The DiPTox GUI requires NiceGUI. Install the GUI dependencies "
                "with `pip install nicegui` or reinstall DiPTox."
            ) from exc
        raise
    return main()


if __name__ == "__main__":
    raise SystemExit(run_gui())
