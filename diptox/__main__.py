"""Run the non-interactive DiPTox CLI with ``python -m diptox``."""


if __name__ == "__main__":
    from multiprocessing import freeze_support

    freeze_support()
    from .command_line import main

    raise SystemExit(main())
