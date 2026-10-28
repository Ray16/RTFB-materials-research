"""Module entry point, identical to the installed redox command."""
from redox.cli import COMMANDS, main, runpy

if __name__ == "__main__":
    raise SystemExit(main())
