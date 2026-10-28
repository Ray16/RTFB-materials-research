"""Compatibility entry point for installations made before the clean refactor."""
from redox.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
