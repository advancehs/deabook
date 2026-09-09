"""Command-line entry point for the ``deabook`` package.

This module backs the ``deabook`` console script declared in
``pyproject.toml`` (``[project.entry-points."console_scripts"]``).
"""

import sys

from deabook import __version__


def main(argv=None):
    """Print the installed ``deabook`` version."""
    print(f"deabook {__version__}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
