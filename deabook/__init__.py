"""Top-level package for deabook."""

__author__ = """advancehs"""
__email__ = "1019753743@qq.com"
__version__ = "0.0.11"

from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace


def data_path(name):
    """Return the absolute path to a bundled example data file.

    Resolves to the installed package location (works for both
    ``pip install deabook`` and editable installs), independent of the
    current working directory. Manifests call this via
    ``pd.read_excel(deabook.data_path("oecd data.xlsx"))``.
    """
    try:
        return str(files("deabook").joinpath("example", "book_cases", "data", name))
    except (ModuleNotFoundError, FileNotFoundError, TypeError):
        # Fallback for running directly from the source tree
        return str(Path(__file__).resolve().parent / "example" / "book_cases" / "data" / name)
# from . import CNLS
# from . import CNLSSDFDDF
# from . import CNLSSDFDDFweak
# from . import CNLSDDFG
# from . import CNLSG
from . import constant
# from . import CQER
# from . import CQERDDF
# from . import CQERG
# from . import CSVR
# from . import dataset
from . import DEA
from . import DEAweak
from . import DEAmb
from . import panel
from .panel import Panel
from .panel import StoNEDPanel
# Canonical panel namespace.  ``PANEL.tech`` is a ``Panel`` subclass used directly as
# the panel-solver constructor (``PANEL.tech(...)``), replacing the older
# ``PANEL.Panel(...)``.  The technology-type constants (``GLOBAL`` / ``CONTEMPORARY``
# / ``SEQUENTIAL`` / ``WINDOW``) are NOT exposed under ``PANEL.tech`` nor as bare names
# on ``deabook``; they live in ``deabook.constant`` and are imported as
# ``from deabook.constant import CONTEMPORARY`` (etc.), then passed via ``tech=``.


class _TechPanel(Panel):
    """``Panel`` subclass exposed as ``PANEL.tech`` so that ``PANEL.tech(...)``
    constructs a panel solver.  Technology-type constants live in
    ``deabook.constant`` (``GLOBAL`` / ``CONTEMPORARY`` / ``SEQUENTIAL`` /
    ``WINDOW``), not here."""


PANEL = SimpleNamespace(
    Panel=Panel,
    StoNEDPanel=StoNEDPanel,
    tech=_TechPanel,
)
# from . import FDH
# from . import plot
from . import StoNED
# from . import weakCNLS
# from . import weakCNLSb
# from . import weakCNLSbG
# from . import weakCNLSDDF
# from . import weakCNLSDDFG
# from . import weakCNLSG
# from . import weakCNLSx
# from . import weakCNLSxG
__all__ = [
    # 'CNLS',
    'CNLSSDFDDF',
    'CNLSSDFDDFweak',
    # 'CNLSDDFG',
    # 'CNLSG',
    'constant',
    # 'CQER',
    # 'CQERDDF',
    # 'CQERG',
    # 'CSVR',
    # 'dataset',
    'DEA',
    'DEAweak',  
    'DEAmb',
    'panel',
    'Panel',
    'StoNEDPanel',
    'PANEL',
    # 'FDH',
    # 'plot',
    'StoNED',
    # 'weakCNLS',
    # 'weakCNLSb',
    # 'weakCNLSbG',
    # 'weakCNLSDDF',
    # 'weakCNLSDDFG',
    # 'weakCNLSG',
    # 'weakCNLSx',
    # 'weakCNLSxG',
]
