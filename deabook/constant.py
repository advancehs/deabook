# Composite error term
CET_ADDI = "addi"
"""
CET_ADDI: Additive composite error term.
"""

CET_MULT = "mult"
"""
CET_MULT: Multiplicative composite error term.
"""

CET_Categories = {
    CET_ADDI: "Additive composite error term",
    CET_MULT: "Multiplicative composite error term"
}

CET_Model_Categories = {
    CET_ADDI: "additive model",
    CET_MULT: "multiplicative model"
}

# Frontier
FUN_PROD = "prod"
"""
FUN_PROD: Production frontier.
"""

FUN_COST = "cost"
"""
FUN_COST: Cost frontier.
"""

FUN_Categories = {
    FUN_PROD: "Production frontier",
    FUN_COST: "Cost frontier"
}

# Return to scale
RTS_VRS1 = "vrs1"
"""
RTS_VRS1: Variable returns to scale with same reduction factor.
"""
RTS_VRS2 = "vrs2"
"""
RTS_VRS2: Variable returns to scale with different reduction factor.
"""
RTS_CRS = "crs"
"""
RTS_CRS: Constant returns to scale.
"""

# Generic VRS (single reduction factor) used by the classical DEA/DDF
# interfaces.  ``RTS_VRS1``/``RTS_VRS2`` remain the explicit weak-disposability
# variants (same vs. different reduction factor).
RTS_VRS = "vrs"

RTS_Categories = {
    RTS_VRS: "Variable returns to scale",
    RTS_VRS1: "Variable returns to scale with same reduction factor",
    RTS_VRS2: "Variable returns to scale with different reduction factor",
    RTS_CRS: "Constant returns to scale"
}

RTS_ALLOWED = frozenset(RTS_Categories)

# Orientation
ORIENT_IO = "io"
"""
ORIENT_IO: Input orientation.
"""

ORIENT_OO = "oo"
"""
ORIENT_OO: Output orientation.
"""

ORIENT_BO = "bo"
"""
ORIENT_BO: Undesirable Output orientation.
"""
ORIENT_Categories = {
    ORIENT_IO: "Input orientation",
    ORIENT_OO: "Output orientation",
    ORIENT_BO: "Undesirable Output orientation"
}

# left-hand or right-hand Derivative
LEFT = "left-hand Derivative"
"""
LEFT: left-hand Derivative.
"""

RIGHT = "right-hand Derivative"
"""
RIGHT: right-hand Derivative.
"""

Derivative_Categories = {
    LEFT: "left-hand Derivative",
    RIGHT: "right-hand Derivative"
}


# Residual decomposition
RED_MOM = "MOM"
"""
RED_MOM: Method of moments.
"""

RED_QLE = "QLE"
"""
RED_QLE: Quassi-likelihood estimation.
"""

RED_KDE = "KDE"
"""
RED_KDE: Kernel deconvolution estimation.
"""

RED_Categories = {
    RED_MOM: "Method of moments",
    RED_QLE: "Quassi-likelihood estimation",
    RED_KDE: "Kernel deconvolution estimation"
}

# Optimization
OPT_LOCAL = "local"
OPT_DEFAULT = None



# technology
GLOBAL = "Global production technology"
"""
    window(#)                   use window production technology with the #-period bandwidth
    biennial                    use biennial production technology
    sequential                  use sequential production technology
    GLOBAL:  Global production technology.
"""
CONTEMPORARY = "Contemporary production technology"
"""
CONTEMPORARY:  Contemporary production technology.
"""
SEQUENTIAL = "sequential"
"""
SEQUENTIAL:  Sequential (cumulative) production technology -- the reference set for
period t pools every observation with year <= t.
"""
WINDOW = "window"
"""
WINDOW:  Window production technology -- the reference set for period t is a
centered local window of width ``window`` around t (e.g. [t-1, t, t+1]).
"""

TECH_Categories = {
    GLOBAL: "Global production technology",
    CONTEMPORARY:"Contemporary production technology",
    SEQUENTIAL: "Sequential production technology",
    WINDOW: "Window production technology",
}


# dynamic productivity index
MAL = " malquist prodcutivity index or malquist-luenberger prodcutivity index"

"""
    MAL              malquist prodcutivity index or malquist-luenberger prodcutivity index
"""
LUE = "luenberger prodcutivity index"
"""
    LUE              luenberger prodcutivity index
"""
DYNAMIC_Categories = {
    MAL : " malquist prodcutivity index or malquist-luenberger prodcutivity index",
    LUE : "luenberger prodcutivity index"
}


def validate_category(value, categories, name):
    """Validate a constant against a category dictionary and return it."""
    if value not in categories:
        raise ValueError(f"{name} must be one of {list(categories)}; got {value!r}.")
    return value


def validate_rts(value):
    """Validate returns-to-scale constants used throughout deabook."""
    return validate_category(value, RTS_Categories, "rts")


def validate_technology(value):
    """Validate dynamic technology constants."""
    return validate_category(value, TECH_Categories, "tech")


def validate_dynamic(value):
    """Validate dynamic-productivity index constants."""
    return validate_category(value, DYNAMIC_Categories, "dynamic")
