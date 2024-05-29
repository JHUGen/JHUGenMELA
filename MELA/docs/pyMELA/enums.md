# Stored Enumerations {#py_enums}

There are a variety of enumerated values stored in PyMELA in order to facilitate running MELA successfully. All enumerated values are stored as PyBind11 enumerations, which are documented [here](https://pybind11.readthedocs.io/en/stable/classes.html#enumerations-and-internal-types).

With every enumerated type one can look at all the members and their values as a dictionary via the `enum.__members__` function for a given enumeration. See the example below for more information.

~~~~~~~~~~~~~{.py}
import Mela
print(Mela.MatrixElement.__members__)
~~~~~~~~~~~~~

## Verbosity Level {#verbosity_enum}

`Mela.VerbosityLevel` controls how verbose MELA can be. 
Every verbosity level is a subset of the higher one (i.e. `Mela.VerbosityLevel.ERROR` contains a subset of the output from `Mela.VerbosityLevel.INFO`).

Its values are:
| Name | Value | Summary |
| ---- | ----- | ------- |
| `Mela.VerbosityLevel.SILENT` | 0 | Only *required* information |
| `Mela.VerbosityLevel.ERROR` | 1 | Only outputs *unexpected* behavior |
| `Mela.VerbosityLevel.INFO` | 2 | Outputs out useful information as well |
| `Mela.VerbosityLevel.DEBUG` | 3 | Outputs some barebones debugging information |
| `Mela.VerbosityLevel.DEBUG_VERBOSE` | 4 | Outputs more debugging information |
| `Mela.VerbosityLevel.DEBUG_MECHECK` | 5 | Outputs information directly relating to the matrix element |

## Matrix Element {#decay_enum}
