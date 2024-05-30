# Stored Enumerations {#py_enums}

There are a variety of enumerated values stored in PyMELA in order to facilitate running MELA successfully. All enumerated values are stored as PyBind11 enumerations, which are documented [here](https://pybind11.readthedocs.io/en/stable/classes.html#enumerations-and-internal-types).

With every enumerated type one can look at all the members and their values as a dictionary via the `enum.__members__` function for a given enumeration. See the example below for more information.

~~~~~~~~~~~~~{.py}
import Mela
print(Mela.MatrixElement.__members__)
~~~~~~~~~~~~~

There are 2 types of enumerations - those stored in TVar.hh, and those stored in TCouplingsBase.hh.
Generally, the former deals with values used to set up MELA, and the latter deals with defining where couplings exist.
These are tabulated in the following subpages:

- @subpage tvar_enums "TVar Enumerations"
- @subpage tcoupling_enums "TCouplingsBase Enumerations"
