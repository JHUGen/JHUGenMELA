# PyMela       {#PyMela_page}

There exists a Python implementation of MELA functions created via PyBind11. These functions are written in mela_binding.cpp.
There are three groups of functions that are implemented:

- Functions that are implemented as pass-by-reference in C++ that need to be converted to functions that return a value
- Functions that don't translate well between Python and C++ that need to be rewritten in a more compatible format
- Functions that just need to be referenced within the PyBind framework without any extra overhead

The first two categories of functions are detailed in the @ref Pychanges section. These changes vary, but generally exist for
compatibility or ease-of-use cases.

The following classes in Mela are encompassed by the PyBind11 module, and are documented in the following subpages:

- @subpage SimpleParticle_t "Mela.SimpleParticle_t"
- @subpage SimpleParticleCollection_t "Mela.SimpleParticleCollection_t"
- @subpage py_enums "Enumerations"
