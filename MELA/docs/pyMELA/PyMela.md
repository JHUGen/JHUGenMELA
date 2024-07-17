# PyMela       {#PyMela_page}

There exists a Python implementation of MELA functions created via PyBind11. These functions are written in mela_binding.cpp.
There are three groups of functions that are implemented:

- Functions that are implemented as pass-by-reference in C++ that need to be converted to functions that return a value
- Functions that don't translate well between Python and C++ that need to be rewritten in a more compatible format
- Functions that just need to be referenced within the PyBind framework without any extra overhead

The first two categories of functions are detailed in the @ref Pychanges section. These changes vary, but generally exist for
compatibility or ease-of-use cases.

## Other Classes Encompassed

The following classes in the C++ implementation are encompassed by the PyBind11 module, and are documented in the following subpages:

- @subpage SimpleParticle_t "Mela.SimpleParticle_t"
- @subpage SimpleParticleCollection_t "Mela.SimpleParticleCollection_t"
- @subpage event_scales_type "Mela.event_scales_type"
- @subpage MELAParticle "Mela.MELAParticle"
- @subpage MELAThreeBodyDecayCandidate "Mela.MELAThreeBodyDecayCandidate"
- @subpage MELACandidate "Mela.MELACandidate"
- @subpage py_enums "Enumerations"

These classes serve to help the functionality of MELA such as providing an avenue to input vectors, index values through enumerations, etc.

## MELA itself {#PyMela_doc}

MELA in Python is a (quite large!) class within the PyBind11 module stated at the top of the page. All the functions are either defined at @ref Pychanges "the top of the file", or is taken directly from the header, Mela.h. They are defined in the @subpage MELA_python "following subpage".
