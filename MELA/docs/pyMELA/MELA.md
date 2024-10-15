# MELA {#MELA_python}

This is what contains the bulk of what one desires to do with MELA within Python. The baseplate of this, is of course, the Mela class in `C++`.

## Constructors {#MELA_constructor}

The constructor inherits from the constructor in the @ref MELA_constructor_cpp "C++ implementation". 
All the default arguments work,
so one could create the Mela object with just a call to `Mela()`.

## Couplings {#MELA_couplings}

Couplings in MELA are what define the interactions between
different particles.

These are couplings using the amplitude basis for the Higgs Boson.
These are described in the following papers on the
main page for JHUGen linked [here](https://spin.pha.jhu.edu/).

In MadMELA (or the Madgraph matrix elements used within MELA),
the values set will colloquially be called "couplings" but are
actually the Wilson Coefficients that are used in Standard
Model Effective Field Theory (SMEFT). Under the
SMEFTSim framework, all the Wilson Coefficients being set to
0 corresponds to the Standard Model. The SMEFTSim
formulation is described in the following paper
[here](https://arxiv.org/abs/2012.11343).

### Standard Model Definition

In MELA, the Standard Model is defined as:

* `ghz1=1, ghg2=1` for JHUGen
* `ghz1=1, kappa_top=1, kappa_bot=1` for MCFM
* No couplings set for MadMELA (SMEFTSim default is Standard Model)

### Coupling Values

MELA in Python relies upon named couplings, rather than couplings
done through indices as in `C++`. See the following
table of all couplings @subpage MELA_couplings_table "here".
All of the possible couplings/Wilson Coefficients
you would want to set are defined by macros in mela_binding.cpp.
The definition for those macros can be found @ref py_macros "here".

**Most** couplings can be both real and imaginary. As a consequence,
**most** couplings are assigned via a Python iterable of size 2.

Some couplings are different. For example, all the
SMEFTSim couplings for usage within MadMELA are
single-valued real inputs.

~~~~~~~~~~~~~{.py}
import Mela
m = Mela.Mela()
m.ghz1 = [1,0] #sets ghz1 to 1
m.ghz4 = [0,1] #sets ghz4 to i
m.mdl_chwb = 1 #MadMELA coupling set to a real value!
~~~~~~~~~~~~~

### Coupling Arrays

One can also access the arrays for each of the coupling arrays
in the table above. This is sometimes useful
when you want to see the entire array at once.

Due to the way the Python bindings were made,
the calls to arrays need to be formulated as
function calls. This would look like the following:

~~~~~~~~~~~~~{.py}
import Mela
m = Mela.Mela()
print(
    m.selfDHggcoupl()
)

# You can also set other variables to the list to manipulate.
# Should you change the values of the list assigned
# to the variable, the actual values
# in MELA will change as well
selfDHggcoupl = m.selfDHggcoupl()
~~~~~~~~~~~~~

## Angles {#MELA_angles}

The MELA package can also be used to compute the angles of
a scattering process. Simply input a given MELA event
using Mela::SetInputEvent and run one of the following
three functions:

* Mela::computeDecayAngles
* Mela::computeVBFAngles
* Mela::computeVHAngles
* Mela::computeTTHAngles

Here is an example:

~~~~~~~~~~~~~{.py}
import Mela
m = Mela.Mela()
#Imagine making three SimpleParticle_t objects here
#to represent the mothers, daughters, and associated particles
m.SetInputEvent(daughters, associated, mothers)
m.setProcess(PROCESS, Mela.MatrixElement.JHUGen, PRODUCTION)
m4L, mZ1, mZ2, costheta1, costheta2, phi, costhetastar, phi1 = m.computeDecayAngles()
~~~~~~~~~~~~~
