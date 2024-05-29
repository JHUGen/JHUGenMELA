# SimpleParticle_t {#SimpleParticle_t}

The SimpleParticle_t for Python is a type definition defined in TVar.hh. These are defined as pairs of an integer and a 4-vector, and are the
primary input object used for every particle.

Due to difficulties with interfacing between ROOT's pyROOT implemenetation of
[TLorentzVector](https://root.cern.ch/doc/master/classTLorentzVector.html), there are new functions that operate as constructors.

## Constructor {#simpleparticle_constructor}

There is one constructor defined in the @ref particle_initializer "particle_initializer" function. This was required due to the incompatibility of pyROOT and C++ based ROOT. The function takes in 5 required and 1 optional variable that are named as
`x`, `y`, `z`, `e`, and `ptEtaPhi`.

Example:

~~~~~~~~~~~~~{.py}
import Mela
gluon = Mela.SimpleParticle_t(id=21, x=0, y=0, z=12, e=12)
lepton = Mela.SimpleParticle_t(id=25, x=1, pt=0, eta=0, m=125, ptEtaPhi=True)
~~~~~~~~~~~~~

## Attributes {#simpleparticle_getter}

There are 2 read-only attributes for a SimpleParticle_t. These cannot be edited; make a new object if you want to make a new vector.

### id {#simpleparticle_id}

The id returns the PDG ID of the particle.

Example:

~~~~~~~~~~~~~{.py}
import Mela
gluon = Mela.SimpleParticle_t(id=21, x=0, y=0, z=12, e=12)
print(gluon.id) #will print 21
~~~~~~~~~~~~~

### vector {#simpleparticle_vector}

The vector attribute returns a 4-vector as a tuple of 4 values, the \f$P_x\f$, \f$P_y\f$, \f$P_z\f$, \f$E\f$ of the particle.
This representation of the vector will always be printed, even if the SimpleParticle_t was instantiated using \f$P_t\f$, \f$\eta\f$, \f$\phi\f$ and mass.

Example:

~~~~~~~~~~~~~{.py}
import Mela
gluon = Mela.SimpleParticle_t(id=21, x=0, y=0, z=12, e=12)
print(gluon.vector) #will print (0, 0, 12, 12)
~~~~~~~~~~~~~
