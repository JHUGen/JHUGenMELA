# SimpleParticleCollection_t {#SimpleParticleCollection_t}

The type `SimpleParticleCollection_t` is a typedef described in TVar.hh. It is a vector of `SimpleParticle_t` objects.
These correspond to the particle lists that are submitted to Mela::setInputEvent.
There are a few constructors for this object, each of which depend on what suits your needs.

## Constructors {#simpleparticlecollection_constructor}

There are 4 different ways to construct the `SimpleParticleCollection_t` object. This is due to both the aforementioned pyROOT and C++ ROOT
compatibility issues alongside the usage of syntactic sugar to allow for easier columnar data entry.

## collection_initializer_from_column {#collection_initializer_from_column_doc}

This initializes a SimpleParticleCollection_t from a list of PDG ids, and either a set of \f$P_x\f$, \f$P_y\f$, \f$P_z\f$, \f$E\f$ *OR* \f$P_t\f$, \f$\eta\f$, \f$\phi\f$, mass depending on the boolean
flag set (just like in @ref particle_initializer "particle_initializer"), and is documented @ref collection_initializer_from_column "here".

Example:

~~~~~~~~~~~~~{.py}
import Mela
# Let's make a pair of opposing gluons and 4 leptons
gluon_Collection = Mela.SimpleParticleCollection_t([21, 21], [0, 0], [0, 0], [12, -12], [12, 12]) #px, py, pz, E

lepton_Collection = Mela.SimpleParticleCollection_t(
     [13, -13, 11, -11],
     [pt1, pt2, pt3, pt4],
     [eta1, eta2, eta3, eta4],
     [phi1, phi2, phi3, phi4],
     [mu_M, mu_M, e_M, e_M],
     ptEtaPhi=True
)
~~~~~~~~~~~~~

### collection_initializer {#collection_initializer_doc}

This initializes a SimpleParticleCollection_t from a list of SimpleParticle_t objects, and is documented @ref collection_initializer "here".
One benefit of operating this way is that you can mix and match particles initialized via {\f$P_x\f$, \f$P_y\f$, \f$P_z\f$, \f$E\f$}and those via {\f$P_t\f$, \f$\eta\f$, \f$\phi\f$, \f$m\f$}.

Example:

~~~~~~~~~~~~~{.py}
import Mela
# Let's make a pair of opposing gluons and 4 leptons
g1 = Mela.SimpleParticle_t(21, 0, 0, 12, 12)
g2 = Mela.SimpleParticle_t(21, 0, 0, -12, 12)
gluon_Collection = Mela.SimpleParticleCollection_t([g1, g2])

l1 = Mela.SimpleParticle_t(13, px1, py1, pz1, E1)
l2 = Mela.SimpleParticle_t(-13, px2, py2, pz2, E2)
l3 = Mela.SimpleParticle_t(11, pt3, eta3, phi3, e_M) #mixing bases!
l4 = Mela.SimpleParticle_t(-11, pt4, eta4, phi4, e_M)
lepton_Collection = Mela.SimpleParticleCollection_t([l1, l2, l3, l4])
~~~~~~~~~~~~~

### Default Constructor {#simpleparticlecollection_empty}

This is an empty initializer that simply creates a "list-like" object you can add particles to.
You can use the @ref add_particle "Mela.SimpleParticleCollection_t.add_particle" method to do this. Just like [`collection_intializer`](collection_initializer_doc), one can mix and match bases of vectors between {\f$P_x\f$, \f$P_y\f$, \f$P_z\f$, \f$E\f$}and {\f$P_t\f$, \f$\eta\f$, \f$\phi\f$, \f$m\f$}.

Example:

~~~~~~~~~~~~~{.py}
import Mela
# Let's make a pair of opposing gluons and 4 leptons
gluon_Collection = Mela.SimpleParticleCollection_t()
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, 12, 12))
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, -12, 12))

lepton_Collection = Mela.SimpleParticleCollection_t([l1, l2, l3, l4])
lepton_Collection.add_particle(Mela.SimpleParticle_t(13, px1, py1, pz1, E1))
lepton_Collection.add_particle(Mela.SimpleParticle_t(-13, px2, py2, pz2, E2))
lepton_Collection.add_particle(Mela.SimpleParticle_t(11, pt3, eta3, phi3, e_M))
lepton_Collection.add_particle(Mela.SimpleParticle_t(-11, pt4, eta4, phi4, e_M))
~~~~~~~~~~~~~

## Methods {#simpleparticlecollection_methods}

The SimpleParticleCollection_t object has a few methods that you can use alongside iteration capabilities.

### add_particle {#add_particle}

Simply call `Mela.SimpleParticleCollection_t.add_particle(Mela.SimpleParticle_t P)` to add a particle to the collection.

Example in the @ref simpleparticlecollection_empty "Default Constructor" section.

### toList {#toList}

Simple call `Mela.SimpleParticleCollection_t.toList()` to return a normal Python list of SimpleParticle objects.

Example:

~~~~~~~~~~~~~{.py}
import Mela
# Let's make a pair of opposing gluons
gluon_Collection = Mela.SimpleParticleCollection_t()
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, 12, 12))
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, -12, 12))

print(gluon_Collection) #prints an ugly repr function
print(gluon_Collection.toList) #prints a list of SimpleParticles
~~~~~~~~~~~~~

### Iteration {#simpleparticlecollection_iter}

**Even without converting a SimpleParticleCollection_t into a list, one can still iterate over it!** Iteration over the `Mela.SimpleParticle_t` objects stored within can be accessed in a for loop. See the example below:

~~~~~~~~~~~~~{.py}
import Mela
# Let's make a pair of opposing gluons
gluon_Collection = Mela.SimpleParticleCollection_t()
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, 12, 12))
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, -12, 12))

for n, particle in enumerate(gluon_Collection):
     print(f"particle {n} has an id of {particle.id}")
     print(f"with P_x, P_y, P_z, E = {particle.vector}")
~~~~~~~~~~~~~

### Indexing {#simpleparticlecollection_indexing}

One can also index the underlying vector for the
`simpleParticleCollection_t` class. Indexing works as normal.

~~~~~~~~~~~~~{.py}
import Mela
gluon_Collection = Mela.SimpleParticleCollection_t()
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, 12, 12))
gluon_Collection[0] = Mela.SimpleParticle_t(21, 0, 0, -12, 12)
print(gluon_Collection[0].PxPyPzE_vector)
~~~~~~~~~~~~~

### Sum {#simpleparticelcollection_sum}

This method returns the 4-vector that is the sum
of all the 4-vectors stored inside the class. It returns
them as a tuple of 4 values: 
\f$P_x\f$, \f$P_y\f$, \f$P_z\f$, \f$E\f$.

~~~~~~~~~~~~~{.py}
import Mela
gluon_Collection = Mela.SimpleParticleCollection_t()
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, 12, 12))
gluon_Collection.add_particle(Mela.SimpleParticle_t(21, 0, 0, -12, 12))
print(gluon_Collection.Sum()) #should print (0,0,0,0)
~~~~~~~~~~~~~

### MTotal {#simpleparticlecollection_mtotal} 

This function is just like `SimpleParticleCollection_t.Sum`, however
it returns the total mass of the resultant 4-vector.

## Pickling Support

Should you desire to pickle a `SimpleParticleCollection_t` object,
this is supported through the Python MELA bindings. This support
for pickling also means that `SimpleParticleCollection_t`
supports the python `multiprocessing` package.
