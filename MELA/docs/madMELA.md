# madMELA {#madMELA}

Our approach to matrix elements in Madgraph is different to that of EFT2Obs.
In keeping with the way that MELA is structured, there are a select few
processes/production modes that are supported through the native MELA ecosystem.
All the final states are 4-leptons. This is colloquially referenced
as "MadMELA" throughout the rest of the text, and is managed by MadMela.hh.

So far, these production modes (and their respective processes) are:

| Production | Process |
| ---------- | ------- |
| quark-antiquark production | Signal  |
| quark-antiquark production | Background  |
| quark-antiquark production | BSI  |
| gluon-fusion | Signal |

Where

* "Signal" indicates the presence of a Higgs Boson,
* "Background" a process with the same initial and final state
but without the presence of a Higgs
* "BSI" the sum of signal, background, and interference between the two

## Configuration {#madMELA_configure}

MadMELA is configured through a series of Wilson Coefficients
and functions. These coefficients are set in the same manner as
the couplings used in JHUGen or MCFM, but are all
real-valued.

### C++

In `C++`, couplings are set using the `SelfDSMEFTSimCoupl`
array (Mela::SelfDSMEFTSimCoupl), with indices corresponding to the indices shown @ref table_SMEFT "here".

### Python

In Python, couplings *can* be set using `SelfDSMEFTSimCoupl`, similar
to `C++`. However, Python utilizes named couplings to self-configure. These can
also be seem in the same table @ref table_SMEFT "referenced above".

The documentation for setting various items in Python can be seen in @ref PyMela_page "here"

## Compilation {#madMELA_compile}

MadMELA is compiled through combining multiple different
Madgraph stand-alone reweighting areas into one
shareable library that is used by the C++ Mela
interface through MadMela.hh.

The steps are as follows:

1. Each different sub-area is compiled to form a static library for 
that process.
2. While each area is being compiled, a large common block that contains
all the possible parameters for each area is created. This large
common block contains all the numbered general couplings
that are possible for each Madgraph reweighting area.
3. This large common block replaces the respective
common blocks for each area through symbolic links
pointing to the large block placed in a separate dedicated
directory.
4. All the compute functions and other commonly-named items
for each given area are renamed in the symbol table using `objcopy`
to something else
unique and appropriate.
5. All of these are collated together and compiled
