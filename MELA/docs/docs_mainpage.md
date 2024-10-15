# Introduction {#mainpage}

MELA is an important tool that was used for the Higgs boson discovery and for precise measurements of its structure and interactions. Please see the [website](https://spin.pha.jhu.edu/) and papers cited there for more details, and kindly cite those papers when using this code.

@authors I. Anderson, S. Bolognesi, F. Caola, J. Davis, Y. Gao, A. V. Gritsan,
@authors L. S. Mandacaru Guerra, Z. Guo, L. Kang, S. Kyriacou, C. B. Martin, T. Martini,
@authors K. Melnikov, R. Pan, M. Panagiotou, R. Rontsch, J. Roskes, U. Sarica,
@authors M. Schulze, M. V. Srivastav, N. V. Tran, A. Whitbeck, M. Xiao, Y. Zhou

This is the documentation for JHUGenMELA, a package that produces Matrix Element Calculations as used in JHUGen, MCFM-JHUGen, and Madgraph.

The main user-interface in C++ can be found at the documentation for the Mela class,
 while the Python can be found in the @ref PyMela_page "PyMela" section.

**Documentation written by Mohit Srivastav in 2024.**

## Installation {#install_sec}

Mela is a part of the JHUGen package, and is installed by downloading a tar file located on the [website](https://spin.pha.jhu.edu/). Once you have un-tar'd the file, you will have to compile JHUGenMELA (alongside JHUGen, MCFM, and anything else you wish to use within the JHUGen package).

Compile using

```{.sh}
./setup.sh -j N # (N being the number of cores to compile, or blank for max. allowed)
```

* The script decides whether to use integration with the experiments' software.
* Either the one in this folder, or that in the main JHUGenMELA folder works; they do the same thing.
* Feedback on how this script works in different environments is appreciated.

Once you compile, make sure to follow the on-screen instructions to set up all extra environment variables, i.e., run

```{.sh}
eval $(./setup.sh env) # or equivalents in non-bash shells
```

## Previous Tutorials {#prev_tutorials}

There are a few tutorials lying around for MELA. They are listed below:

* [LPC EFT Workshop Notre Dame](https://indico.cern.ch/event/1378665/timetable/?view=standard#30-mela-tools)
* [LPC Offshell Workshop at Fermilab](https://indico.cern.ch/event/1375252/timetable/#8-mc-generators-2-mela-tools)

## Other Presentations involving MELA {#prev_presentations}

There have also been a few presentations regarding improvements to MELA,
and about the MELA package in general. Refer to these presentations
for more information.

* [CMS GEN Meeting September 2023](https://indico.cern.ch/event/1316531/#42-new-features-in-the-jhugen)

## Code Examples {#MELA_examples}

There are also examples that may prove useful in your own MELA endeavors.
The repository for the aforementioned MELA tutorials are located below:

* [LPC EFT Workshop Repo](https://github.com/MohitS704/EFT-Workshop-JHUGenMELA-tutorial/blob/main/)
* [LPC Off-Shell Higgs Workshop Repo **WARNING BUGS**](https://github.com/Offshell-Workshop-LPC/JHUGen-MELA-tutorial)

On a more general note - the steps to begin re-weighting events with MELA are as follows:

1. Set the process, production, and matrix element. These values are defined
for both @ref tvar_enums "Python" and in the `C++` in the TVar namespace.
Use Mela::SetProcess to do this step.
1. Set the input event using Mela::SetInputEvent.
2. Set coupling values for MELA (examples for Python are available
in documentation found @ref MELA_python "here")
1. Call a compute function of some kind. This is either:
    * Mela::computeP
    * Mela::computeProdP
    * Mela::computeProdDecP

## Known Interfaces {#MELA_interface}

There are also two MELA interfaces.
These interfaces are two different methods of configuring MELA
and utilizing it, usually starting from a ROOT file.

### MELACalc {#MELAcalc}

MELAcalc is a fully Pythonic interface to using ROOT files
reweighted with MELA. One does not need `CERN ROOT`, only
`uproot`, `numpy`, and `MELA` itself.
MELAcalc is currently on version 3 - you can find it
[here](https://github.com/hexutils/MELAcalc/blob/main/).
The configuration for the interface is done through `.json` files.
You can find the documentation for MELAcalc here: **TODO: MELACALC DOCS**.

### MELAAnalytics {#MELAAnalytics}

MELAAnalytics is a `C++` based interface to MELA that utilizes
csv files as configurations to MELA. MELAAnalytics has
been used for analyses in the past. You can find it
[here](https://github.com/MELALabs/MelaAnalytics/)

## Structure {#MELA_structure}

The root of MELA is underlying FORTRAN code. From there, there is a
`C++` interface that adds a variety of functionality such as
angle calculation, and the ability to interface to other
different matrix-element calculators like `MCFM` or `Madgraph`. The
functions in the `C++` implementation are documented throughout
this page, but especially in the source code documentation for Mela.

From there, there is a Python implementation that utilizes `PyBind11`
to generate the Python package. It is somewhat unique in nature,
and so is documented differently @subpage PyMela_page "here"

## Acknowledgements {#Other_Packages}

JHUGenMELA is a tool that also builds off of the work of others.
Their work is credited here.

### MCFM {#MCFM}

MCFM-JHUGen is built off of MCFM, and its printout is shown below:

```text
************** MCFM - version 7.0 ****************
*                                                *
* MCFM, v7.0                   March 20th, 2015  *
*                                                *
* Authors: John Campbell, Keith Ellis,           *
*          Walter Giele, Ciaran Williams         *
*         (johnmc@fnal.gov, ellis@fnal.gov,      *
*          giele@fnal.gov,ciaran@fnal.gov)       *
*                                                *
* For details see:                               *
*                                                *
*  arXiv:1502.02990 (VBF and VBS Higgs)          *
*  arXiv:1403.2641  (Triphoton production)       *
*  arXiv:1312.1628  (gg->WW, Higgs interference) *
*  arXiv:1311.3589  (gg->ZZ, Higgs interference) *
*  Phys.Rev.D87:114006, arXiv:1302.3856          *
*  (tZ, tH -- with R. Rontsch)                   *
*  arXiv:1211.6390 (DM, P. Fox and C. Williams)  *
*  JHEP 1211:162 (2012), arXiv:1208.0566         *
*  (Z+gam+jet,Z+gam+gam -- with H. Hartanto)     *
*  arXiv:1204.1513 (top production+decay)        *
*  JHEP 1207:052 (2012), arXiv:1204.5678 (ttW)   *
*  JHEP 1110:005 (2011), arXiv:1107.5569         *
*         (gg->WW,Higgs intference)              *
*  JHEP 1107:018 (2011), arXiv:1105.0020         *
*         (diboson update)                       *
*  JHEP 1103:027 (2011), arXiv:1011.6647         *
*         (Wbb for mb>0, with S. Badger)         *
*  Phys.Rev.D81:074023, arXiv:1001.4495 (H+2jet) *
*                                                *
*  P.R.L. 102:142001, arXiv:0903.0005 [hep-ph]   *
*    (t-channel single top + explicit b,         *
*      JC, R.Frederix, F.Maltoni, F.Tramontano)  *
*  N.P.B 726:109(2005), hep-ph/0506289 (W+t)     *
*  Phys.Rev.D70:094012, hep-ph/0408158 (Sngl Top)*
*       (with Francesco Tramontano)              *
*                                                *
*  Phys.Rev.D65:113007, hep-ph/0202176 (W,Z+2j)  *
*  Phys.Rev.D62:114012, hep-ph/0006304 (W,Z+bb)  *
*  Phys.Rev.D60:113006, hep-ph/9905386 (diboson) *
*                                                *
* On the web:  http://mcfm.fnal.gov/             *
*                                                *
*  With modifications from the JHUGen Team:      *
*  https://spin.pha.jhu.edu/                     *
*                                                *
*  Special thanks from JHUGen Team               *
*  to Oscar Eboli for anomalous Zqq couplings    *
*  in ggZZ process                               *
**************************************************
```

MCFM matrix elements are used in a variety of different calculations, and can be compiled in the JHUGen package by calling `./compile.sh` within the MCFM-JHUGen folder.

### Madgraph {#Madgraph}

There are certain processes that are use matrix elements generated with Madgraph utilizing the SMEFTsim package. MG5_aMC@nlo is used to generate these matrix elements, with all of them being at leading order. The printout for Madgraph with SMEFTsim is located below:

```text
**********************************************************
*                   MadGraph5_aMC@NLO                    *
**********************************************************
*                                                        *
*                      Going Beyond                      *
*                                                        *
*              http://madgraph.hep.uiuc.edu              *
*             http://madgraph.phys.ucl.ac.be             *
*                http://amcatnlo.cern.ch                 *
*                                                        *
*               The MadGraph5_aMC@NLO team               *
*                                                        *
*             Utilizing the SMEFTsim Package             *
*    I. Brivio, Y. Jiang, M. Trott, arXiv: 1709.06492    *
*              I. Brivio, arXiv: 2012.11343              *
*                                                        *
**********************************************************
```

Madgraph matrix elements are currently used to reweight processes
that were generated via Madgraph, colloquially
referred to as "MadMELA", (which
will be the preferred nomenclature moving 
forward). The Madgraph matrix elements
that we provide are generated using the `SMEFTsim` package.

Information about "MadMELA" can be found @subpage madMELA "here".