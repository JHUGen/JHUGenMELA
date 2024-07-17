# TVar Enumerations {#tvar_enums}

These are the enumerated values set in TVar.hh. They are used to set up the settings in MELA.

## Verbosity Level {#verbosity_enum}

`Mela.VerbosityLevel` controls how verbose MELA can be. These are originally defined in TVar::VerbosityLevel.
Every verbosity level is a subset of the higher one (i.e. `Mela.VerbosityLevel.ERROR` contains a subset of the output from `Mela.VerbosityLevel.INFO`). The values are tabulated below:

| Name | Value | Summary |
| ---- | ----- | ------- |
| `Mela.VerbosityLevel.SILENT` | 0 | Only *required* information |
| `Mela.VerbosityLevel.ERROR` | 1 | Only outputs *unexpected* behavior |
| `Mela.VerbosityLevel.INFO` | 2 | Outputs out useful information as well |
| `Mela.VerbosityLevel.DEBUG` | 3 | Outputs some barebones debugging information |
| `Mela.VerbosityLevel.DEBUG_VERBOSE` | 4 | Outputs more debugging information |
| `Mela.VerbosityLevel.DEBUG_MECHECK` | 5 | Outputs information directly relating to the matrix element |

One can set the verbosity in the @ref MELA_constructor "MELA constructor" to begin with, or use `Mela.SetVerbosity(Mela.VerbosityLevel)`, which is originally defined in Mela::setVerbosity.

## Matrix Element {#matel_enum}

`Mela.MatrixElement` controls which matrix element MELA is using for its calculation. These are originally defined in TVar::MatrixElement.
The values are tabulated below:

| Name | Value | Summary |
| ---- | ----- | ------- |
| `Mela.MatrixElement.MCFM` | 0 | Evaluates via JHUGen-MCFM |
| `Mela.MatrixElement.JHUGen` | 1 | Evaluates via pure JHUGen |
| `Mela.MatrixElement.ANALYTICAL` | 2 | Evaluates via analytic methods |

One sets the matrix element using `Mela.setProcess(Mela.Process, Mela.MatrixElement, Mela.Production)`, which is originally defined in Mela::setProcess.

## Production {#prod_enum}

`Mela.Production` controls what production mode MELA is using for its calculations. These are originally defined in TVar::Production.
The values are tabulated below:

| Name | Value | Summary |
| ---- | ----- | ------- |
| `Mela.Production.ZZGG` | 0 | Gluon Fusion production |
| `Mela.Production.ZZQQB` | 1 | Quark-Antiquark production |
| `Mela.Production.ZZQQB_STU` | 2 | Quark-Antiquark production |
| `Mela.Production.ZZINDEPENDENT` | 3 | Production-Independent Calculation |
| `Mela.Production.ttH` | 4 | Explicitly only \f$ t \bar{t} \f$ production |
| `Mela.Production.bbH` | 5 | Explicitly only \f$ b \bar{b} \f$ production |
| `Mela.Production.JQCD` | 6 | Single-Jet QCD production |
| `Mela.Production.JJQCD` | 7 | **IDK** |
| `Mela.Production.JJVBF` | 8 | Double-Jet VBF production (ZZ/WW Fusion) |
| `Mela.Production.JJEW` | 9 | Combination of VBF and ZH/WH production |
| `Mela.Production.JJEWQCD` | 10 | Combination of JJEW and JJQCD |
| `Mela.Production.Had_ZH` | 11 | Hadronic ZH production |
| `Mela.Production.Had_WH` | 12 | Hadronic WH production |
| `Mela.Production.Lep_ZH` | 13 | Leptonic ZH production (i.e. \f$ e^+ e^- \f$ collisions) |
| `Mela.Production.Lep_WH` | 14 | Leptonic WH production |
| `Mela.Production.ZZQQB_S` | 15 | **IDK** |
| `Mela.Production.JJQCD_S` | 16 | **IDK** |
| `Mela.Production.JJVBF_S` | 17 | **IDK** |
| `Mela.Production.JJEW_S` | 18 | **IDK** |
| `Mela.Production.JJEWQCD_S` | 19 | **IDK** |
| `Mela.Production.Had_ZH_S` | 20 | **IDK** |
| `Mela.Production.Had_WH_S` | 21 | **IDK** |
| `Mela.Production.Lep_ZH_S` | 22 | **IDK** |
| `Mela.Production.Lep_WH_S` | 23 | **IDK** |
| `Mela.Production.ZZQQB_TU` | 24 | **IDK** |
| `Mela.Production.JJQCD_TU` | 25 | **IDK** |
| `Mela.Production.JJVBF_TU` | 26 | **IDK** |
| `Mela.Production.JJEW_TU` | 27 | **IDK** |
| `Mela.Production.JJEWQCD_TU` | 28 | **IDK** |
| `Mela.Production.Had_ZH_TU` | 29 | **IDK** |
| `Mela.Production.Had_WH_TU` | 30 | **IDK** |
| `Mela.Production.Lep_ZH_TU` | 31 | **IDK** |
| `Mela.Production.Lep_WH_TU` | 32 | **IDK** |
| `Mela.Production.GammaH` | 33 | **IDK** |
| `Mela.Production.nProductions` | 34 | This is the total number of production modes that exist within MELA |

One sets the production mode using `Mela.setProcess(Mela.Process, Mela.MatrixElement, Mela.Production)`, which is originally defined in Mela::setProcess.

## Process {#proc_enum}

`Mela.Process` controls what process MELA is running. This is usually the spin of the particle for JHUGen, or signal/background for MCFM.

**By and large, many of the processes are redundant, and can be re-created using some combination of couplings and production modes.**
The "minimal basis" for the processes are as follow.

* JHUGen Processes
  * SelfDefine_spin0
  * SelfDefine_spin1
  * SelfDefine_spin2
* JHUGen-MCFM Processes
  * HSMHiggs
  * bkgGammaGamma
  * bkgZGamma
  * bkgZJets
  * bkgZZ
  * bkgWW
  * bkgWWZZ
  * bkgZZ_SMHiggs
  * bkgWW_SMHiggs
  * bkgWWZZ_SMHiggs
  * HSMHiggs_WWZZ

**If you pair a process with the incorrect matrix element there will be an error!**

| Name | Value | Matrix Element | Summary |
| ---- | ----- | ------- | ------ |
| `HSMHiggs` | 0 | MCFM |  |
| `H0_g1prime2` | 0 | JHUGen |  |
| `H0hplus` | 0 | JHUGen |  |
| `H0minus` | 0 | JHUGen |  |
| `H0_Zgsg1prime2` | 0 | JHUGen |  |
| `H0_Zgs` | 0 | JHUGen |  |
| `H0_Zgs_PS` | 0 | JHUGen |  |
| `H0_gsgs` | 0 | JHUGen |  |
| `H0_gsgs_PS` | 0 | JHUGen |  |
| `D_g1g1prime2` | 0 | JHUGen |  |
| `D_g1g2` | 0 | JHUGen |  |
| `D_g1g2_pi_2` | 0 | JHUGen |  |
| `D_g1g4` | 0 | JHUGen |  |
| `D_g1g4_pi_2` | 0 | JHUGen |  |
| `D_zzzg` | 0 | JHUGen |  |
| `D_zzgg` | 0 | JHUGen |  |
| `D_zzzg_PS` | 0 | JHUGen |  |
| `D_zzgg_PS` | 0 | JHUGen |  |
| `D_zzzg_g1prime2` | 0 | JHUGen |  |
| `D_zzzg_g1prime2_pi_2` | 0 | JHUGen |  |
| `H1minus` | 0 | JHUGen |  |
| `H1plus` | 0 | JHUGen |  |
| `H2_g1` | 0 | JHUGen |  |
| `H2_g2` | 0 | JHUGen |  |
| `H2_g3` | 0 | JHUGen |  |
| `H2_g4` | 0 | JHUGen |  |
| `H2_g5` | 0 | JHUGen |  |
| `H2_g1g5` | 0 | JHUGen |  |
| `H2_g6` | 0 | JHUGen |  |
| `H2_g7` | 0 | JHUGen |  |
| `H2_g8` | 0 | JHUGen |  |
| `H2_g9` | 0 | JHUGen |  |
| `H2_g10` | 0 | JHUGen |  |
| `bkgGammaGamma` | 0 | MCFM |  |
| `bkgZGamma` | 0 | MCFM |  |
| `bkgZJets` | 0 | MCFM |  |
| `bkgZZ` | 0 | MCFM |  |
| `bkgWW` | 0 | MCFM |  |
| `bkgWWZZ` | 0 | MCFM |  |
| `bkgZZ_SMHiggs` | 0 | MCFM |  |
| `bkgWW_SMHiggs` | 0 | MCFM |  |
| `bkgWWZZ_SMHiggs` | 0 | MCFM |  |
| `HSMHiggs_WWZZ` | 0 | MCFM |  |
| `D_gg10` | 0 | JHUGen |  |
| `SelfDefine_spin0` | 0 | JHUGen |  |
| `SelfDefine_spin1` | 0 | JHUGen |  |
| `SelfDefine_spin2` | 0 | JHUGen |  |
| `nProcesses` | 0 |  | The total number of processes available to you |

One sets the process using `Mela.setProcess(Mela.Process, Mela.MatrixElement, Mela.Production)`, which is originally defined in Mela::setProcess.

## Resonance Propagator Scheme {#reso_enum}

`Mela.ResonancePropagatorScheme` controls the scheme that resonances are defined by. This is used to control mass shapes as well as conduct POWHEG high-mass reweighting. This is originally defined in TVar::ResonancePropagatorScheme. **`Mela.ResonancePropagatorScheme.FixedWidth` is set by default**.
The values are tabulated below:

| Name | Value | Summary |
| ---- | ----- | ------- |
| NoPropagator | 0 | No propagator is applied to the resonance |
| RunningWidth | 1 | A running width scheme is applied to the resonance |
| FixedWidth | 2 | A fixed width scheme is applied to the resonance |
| CPS | 3 | A Complex-Pole scheme is applied to the resonance |
| AltRunningWidth | 4 | An S-Wave Breit-Wigner is applied to the resonance |

These are provided as inputs to `Mela.getXPropagator(Mela.ResonancePropagatorScheme)` to get a value for mass shape reweighting.

## EventScaleScheme {#scale_enum}

`Mela.EventScaleScheme` controls the scaling scheme of the event. You can set it with `Mela.setRenFacScaleMode` (originally defined in Mela::setRenFacScaleMode), and get the value from `Mela.getRenFacScaleMode` (originally defined in Mela::getRenFacScaleMode).

| Name | Value | Summary |
| ---- | ----- | ------- |
| DefaultScaleScheme | 0 | The default scale factor scheme for the process/production. |
| Fixed_mH | 1 | A scale based off of the pole mass of \f$m_H\f$ |
| Fixed_mW | 2 | A scale based off of the pole mass of \f$m_Z\f$ |
| Fixed_mZ | 3 | A scale based off of the pole mass of \f$m_W\f$ |
| Fixed_mWPlusmH | 4 | A scale based off of the pole value \f$m_W + m_H\f$ |
| Fixed_mZPlusmH | 5 | A scale based off of the pole value \f$m_Z + m_H\f$ |
| Fixed_TwomtPlusmH | 6 | A scale based off of the pole value \f$2m_t + m_H\f$|
| Fixed_mtPlusmH | 7 | A scale based off of the pole value \f$m_t + m_H\f$ |
| Dynamic_qH | 8 | A scale based off of the \f$q^2\f$ of the Higgs |
| Dynamic_qJJH | 9 | A scale based off of the \f$q^2\f$ of H+2 jets |
| Dynamic_qJJ_qH | 10 | **IDK** |
| Dynamic_qJ_qJ_qH | 11 | **IDK** |
| Dynamic_HT | 12 | **IDK** |
| Dynamic_Leading_pTJ | 13 | **IDK** |
| Dynamic_Softest_pTJ | 14 | **IDK** |
| Dynamic_RandomUniform_Constrained | 15 | **IDK** |
| nEventScaleSchemes | 16 | The number of event scale schemes available to you |

## CandidateDecayMode {#decmode_enum}

`Mela.CandidateDecayMode` quantifies the decay mode of the candidate in MELA. This can be a range of things 
