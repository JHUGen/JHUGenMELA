/**
 * @file Mela.h
 * @brief This is the "MELA" object that interfaces with the Fortran code in both MCFM-JHUGen and pure JHUGen
 * 
 * @authors I. Anderson, S. Bolognesi, F. Caola, J. Davis, Y. Gao, A. V. Gritsan, 
 * @authors L. S. Mandacaru Guerra, Z. Guo, L. Kang, S. Kyriacou, C. B. Martin, T. Martini, 
 * @authors K. Melnikov, R. Pan, M. Panagiotou, R. Rontsch, J. Roskes, U. Sarica, 
 * @authors M. Schulze, M. V. Srivastav, N. V. Tran, A. Whitbeck, M. Xiao, Y. Zhou
*/

/*
************* HEADER: CMS MELA interface to MCFM/JHUGen-MELA *************
Please see the ../src/Mela.cc file for the instructions.
*/

#ifndef MELA_Mela_h
#define MELA_Mela_h

#include <vector>
#include "TLorentzVector.h"
#include "TRandom3.h"


class TFile; 
class TGraph;
class TH1F;
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class ScalarPdfFactory_HVV;
class VectorPdfFactory;
class TensorPdfFactory;
class RooqqZZ_JHU_ZgammaZZ_fast;
class ZZMatrixElement;
class SuperMELA;

#include "TVar.hh"
#include "TEvtProb.hh"
#include "MelaPConstant.h"
#include "SuperDijetMela.h"
#include "ScalarPdfFactory_HVV.h"
#include "VectorPdfFactory.h"
#include "TensorPdfFactory_ppHVV.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"

class Mela{

public:
  /**
   * @brief the MELA constructor
   * 
   * @param[in] LHCsqrts_ The luminosity of your collider in TeV. Default is LHC at 13 TeV.
   * @param[in] mh_ The mass of the Higgs in GeV. Default is 125 GeV
   * @param[in] verbosity_ The verbosity of MELA that you desire, as defined in TVar::ERROR
  */
  Mela(double LHCsqrts_=13., double mh_=125., TVar::VerbosityLevel verbosity_=TVar::ERROR); // Higgs mass for supermela

  /**
   * @brief Copy constructor for MELA
   * 
   * @param[in] other another MELA instance
  */
  Mela(const Mela& other);

  /**
   * @brief MELA destructor
  */
  ~Mela();

  /**
   * @brief This is the actual building of the tool that occurs in each instance of the Mela::Mela constructor.
   * 
   * @param mh_ This is the mass of the Higgs in GeV
  */
  void build(double mh_);

  /**
   * @brief Sets the process, matrix element, and production that MELA is to use for this event. Calls ZZMatrixElement::set_Process, which calls TEvtProb::SetProcess.
   * 
   * @attention Remember to set the process for each event, otherwise the MELA event loop will throw a segmentation error.
   * @param[in] myModel a TVar for the Process you would like, as defined in TVar::Process
   * @param[in] myME a TVar for the matrix element you would like, as defined in TVar::MatrixElement
   * @param[in] myProduction a TVar for the production mode you would like, as defined in TVar::Production
  */
  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);

  /**
   * @brief Sets the verbosity for MELA outside of the initial constructor.
   * 
   * @param[in] verbosity_ The verbosity of MELA that you desire, as defined in TVar::ERROR.
  */
  void setVerbosity(TVar::VerbosityLevel verbosity_=TVar::ERROR);

  /**
   * @brief Gets the current verbosity level for MELA.
   * 
   * @return a TVar::VerbosityLevel describing the verbosity level for MELA. This is a number from 0 to 5, and corresponds to the values in TVar::VerbosityLevel.
  */
  TVar::VerbosityLevel getVerbosity();

  /**
   * @brief Sets the MELA Lepton Interference
   * 
   * @param myLepInterf sets the myLepInterf_ variable to one found in TVar::LeptonInterface through ZZMatrixElement::set_LeptonInterface. TVar::DefaultLeptonInterf by default.
  */
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);

  /**
   * @brief either permits or forbids massive leptons.
   * @sa Wrapper for the function TUtil::applyLeptonMassCorrection
   * 
   * @param MasslessLeptonSwtich Whether you would like your leptons to be massive, by default true.
  */
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);

  /**
   * @brief either permits or forbids massive jets.
   * @sa Wrapper for the function TUtil::applyJetMassCorrection
   * 
   * @param MasslessLeptonSwitch Whether you would like your jets to be massive, by default true.
  */
  void setRemoveJetMasses(bool MasslessLeptonSwitch=true);

  /**
   * @brief Sets the mass of the "primary" higgs.
   * @sa Wrapper for the function ZZMatrixElement::set_PrimaryHiggsMass, which is a wrapper for TEvtProb::SetPrimaryHiggsMass
   * 
   * @param myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * 
   * @attention The primary Higgs is the first resonance - nominally MELA can have 2 resonances when working through JHUGen-MCFM
   * @remark This function is effectively the same as setMelaHiggsMass(mass, 0)
  */
  void setMelaPrimaryHiggsMass(double myHiggsMass);

  /**
   * @brief Sets the mass of your chosen Higgs. 
   * @sa Wrapper for ZZMatrixElement::set_mHiggs.
   * 
   * @param myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * @param index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsMass(double myHiggsMass, int index=0);

  /**
   * @brief Sets the width of your chosen Higgs. 
   * @sa Wrapper for ZZMatrixElement::set_wHiggs
   * 
   * @param myHiggsWidth This is the mass of the Higgs that you would like, in GeV
   * @param index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);

  /**
   * @brief a combination of setMelaHiggsMass and setMelaHiggsWidth. 
   * @sa Wrapper for ZZMatrixElement::set_mHiggs_wHiggs.
   * 
   * @param myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * @param myHiggsWidth This is the mass of the Higgs that you would like, in GeV
   * @param index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsMassWidth(double myHiggsMass, double myHiggsWidth, int index);

  /**
   * @brief Sets the renormalization and the factorization schemes. 
   * @sa Wrapper for ZZMatrixElement::set_RenFacScaleMode, which is a wrapper for TEvtProb::SetRenFacScaleMode, which edits TEvtProb.event_scales
   * 
   * @param renormalizationSch This is the renormalization scheme that you are picking from TVar::EventScaleScheme
   * @param factorizationSch This is the factorization scheme that you are picking from TVar::EventScaleScheme
   * @param ren_sf This is the renormalization scale factor that you would like
   * @param  fac_sf This is the scale factor for the factorization scale that you would like
  */
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);

  /**
   * @brief Sets the decay mode for your event. 
   * @sa Wrapper for ZZMatrixElement::set_CandidateDecayMode, which is a wrapper for TEvtProb::SetCandidateDecayMode
   * 
   * @param mode The decay mode you would like picked from TVar::CandidateDecayMode
  */
  void setCandidateDecayMode(TVar::CandidateDecayMode mode);

  /**
   * @brief Switches the candidate that you are working on to another candidate based off of an index. 
   * @sa Wrapper for ZZMatrixElement::set_CurrentCandidateFromIndex, which is a wrapper for TEvtProb::SetCurrentCandidateFromIndex.
   * 
   * @param icand The index of the candidate that you would like to switch to
   * 
  */
  void setCurrentCandidateFromIndex(unsigned int icand); // Switches to another candidate

  /**
   * @brief Switches the candidate that you are working on to another candidate object specified
   * @sa Wrapper for ZZMatrixElement::set_CurrentCandidate, which is a wrapper for TEvtProb::SetCurrentCandidate.
   * 
   * @param cand The MELACandidate object you would like to switch to
  */
  void setCurrentCandidate(MELACandidate* cand); // Switches to another candidate

  /**
   * @brief Sets the input event for MELA. MELA cannot run without this.
   * @sa Wrapper for ZZMatrixElement::set_InputEvent, which is a wrapper for TEvtProb::SetInputEvent, which calls TUtil::ConvertVectorFormat
   * @attention An input event must be set for each event in an event loop, otherwise MELA will throw a segmentation error
   * 
   * @param pDaughters A SimpleParticleCollection_t of particle daughters (decay products)
   * @param pAssociated A SimpleParticleCollection_t of associated particles (i.e. jets), by default 0 (no jets)
   * @param pMothers A SimpleParticleCollection_t of particle mothers (i.e. gluons), by default 0 (reco data contains no mother information)
   * @param isGen A boolean signifying whether the event in question is a Gen event or a reco event, by default false (reco)
  */
  void setInputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds another candidate
  
  /**
   * @brief Resets the event in preparation for the next iteration of the event loop.
   * @sa Wrapper for ZZMatrixElement::reset_InputEvent, which is a wrapper for TEvtProb::ResetInputEvent
   * @attention Without resetting the input event at the end of each for loop, behavior could be unexpected!
   * @attention It is important to call this at the end of every event loop iteration to clean up TEvtProb!
  */
  void resetInputEvent(); // Reset the input candidates. Important to call in order to clean up TEvtProb!

  /**
   * @brief Sets a temporary MELA candidate, by setting melaCand in Xcal2 to a temporary candidate without pushing this candidate to the candList of Xcal2
   * @sa Wrapper for ZZMatrixElement::set_TempCandidate 
   * 
   * @param pDaughters A SimpleParticleCollection_t of particle daughters (decay products)
   * @param pAssociated A SimpleParticleCollection_t of associated particles (i.e. jets), by default 0 (no jets)
   * @param pMothers A SimpleParticleCollection_t of particle mothers (i.e. gluons), by default 0 (reco data contains no mother information)
   * @param isGen A boolean signifying whether the event in question is a Gen event or a reco event, by default false (reco)
  */
  void setTempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds a temp. candidate
  
  /**
   * @brief Adds a top quark as a MELA candidate
   * @sa Wrapper for ZZMatrixElement::append_TopCandidate, which is a wrapper for TEvtProb::AppendTopCandidate
   * 
   * @param TopDaughters A SimpleParticleCollection_t of the top's daughters, or decay products
  */
  void appendTopCandidate(SimpleParticleCollection_t* TopDaughters); // Adds a top

  // Function to set EW parameters in MCFM/JHUGen
  /**
   * @brief Resets the mass for a particle that is an electroweak parameter according to its id
   * @sa Wrapper for ZZMatrixElement::reset_Mass, which is a wrapper for TEvtProb::ResetMass, which is a wrapper for TUtil::SetMass
   * @sa This finally interfaces with either the function SetMass in mod_parameters.F90 for JHUGen, accessed through .TModParameters.hh
   * @sa in mod_parameters.F90 the conversion from ipart to the values in JHUGen's internal code happens with the function convertLHEreverse(Part)
   * @sa or the masses in the JHUGen-MCFM library at /src/Inc/masses.F or at TMCFM::spinzerohiggs_anomcoupl (for tPrime and bPrime)
   * 
   * @warning
   * It is not recommended to use this function to edit the mass of the Higgs in MCFM.
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * The particles that you can change when using the MCFM matrix element are as as follow:
   * absolute value of ipart | Particle Name
   * :---------: | :-----------:
   * 8 | t Prime quark (4th generation)
   * 7 | b Prime quark (4th generation)
   * 6 | Top Quark
   * 5 | Bottom Quark
   * 4 | Charm Quark
   * 3 | Strange Quark
   * 2 | Up Quark
   * 1 | Down Quark
   * 11 | Electron
   * 13 | Muon
   * 15 | Tau
   * 23 | Z Boson
   * 24 | W Boson
   * 25 | Higgs Boson
   * 
   * @attention The particles that you can change when using the JHUGen matrix element are as follow:
   * ipart | Particle Name | Sign Sensitive
   * :---------: | :-----------: | :-----------:
   * 0 OR 21 | Gluon | No
   * 1 | Down Quark | Yes
   * 2 | Up Quark | Yes
   * 3 | Strange Quark | Yes
   * 4 | Charm Quark | Yes
   * 5 | Bottom Quark | Yes
   * 6 | Top Quark | Yes
   * 11 | Electron | Yes
   * 22 | Photon | No
   * 23 | Z Boson | No
   * 24 | W Boson | Yes
   * 13 | Muon | Yes
   * 15 | Tau | Yes
   * 12 | Electron Neutrino | Yes
   * 14 | Muon Neutrino | Yes
   * 16 | Tau Neutrino | Yes
   * 25 | Higgs Boson | No
   * 32 | Z Prime | No
   * 33 | Z Prime 2 | No
   * 34 | W Prime | Yes
   * 
   * @param inmass the mass that you want in GeV
   * @param ipart the particle whose mass you wish to change
   * 
  */
  void resetMass(double inmass, int ipart);

  /**
   * @brief Resets the width for a particle that is an electroweak parameter according to its id
   * @sa Wrapper for ZZMatrixElement::reset_Width, which is a wrapper for TEvtProb::ResetWidth, which is a wrapper for TUtil::SetDecayWidth
   * @sa This finally interfaces with either the function SetDecayWidth in mod_parameters.F90 for JHUGen, accessed through .TModParameters.hh
   * @sa in mod_parameters.F90 the conversion from ipart to the values in JHUGen's internal code happens with the function convertLHEreverse(Part)
   * @sa or the masses in the JHUGen-MCFM library at /src/Inc/masses.F or at TMCFM::spinzerohiggs_anomcoupl (for tPrime and bPrime)
   * 
   * @warning
   * It is not recommended to use this function to edit the mass of the Higgs in MCFM.
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * The particles that you can change when using the MCFM matrix element are as as follow:
   * absolute value of ipart | Particle Name
   * :---------: | :-----------:
   * 15 | Tau
   * 23 | Z Boson
   * 24 | W Boson
   * 25 | Higgs Boson
   * 
   * @attention The particles that you can change when using the JHUGen matrix element are as follow:
   * ipart | Particle Name | Sign Sensitive
   * :---------: | :-----------: | :-----------:
   * 0 OR 21 | Gluon | No
   * 1 | Down Quark | Yes
   * 2 | Up Quark | Yes
   * 3 | Strange Quark | Yes
   * 4 | Charm Quark | Yes
   * 5 | Bottom Quark | Yes
   * 6 | Top Quark | Yes
   * 11 | Electron | Yes
   * 22 | Photon | No
   * 23 | Z Boson | No
   * 24 | W Boson | Yes
   * 13 | Muon | Yes
   * 15 | Tau | Yes
   * 12 | Electron Neutrino | Yes
   * 14 | Muon Neutrino | Yes
   * 16 | Tau Neutrino | Yes
   * 25 | Higgs Boson | No
   * 32 | Z Prime | No
   * 33 | Z Prime 2 | No
   * 34 | W Prime | Yes
   * 
   * @param inmass the mass that you want in GeV
   * @param ipart the particle whose mass you wish to change
   * 
  */
  void resetWidth(double inwidth, int ipart);

  /**
   * @brief Resets the masses of each quark to their original values.
   * @sa Wrapper for ZZMatrixElement::reset_QuarkMasses, which is a wrapper for TEvtProb::ResetQuarkMasses, which calls TEvtProb::ResetMass
   * The quarks that are reset are as follow, with their default masses in GeV listed as well:
   * 
   * Quark | Mass 
   * Down | 0.001
   * Up | 0.005
   * Strange | 0.1
   * Charm | 1.275
   * Bottom | 4.75
   * Top | 173.2
   * B Prime | 1e5
   * T Prime | 1e5
  */
  void resetQuarkMasses();
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  // Function to get current primary EW/QCD parameters from MCFM/JHUGen (notice Higgs mass/width used in the ME could be different)
  double getPrimaryMass(int ipart);
  double getPrimaryWidth(int ipart);
  double getHiggsWidthAtPoleMass(double mass);


  MelaIO* getIORecord(); // Full parton-by-parton ME record
  MELACandidate* getCurrentCandidate();
  int getCurrentCandidateIndex();
  int getNCandidates();
  std::vector<MELATopCandidate_t*>* getTopCandidateCollection();


  void getConstant(float& prob); // <ME> constants
  void getPAux(float& prob); // SuperProb


  RooSpin::modelMeasurables getMeasurablesRRV();


  void computeDecayAngles(
    float& qH,
    float& m1,
    float& m2,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
    );
  void computeVBFAngles(
    float& Q2V1,
    float& Q2V2,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );
  void computeVBFAngles_ComplexBoost(
    float& Q2V1,
    float& Q2V2,
    float& costheta1_real, float& costheta1_imag,
    float& costheta2_real, float& costheta2_imag,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );
  void computeVHAngles(
    float& mVstar,
    float& mV,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );
  void computeTTHAngles(
    int topDecay,

    float& mT1,
    float& mW1,
    float& mT2,
    float& mW2,

    // TTH system
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1,

    // TT system
    float& hbb,
    float& hWW,
    float& Phibb,
    float& Phi1bb,

    // Wplus system
    float& hWplusf,
    float& PhiWplusf,

    // Wminus system
    float& hWminusf,
    float& PhiWminusf
  );

  void computeP_selfDspin0(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZqqcoupl_input[SIZE_ZQQ][2],
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGqqcoupl_input[SIZE_GQQ][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP(
    float& prob,
    bool useConstant=true
    );

  void computeD_CP(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  //****VVH Spin-0****//
  void computeProdDecP(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDaTQGCcoupl_input[SIZE_ATQGC][2],
    double selfDAZffcoupl_input[SIZE_AZff][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );

  //****HJ/HJJ/VBF Spin-0****//
  void computeProdP(
    double selfDHggcoupl_input[SIZE_HGG][2],
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdP(
    float& prob,
    bool useConstant=true
    );

  //****VH Spin-0****//
  void computeProdP_VH(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );
  void computeProdP_VH(
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );

  //***ttH Spin-0****//
  void computeProdP_ttH(
    float& prob,
    int topProcess=2,
    int topDecay=0,
    bool useConstant=true
    );

  // Calculation weight to correct for fermion interference
  void compute4FermionWeight(float& w);

  // Calculation of X propagator
  void getXPropagator(TVar::ResonancePropagatorScheme scheme, float& prop);

  //*** SuperMela ***//
  void computePM4l(
    TVar::SuperMelaSyst syst,
    float& prob
    );

  //*** SuperJJMela ***//
  void computeDijetConvBW(float& prob, bool useTrueBW=false);

  //*** Dgg10 ***//
  void computeD_gg(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  // Access ZZMEs Calculate4Momentum
  std::vector<TLorentzVector> calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  /********************/
  /*** Data members ***/
  /********************/
  TRandom3 melaRandomNumber; // Used in SuperMELA smearing
  RooRealVar* mzz_rrv;
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv;
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* Y_rrv;
  RooRealVar* upFrac_rrv;

  RooAbsPdf* pdf;
  ScalarPdfFactory_HVV* ggSpin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory_ppHVV* spin2Model;
  RooqqZZ_JHU_ZgammaZZ_fast* qqZZmodel;

  SuperMELA* super;

  // Self-define arrays are now members of MELA.
  // There are a lot of them!
  //****Spin-0****//
  // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
  double selfDHggcoupl[nSupportedHiggses][SIZE_HGG][2];
  double selfDHg4g4coupl[nSupportedHiggses][SIZE_HGG][2];
  double selfDHqqcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHbbcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHttcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHb4b4coupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHt4t4coupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHzzLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  double selfDHwwLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  int selfDHzzCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ];
  int selfDHwwCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ];
  bool differentiate_HWW_HZZ;
  double selfDHzzpcoupl[SIZE_HVV][2];
  double selfDHzpzpcoupl[SIZE_HVV][2];
  double selfDZpffcoupl[SIZE_Vpff][2];
  double selfDHwwpcoupl[SIZE_HVV][2];
  double selfDHwpwpcoupl[SIZE_HVV][2];
  double selfDWpffcoupl[SIZE_Vpff][2];
  double selfDM_Zprime;
  double selfDGa_Zprime;
  double selfDM_Wprime;
  double selfDGa_Wprime;
  //****Spin-1****//
  double selfDZqqcoupl[SIZE_ZQQ][2];
  double selfDZvvcoupl[SIZE_ZVV][2];
  //****Spin-2****//
  double selfDGqqcoupl[SIZE_GQQ][2];
  double selfDGggcoupl[SIZE_GGG][2];
  double selfDGvvcoupl[SIZE_GVV][2];
  double selfDGvvpcoupl[SIZE_GVV][2];
  double selfDGvpvpcoupl[SIZE_GVV][2];
  //****aTQGC****//
  double selfDaTQGCcoupl[SIZE_ATQGC][2];
  //****Anomnalous Zff**//
  double selfDAZffcoupl[SIZE_AZff][2];
  // That is a lot of them!

  static void cleanLinkedFiles();

protected:
  /********************/
  /*** Data members ***/
  /********************/
  double LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  TVar::LeptonInterference myLepInterf_;
  TVar::VerbosityLevel myVerbosity_;

  ZZMatrixElement* ZZME;
  SuperDijetMela* superDijet;


  float auxiliaryProb;

  MELACandidate* melaCand; // Pointer to persistent TEvtProb object

  /***** ME CONSTANT HANDLES *****/
  // Constants that vary with sqrts due to application of PDFs
  //
  MelaPConstant* pAvgSmooth_JHUGen_JQCD_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_JJQCD_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_JJVBF_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_Had_ZH_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_Had_WH_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  // Decay ME constants that do not use PDFs
  //
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZJets_2l2q;


  /*****************/
  /*** Functions ***/
  /*****************/
  void printLogo() const;

  void setSpinZeroCouplings();
  void setSpinOneCouplings();
  void setSpinTwoCouplings();
  void setATQGCCouplings();
  void setAZffCouplings();

  bool configureAnalyticalPDFs();
  void reset_SelfDCouplings();
  void reset_PAux(); // SuperProb reset
  void reset_CandRef();

  void constructDggr(
    float bkg_VAMCFM_noscale,
    float ggzz_VAMCFM_noscale,
    float ggHZZ_prob_pure_noscale,
    float ggHZZ_prob_int_noscale,
    float widthScale,
    float& myDggr
    );

  void getPConstantHandles();
  void deletePConstantHandles();
  MelaPConstant* getPConstantHandle(
    TVar::MatrixElement me_,
    TVar::Production prod_,
    TVar::Process proc_,
    TString relpath,
    TString spname,
    const bool useSqrts=false
    );
  void deletePConstantHandle(MelaPConstant*& handle);
  void computeConstant(float& prob);
  void setConstant();
  float getConstant_JHUGenUndecayed();
  float getConstant_4l();
  float getConstant_2l2q();
  float getConstant_4q();
  float getConstant_FourFermionDecay(const int& decid);

};

#endif

