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
   * @param[in] myLepInterf sets the myLepInterf_ variable to one found in TVar::LeptonInterface through ZZMatrixElement::set_LeptonInterface. TVar::DefaultLeptonInterf by default.
  */
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);

  /**
   * @brief either permits or forbids massive leptons.
   * @sa Wrapper for the function TUtil::applyLeptonMassCorrection
   * 
   * @param[in] MasslessLeptonSwtich Whether you would like your leptons to be massive, by default true.
  */
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);

  /**
   * @brief either permits or forbids massive jets.
   * @sa Wrapper for the function TUtil::applyJetMassCorrection
   * 
   * @param[in] MasslessLeptonSwitch Whether you would like your jets to be massive, by default true.
  */
  void setRemoveJetMasses(bool MasslessLeptonSwitch=true);

  /**
   * @brief Sets the mass of the "primary" higgs.
   * @sa Wrapper for the function ZZMatrixElement::set_PrimaryHiggsMass, which is a wrapper for TEvtProb::SetPrimaryHiggsMass
   * 
   * @param[in] myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * 
   * @attention The primary Higgs is the first resonance - nominally MELA can have 2 resonances when working through JHUGen-MCFM
   * @remark This function is effectively the same as setMelaHiggsMass(mass, 0)
  */
  void setMelaPrimaryHiggsMass(double myHiggsMass);

  /**
   * @brief Sets the mass of your chosen Higgs. 
   * @sa Wrapper for ZZMatrixElement::set_mHiggs.
   * 
   * @param[in] myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * @param[in] index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsMass(double myHiggsMass, int index=0);

  /**
   * @brief Sets the width of your chosen Higgs. 
   * @sa Wrapper for ZZMatrixElement::set_wHiggs
   * 
   * @param[in] myHiggsWidth This is the mass of the Higgs that you would like, in GeV
   * @param[in] index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);

  /**
   * @brief a combination of setMelaHiggsMass and setMelaHiggsWidth. 
   * @sa Wrapper for ZZMatrixElement::set_mHiggs_wHiggs.
   * 
   * @param[in] myHiggsMass This is the mass of the Higgs that you would like, in GeV
   * @param[in] myHiggsWidth This is the mass of the Higgs that you would like, in GeV
   * @param[in] index This is either 0 or 1, depending on if you want to change the first or second resonance. By default it is 0.
  */
  void setMelaHiggsMassWidth(double myHiggsMass, double myHiggsWidth, int index);

  /**
   * @brief Sets the renormalization and the factorization schemes. 
   * @sa Wrapper for ZZMatrixElement::set_RenFacScaleMode, which is a wrapper for TEvtProb::SetRenFacScaleMode, which edits the TVar::event_scales_type struct in TEvtProb.event_scales
   * 
   * @param[in] renormalizationSch This is the renormalization scheme that you are picking from TVar::EventScaleScheme
   * @param[in] factorizationSch This is the factorization scheme that you are picking from TVar::EventScaleScheme
   * @param[in] ren_sf This is the renormalization scale factor that you would like
   * @param[in]  fac_sf This is the scale factor for the factorization scale that you would like
  */
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);

  /**
   * @brief Sets the decay mode for your event. 
   * @sa Wrapper for ZZMatrixElement::set_CandidateDecayMode, which is a wrapper for TEvtProb::SetCandidateDecayMode
   * 
   * @param[in] mode The decay mode you would like picked from TVar::CandidateDecayMode
  */
  void setCandidateDecayMode(TVar::CandidateDecayMode mode);

  /**
   * @brief Switches the candidate that you are working on to another candidate based off of an index. 
   * @sa Wrapper for ZZMatrixElement::set_CurrentCandidateFromIndex, which is a wrapper for TEvtProb::SetCurrentCandidateFromIndex.
   * 
   * @param[in] icand The index of the candidate that you would like to switch to
   * 
  */
  void setCurrentCandidateFromIndex(unsigned int icand); // Switches to another candidate

  /**
   * @brief Switches the candidate that you are working on to another candidate object specified
   * @sa Wrapper for ZZMatrixElement::set_CurrentCandidate, which is a wrapper for TEvtProb::SetCurrentCandidate.
   * 
   * @param[in] cand The MELACandidate object you would like to switch to
  */
  void setCurrentCandidate(MELACandidate* cand); // Switches to another candidate

  /**
   * @brief Sets the input event for MELA. MELA cannot run without this.
   * @sa Wrapper for ZZMatrixElement::set_InputEvent, which is a wrapper for TEvtProb::SetInputEvent, which calls TUtil::ConvertVectorFormat
   * @attention An input event must be set for each event in an event loop, otherwise MELA will throw a segmentation error
   * 
   * @param[in] pDaughters A SimpleParticleCollection_t of particle daughters (decay products)
   * @param[in] pAssociated A SimpleParticleCollection_t of associated particles (i.e. jets), by default 0 (no jets)
   * @param[in] pMothers A SimpleParticleCollection_t of particle mothers (i.e. gluons), by default 0 (reco data contains no mother information)
   * @param[in] isGen A boolean signifying whether the event in question is a Gen event or a reco event, by default false (reco)
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
   * @param[in] pDaughters A SimpleParticleCollection_t of particle daughters (decay products)
   * @param[in] pAssociated A SimpleParticleCollection_t of associated particles (i.e. jets), by default 0 (no jets)
   * @param[in] pMothers A SimpleParticleCollection_t of particle mothers (i.e. gluons), by default 0 (reco data contains no mother information)
   * @param[in] isGen A boolean signifying whether the event in question is a Gen event or a reco event, by default false (reco)
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
   * @param[in] TopDaughters A SimpleParticleCollection_t of the top's daughters, or decay products
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
   * It is not recommended to use this function to edit the mass of the Higgs in MCFM. Please use either setMelaHiggsMass or setMelaHiggsMassWidth for that.
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * @attention In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * @attention The particles that you can change when using the MCFM matrix element are as as follow:
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
   * @param[in] inmass the mass that you want in GeV
   * @param[in] ipart the particle whose mass you wish to change
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
   * It is not recommended to use this function to edit the mass of the Higgs in MCFM. Please use either setMelaHiggsWidth or setMelaHiggsMassWidth for that.
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * @attention In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * @attention The particles that you can change when using the MCFM matrix element are as as follow:
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
   * @param[in] inmass the mass that you want in GeV
   * @param[in] ipart the particle whose mass you wish to change
   * 
  */
  void resetWidth(double inwidth, int ipart);

  /**
   * @brief Resets the masses of each quark to their original values.
   * @sa Wrapper for ZZMatrixElement::reset_QuarkMasses, which is a wrapper for TEvtProb::ResetQuarkMasses, which calls TEvtProb::ResetMass
   * @warning You should run this command each time you would like to reset the mass of a quark to mitigate unexpected behavior
   * @attention The quarks that are reset are as follow, with their default masses in GeV listed as well:
   * Quark | Mass 
   * :---------: | :-----------: |
   * Down | 0.001
   * Up | 0.005
   * Strange | 0.1
   * Charm | 1.275
   * Bottom | 4.75
   * Top | 173.2
   * B Prime | 1e5
   * T Prime | 1e5
   * 
  */
  void resetQuarkMasses();

  /**
   * @brief Resets the electroweak parameters back to their defaults
   * @sa Wrapper for ZZMatrixElement::reset_MCFM_EWKParameters, which is a wrapper for TEvtProb::ResetMCFM_EWKParameters, which calls TUtil::SetEwkCouplingParameters
  */
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);

  // Function to get current primary EW/QCD parameters from MCFM/JHUGen (notice Higgs mass/width used in the ME could be different)
  /**
   * @brief A function to get the current primary EW/QCD parameters from MELA
   * @sa Wrapper for ZZMatrixElement::get_PrimaryMass, which is a Wrapper for TEvtProb::GetPrimaryMass, which calls TUtil::GetMass.
   * @warning (Higgs mass/width used in the Matrix Element could be different). Refer to setMelaHiggsMassWidth or setMELAHiggsMass for explicit references to those.
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * @attention In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * @attention The particles that you can change when using the MCFM matrix element are as as follow:
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
   * @param[in] ipart the particle you would like to get the mass of.
  */
  double getPrimaryMass(int ipart);

  /**
   * @brief A function to get the current primary EW/QCD parameters from MELA
   * @sa Wrapper for ZZMatrixElement::get_PrimaryMass, which is a Wrapper for TEvtProb::GetPrimaryMass, which calls TUtil::GetMass.
   * @warning (Higgs mass/width used in the Matrix Element could be different). Refer to setMelaHiggsMassWidth or setMELAHiggsMass for explicit references to those.
   * 
   * 
   * @attention You can input negative ids to edit antiparticles (i.e. 24 -> W+, -24 -> W-)
   * @attention In MCFM both the particle and its antiparticle's mass are changed at the same time. In JHUGen this is not the case - they are separate entries.
   * @attention The particles that you can change when using the MCFM matrix element are as as follow:
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
   * @param[in] ipart the particle you would like to get the width of.
  */
  double getPrimaryWidth(int ipart);

  /**
   * @brief Returns the width of the Higgs at a given pole mass as a calculation
   * @sa Wrapper for ZZMatrixElement::get_HiggsWidthAtPoleMass, which is a wrapper for TEvtProb::GetHiggsWidthAtPoleMass, which is a wrapper for MELAHXSWidth::HiggsWidth
   * 
   * @attention This function does an independent calculation for the mass of the Higgs at a certain mass, and does not rely on any other set MELA values
   * 
   * @param[in] mass The value of the Higgs pole mass you want in GeV
  */
  double getHiggsWidthAtPoleMass(double mass);

  /**
   * @brief Returns the MELAIO object, and by consequence, the entire parton-by-parton matrix element record
   * @sa wrapper for ZZMatrixElement::get_IORecord, which is a wrapper for TEvtProb::GetIORecord, which is a wrapper for MelaIO::getRef
  */
  MelaIO* getIORecord(); // Full parton-by-parton ME record

  /**
   * @brief Gets the current MELA top-level (input) candList object
   * @sa Wrapper for ZZMatrixElement::get_CurrentCandidate, which is a wrapper for TEvtProb::GetCurrentCandidate, which returns the protected element TEvtProb::melaCand
   *  
  */
  MELACandidate* getCurrentCandidate();

  /**
   * @brief Returns the index of the current MELA candidate - returns -1 if there is no candidate to be found
   * @sa Wrapper for ZZMatrixElement::get_CurrentCandidateIndex, which is a wrapper for TEvtProb::GetCurrentCandidateIndex
  */
  int getCurrentCandidateIndex();

  /**
   * @brief Returns the size of the candidate list TEvtProb::candList.
   * @sa Wrapper for ZZMatrixElement::get_NCandidates, which is a wrapper for TEvtProb::GetNCandidates, which returns the size of TEvtProb::candList
  */
  int getNCandidates();

  /**
   * @brief Same as getNCandidates but specifically for Top Quark Candidates
   * @sa Wrapper for ZZMatrixElement::get_TopCandidateCollection, which is a wrapper for TEvtProb::GetTopCandidates, which returns TEvtProb::topCandList
  */
  std::vector<MELATopCandidate_t*>* getTopCandidateCollection();

  /**
   * @brief This function returns a multiplicative constant to the matrix element calculation 
   * @attention if useConstant=True in the following functions the constant will be applied to the matrix element:
   * - computeP_selfDspin0
   * - computeP_selfDspin1
   * - computeP_selfDspin1
   * - computeP_selfDspin2
   * - computeP_selfDspin2
   * - computeP
   * - computeProdDecP
   * - computeProdDecP
   * - computeProdP
   * - computeProdP
   * - computeProdP_VH
   * - computeProdP_VH
   * - computeProdP_ttH
   * - computePM4l
   * 
   * @param[out] prob This is the number that will be multiplied _inplace_ by the constant.
  */
  void getConstant(float& prob); // <ME> constants


  void getPAux(float& prob); // SuperProb

  /**
   * @brief Returns a RooSpin::modelMeasureables object containing many observable quantities
   * 
   * @attention The parameters in question are:
   * - $\cos(\theta_1)$
   * - $\cos(\theta_2)$
   * - $\phi$
   * - mass of the first Z boson ($M_{Z1}$)
   * - mass of the second Z boson ($M_{Z2}$)
   * - mass of the full particle ($M_{ZZ}$)
   * - $\cos(\theta*)$
   * - $\phi_1$
   * - Y
   * 
   * @remark these angles are defined in Figure 1 of [arXiv:1001.3396](https://arxiv.org/abs/1001.3396) as well as (described more generally) in Figure 1 of [arXiv:1208.4018](https://arxiv.org/abs/1208.4018)
  */
  RooSpin::modelMeasurables getMeasurablesRRV();

  /**
   * @brief computes the decay angles for gg -> H -> ZZ as defined in Figure 1 of [arXiv:1001.3396](https://arxiv.org/abs/1001.3396)
   * @sa Calls TUtil::computeAngles to fully calculate angles
   * @sa Use TUtil::computeAngles if you would like to calculate the angles through inputting your own 4-vectors and quantities rather than through the MELA event loop
   * @attention This function requires there to be an input event already defined
   * @warning This function will edit all values *inplace*. Every value is technically an output.
   * @warning The selection for m1 and m2 are based upon which of the leptons is closer to the z mass. $m_1 \leq m_2$.
   * 
   * @param[out] qH The mass of the "Higgs" as reconstructed by the constituent 4 leptons
   * @param[out] m1 The mass of the first decay particle as reconstructed by 2 of the leptons
   * @param[out] m2 The mass of the first decay particle as reconstructed by 2 of the leptons
   * @param[out] costheta1 $\cos(\theta_1)$
   * @param[out] costheta2 $\cos(\theta_2)$
   * @param[out] Phi $\phi$
   * @param[out] costhetastar $\cos(\theta*)$
   * @param[out] Phi1 $\phi_1$
  */
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
  
  /**
   * @brief computes the decay angles for Vector Boson Fusion (VBF) -> H -> ZZ as defined in Figure 1 of [arXiv:1208.4018](https://arxiv.org/abs/1208.4018)
   * @sa Calls TUtil::computeVBFAngles to fully calculate angles
   * @sa Use TUtil::computeVBFAngles if you would like to calculate the angles through inputting your own 4-vectors and quantities rather than through the MELA event loop
   * @attention This function requires there to be an input event already defined
   * @warning This function will edit all values *inplace*. Every value is technically an output.
   * @warning The selection for m1 and m2 are based upon which of the leptons is closer to the z mass. $m_1 \leq m_2$.
   * 
   * @param[out] Q2V1 The mass of the first decay particle as reconstructed by 2 of the leptons
   * @param[out] m2V2 The mass of the first decay particle as reconstructed by 2 of the leptons
   * @param[out] costheta1 $\cos(\theta_1)$
   * @param[out] costheta2 $\cos(\theta_2)$
   * @param[out] Phi $\phi$
   * @param[out] costhetastar $\cos(\theta*)$
   * @param[out] Phi1 $\phi_1$
  */
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

  /**
   * @brief computes the decay angles for Vector Boson Fusion (VBF) -> H -> ZZ as defined in Figure 1 of [arXiv:1208.4018](https://arxiv.org/abs/1208.4018)
   * @sa Calls TUtil::computeVHAngles to fully calculate angles
   * @sa Use TUtil::computeVHAngles if you would like to calculate the angles through inputting your own 4-vectors and quantities rather than through the MELA event loop
   * @attention This function requires there to be an input event already defined
   * @warning This function will edit all values *inplace*. Every value is technically an output.
   * @warning The selection for m1 and m2 are based upon which of the leptons is closer to the z mass. $m_1 \leq m_2$.
   * 
   * @param[out] mVstar The mass of the virtual Z boson in ZH production
   * @param[out] mV The mass of the Z in decay
   * @param[out] costheta1 $\cos(\theta_1)$
   * @param[out] costheta2 $\cos(\theta_2)$
   * @param[out] Phi $\phi$
   * @param[out] costhetastar $\cos(\theta*)$
   * @param[out] Phi1 $\phi_1$
  */
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

  /**
   * @brief This is a function that calls Mela::computeP with preset quark, gluon, and vector boson couplings for spin 0
   * @sa Mela::computeP, Mela::selfDHqqcoupl, Mela::selfDHggcoupl, Mela::selfDHttcoupl, Mela::selfDHbbcoupl, Mela::selfDHzzcoupl, Mela::selfDHwwcoupl
   * @warning I (Mohit Srivastav) would personally avoid using this function, as it obfuscates the couplings you are <em> really using </em>.
   * @warning It is always a good idea to set couplings yourself so that you know what you're really doing.
   * @param[in] selfDHvvcoupl_input an array corresponding to the couplings of Mela::selfDHzzcoupl (between the Higgs and the Z Boson), Mela::selfDHwwcoupl using the indices in TCouplingsBase.hh
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeP_selfDspin0(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  
  /**
   * @brief This is a function that calls Mela::computeP with preset quark, gluon, and vector boson couplings for spin 1
   * @sa Mela::computeP, Mela::selfDHqqcoupl, Mela::selfDHggcoupl, Mela::selfDHttcoupl, Mela::selfDHbbcoupl, Mela::selfDHzzcoupl, Mela::selfDHwwcoupl
   * @warning I (Mohit Srivastav) would personally avoid using this function, as it obfuscates the couplings you are <em> really using </em>.
   * @warning It is always a good idea to set couplings yourself so that you know what you're really doing.
   * @param[in] selfDZqqcoupl_input an array corresponding to the couplings of Mela::selfDZqqcoupl (between the Higgs and quarks) using the indices in TCouplingsBase.hh
   * @param[in] selfDZvvcoupl_input an array corresponding to the couplings of Mela::selfDZvvcoupl (between the Higgs and the Z/W Bosons) using the indices in TCouplingsBase.hh
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeP_selfDspin1(
    double selfDZqqcoupl_input[SIZE_ZQQ][2],
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  
  /**
   * @brief This is a function that calls Mela::computeP with preset quark, gluon, and vector boson couplings for spin 1, but sets default values for Zqq couplings
   * @sa Mela::computeP, Mela::selfDHqqcoupl, Mela::selfDHggcoupl, Mela::selfDHttcoupl, Mela::selfDHbbcoupl, Mela::selfDHzzcoupl, Mela::selfDHwwcoupl
   * @warning I (Mohit Srivastav) would personally avoid using this function, as it obfuscates the couplings you are <em> really using </em>.
   * @warning It is always a good idea to set couplings yourself so that you know what you're really doing.
   * @param[in] selfDZvvcoupl_input an array corresponding to the couplings of Mela::selfDZvvcoupl (between the Higgs and the Z/W Bosons) using the indices in TCouplingsBase.hh
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeP_selfDspin1(
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );

  /**
   * @brief This is a function that calls Mela::computeP with preset quark, gluon, and vector boson couplings for spin 2
   * @sa Mela::computeP, Mela::selfDGggcoupl, Mela::selfDGqqcoupl, Mela::selfDGvvcoupl
   * @warning I (Mohit Srivastav) would personally avoid using this function, as it obfuscates the couplings you are <em> really using </em>.
   * @warning It is always a good idea to set couplings yourself so that you know what you're really doing.
   * @param[in] selfDGggcoupl_input an array corresponding to the couplings of Mela::selfDGggcoupl (between the "graviton" and gluons) using the indices in TCouplingsBase.hh
   * @param[in] selfDGqqcoupl_input an array corresponding to the couplings of Mela::selfDGqqcoupl (between the "graviton" and quarks) using the indices in TCouplingsBase.hh
   * @param[in] selfDGvvcoupl_input an array corresponding to the couplings of Mela::selfDGvvcoupl (between the "graviton" and the W/Z bosons) using the indices in TCouplingsBase.hh
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGqqcoupl_input[SIZE_GQQ][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );

  /**
   * @brief This is a function that calls Mela::computeP with preset quark, gluon, and vector boson couplings for spin 1, but sets default values for Zqq couplings
   * @sa Mela::computeP, Mela::selfDGggcoupl, Mela::selfDGqqcoupl, Mela::selfDGvvcoupl
   * @warning I (Mohit Srivastav) would personally avoid using this function, as it obfuscates the couplings you are <em> really using </em>.
   * @warning It is always a good idea to set couplings yourself so that you know what you're really doing.
   * @param[in] selfDGggcoupl_input an array corresponding to the couplings of Mela::selfDGggcoupl (between the "graviton" and gluons) using the indices in TCouplingsBase.hh
   * @param[in] selfDGvvcoupl_input an array corresponding to the couplings of Mela::selfDGvvcoupl (between the "graviton" and the W/Z bosons) using the indices in TCouplingsBase.hh
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  
  /**
   * @brief Computes the probability for the probabilities on the decay side of things using the constituent daughter 4 vectors along with any jets that could exist
   * @sa Mela::computeDecayAngles
   * @attention All JHUGen matrix elements, productions, and processes can be used here
   * @attention The following production modes are supported using the MCFM Matrix Element (as named in TVar::Production)
   * TVar::ZZGG
   * TVar::ZZQQB
   * TVar::ZZQQB_STU
   * TVar::ZZQQB_S
   * TVar::ZZQQB_TU
   * TVar::ZZINDEPENDENT
   * TVar::JJQCD && process == TVar::bkgZJets
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
   * 
  */
  void computeP(
    float& prob,
    bool useConstant=true
    );

  /**
   * @brief computes the value of D_CP
   * @sa calls Mela::computeP_selfDspin0
   * @param[in] myME the matrix element you want to use, set as a TVar::MatrixElement
   * @param[in] myType the process you want to use, set as a TVar::Process
   * @param[out] prob This is the value of the output probability - edited by reference
  */
  void computeD_CP(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  //****VVH Spin-0****//
  /**
   * @brief computes the combined production and decay probability while taking in coupling arrays
   * @attention The JHUGen Matrix element is not supported
   * @sa Calls Mela::computeProdDecP(prob, useConstant)
   * @sa Mela::selfDHvvcoupl, Mela::selfDHwwcoupl, Mela::selfDaTQGCcoupl, Mela::selfDAZffcoupl
   * @param[in] selfDHvvcoupl_input input for the couplings set in Mela::selfDHvvcoupl
   * @param[in] selfDHwwcoupl_input input for the couplings set in Mela::selfDHwwcoupl
   * @param[in] selfDaTQGCcoupl_input input for the couplings set in Mela::selfDaTQGCcoupl
   * @param[in] selfDAZffcoupl_input input for the couplings set in Mela::selfDAZffcoupl
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
   * 
  */
  void computeProdDecP(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDaTQGCcoupl_input[SIZE_ATQGC][2],
    double selfDAZffcoupl_input[SIZE_AZff][2],
    float& prob,
    bool useConstant=true
    );

  /**
   * @brief computes the combined production and decay probability, with couplings set beforehand
   * @attention The JHUGen Matrix element is not supported
   * @attention The following production modes are supported using the MCFM Matrix Element (as named in TVar::Production)
   * TVar::Had_WH
   * TVar::Had_ZH
   * TVar::Had_WH_S
   * TVar::Had_ZH_S
   * TVar::Had_WH_TU
   * TVar::Had_ZH_TU
   * TVar::Lep_ZH
   * TVar::Lep_WH
   * TVar::Lep_ZH_S
   * TVar::Lep_WH_S
   * TVar::Lep_ZH_TU
   * TVar::Lep_WH_TU
   * TVar::JJVBF
   * TVar::JJEW
   * TVar::JJEWQCD
   * TVar::JJQCD
   * TVar::JJVBF_S
   * TVar::JJEW_S
   * TVar::JJEWQCD_S
   * TVar::JJQCD_S
   * TVar::JJVBF_TU
   * TVar::JJEW_TU
   * TVar::JJEWQCD_TU
   * TVar::JJQCD_TU
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );

  //****HJ/HJJ/VBF Spin-0****//
  /**
   * @brief computes Production side probabilities while taking in coupling arrays
   * @attention The MCFM Matrix element is not supported
   * @sa Calls Mela::computeProdP(prob, useConstant)
   * @sa Mela:selfDHggcoupl, Mela::selfDHvvcoupl, Mela::selfDHwwcoupl
   * @param[in] selfDHggcoupl_input input for the couplings set in Mela::setlfDHggcoupl
   * @param[in] selfDHvvcoupl_input input for the couplings set in Mela::selfDHvvcoupl
   * @param[in] selfDHwwcoupl_input input for the couplings set in Mela::selfDHwwcoupl
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
  void computeProdP(
    double selfDHggcoupl_input[SIZE_HGG][2],
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );

    /**
   * @brief computes Production side probabilities using couplings set beforehand
   * @attention The MCFM Matrix element is not  supported
   * @param[out] prob This is the value of the output probability - edited by reference
   * @param[in] useConstant This turns on the calculation of a corrective constant to different probabilities through Mela::getConstant. If you would like the "pure" MELA calculation to be run, set useConstant to false. By default true.
  */
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

  /**
   * @defgroup coupl Coupling arrays
   * @brief Self-define arrays for MELA
   * @attention The first dimension (of size [TMCFM::nSupportedHiggses=2]) supports a second resonance present in MCFM
   * @{*/
  /** Couplings between the higgs and gluons. */
  double selfDHggcoupl[nSupportedHiggses][SIZE_HGG][2]; 
  /** Couplings between the Higgs and **PLACEHOLDER** */
  double selfDHg4g4coupl[nSupportedHiggses][SIZE_HGG][2];
  /** Couplings between the Higgs and quarks */
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

