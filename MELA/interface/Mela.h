/*
************* HEADER: CMS MELA interface to MCFM/JHUGen-MELA *************
Please see the ../src/Mela.cc file for the instructions.
*/

#ifndef MELA_Mela_h
#define MELA_Mela_h

#include <vector>
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "MadMela.h"


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

  Mela(double LHCsqrts_=13., double mh_=125., TVar::VerbosityLevel verbosity_=TVar::ERROR); // Higgs mass for supermela
  Mela(const Mela& other);
  ~Mela();

  // Constructor wrapper
  void build(double mh_);

  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  void setVerbosity(TVar::VerbosityLevel verbosity_=TVar::ERROR);
  TVar::VerbosityLevel getVerbosity();
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);
  void setRemoveJetMasses(bool MasslessLeptonSwitch=true);
  void setMelaPrimaryHiggsMass(double myHiggsMass);
  void setMelaHiggsMass(double myHiggsMass, int index=0);
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);
  void setMelaHiggsMassWidth(double myHiggsMass, double myHiggsWidth, int index);
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void setCandidateDecayMode(TVar::CandidateDecayMode mode);
  void setCurrentCandidateFromIndex(unsigned int icand); // Switches to another candidate
  void setCurrentCandidate(MELACandidate* cand); // Switches to another candidate
  void setInputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false,
    bool madMela=false
    ); // Adds another candidate
  void resetInputEvent(); // Reset the input candidates. Important to call in order to clean up TEvtProb!
  void setTempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds a temp. candidate
  void appendTopCandidate(SimpleParticleCollection_t* TopDaughters); // Adds a top

  // Function to set EW parameters in MCFM/JHUGen
  void resetMass(double inmass, int ipart);
  void resetWidth(double inwidth, int ipart);
  void resetYukawaMass(double inmass, int ipart);
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

  const TVar::event_scales_type& getRenFacScaleMode() const;

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
void computeP(
    double &prob, 
    int nhel=-1
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

  //madMELA couplings
  double& mdl_ch = madMela::params_r_.mdl_ch;
  double& mdl_chbox = madMela::params_r_.mdl_chbox;
  double& mdl_chdd = madMela::params_r_.mdl_chdd;
  double& mdl_chg = madMela::params_r_.mdl_chg;
  double& mdl_chw = madMela::params_r_.mdl_chw;
  double& mdl_chb = madMela::params_r_.mdl_chb;
  double& mdl_chwb = madMela::params_r_.mdl_chwb;
  double& mdl_cehre = madMela::params_r_.mdl_cehre;
  double& mdl_cuhre = madMela::params_r_.mdl_cuhre;
  double& mdl_cdhre = madMela::params_r_.mdl_cdhre;
  double& mdl_cewre = madMela::params_r_.mdl_cewre;
  double& mdl_cebre = madMela::params_r_.mdl_cebre;
  double& mdl_cugre = madMela::params_r_.mdl_cugre;
  double& mdl_cuwre = madMela::params_r_.mdl_cuwre;
  double& mdl_cubre = madMela::params_r_.mdl_cubre;
  double& mdl_cdgre = madMela::params_r_.mdl_cdgre;
  double& mdl_cdwre = madMela::params_r_.mdl_cdwre;
  double& mdl_cdbre = madMela::params_r_.mdl_cdbre;
  double& mdl_chl1 = madMela::params_r_.mdl_chl1;
  double& mdl_chl3 = madMela::params_r_.mdl_chl3;
  double& mdl_che = madMela::params_r_.mdl_che;
  double& mdl_chq1 = madMela::params_r_.mdl_chq1;
  double& mdl_chq3 = madMela::params_r_.mdl_chq3;
  double& mdl_chu = madMela::params_r_.mdl_chu;
  double& mdl_chd = madMela::params_r_.mdl_chd;
  double& mdl_chudre = madMela::params_r_.mdl_chudre;
  double& mdl_cll = madMela::params_r_.mdl_cll;
  double& mdl_cll1 = madMela::params_r_.mdl_cll1;
  double& mdl_cqq1 = madMela::params_r_.mdl_cqq1;
  double& mdl_cqq11 = madMela::params_r_.mdl_cqq11;
  double& mdl_cqq3 = madMela::params_r_.mdl_cqq3;
  double& mdl_cqq31 = madMela::params_r_.mdl_cqq31;
  double& mdl_clq1 = madMela::params_r_.mdl_clq1;
  double& mdl_clq3 = madMela::params_r_.mdl_clq3;
  double& mdl_cee = madMela::params_r_.mdl_cee;
  double& mdl_cuu = madMela::params_r_.mdl_cuu;
  double& mdl_cuu1 = madMela::params_r_.mdl_cuu1;
  double& mdl_cdd = madMela::params_r_.mdl_cdd;
  double& mdl_cdd1 = madMela::params_r_.mdl_cdd1;
  double& mdl_ceu = madMela::params_r_.mdl_ceu;
  double& mdl_ced = madMela::params_r_.mdl_ced;
  double& mdl_cud1 = madMela::params_r_.mdl_cud1;
  double& mdl_cud8 = madMela::params_r_.mdl_cud8;
  double& mdl_cle = madMela::params_r_.mdl_cle;
  double& mdl_clu = madMela::params_r_.mdl_clu;
  double& mdl_cld = madMela::params_r_.mdl_cld;
  double& mdl_cqe = madMela::params_r_.mdl_cqe;
  double& mdl_cqu1 = madMela::params_r_.mdl_cqu1;
  double& mdl_cqu8 = madMela::params_r_.mdl_cqu8;
  double& mdl_cqd1 = madMela::params_r_.mdl_cqd1;
  double& mdl_cqd8 = madMela::params_r_.mdl_cqd8;
  double& mdl_cledqre = madMela::params_r_.mdl_cledqre;
  double& mdl_cquqd1re = madMela::params_r_.mdl_cquqd1re;
  double& mdl_cquqd11re = madMela::params_r_.mdl_cquqd11re;
  double& mdl_cquqd8re = madMela::params_r_.mdl_cquqd8re;
  double& mdl_cquqd81re = madMela::params_r_.mdl_cquqd81re;
  double& mdl_clequ1re = madMela::params_r_.mdl_clequ1re;
  double& mdl_clequ3re = madMela::params_r_.mdl_clequ3re;
  double& mdl_cgtil = madMela::params_r_.mdl_cgtil;
  double& mdl_cwtil = madMela::params_r_.mdl_cwtil;
  double& mdl_chgtil = madMela::params_r_.mdl_chgtil;
  double& mdl_chwtil = madMela::params_r_.mdl_chwtil;
  double& mdl_chbtil = madMela::params_r_.mdl_chbtil;
  double& mdl_chwbtil = madMela::params_r_.mdl_chwbtil;
  double& mdl_cewim = madMela::params_r_.mdl_cewim;
  double& mdl_cebim = madMela::params_r_.mdl_cebim;
  double& mdl_cugim = madMela::params_r_.mdl_cugim;
  double& mdl_cuwim = madMela::params_r_.mdl_cuwim;
  double& mdl_cubim = madMela::params_r_.mdl_cubim;
  double& mdl_cdgim = madMela::params_r_.mdl_cdgim;
  double& mdl_cdwim = madMela::params_r_.mdl_cdwim;
  double& mdl_cdbim = madMela::params_r_.mdl_cdbim;
  double& mdl_chudim = madMela::params_r_.mdl_chudim;
  double& mdl_cehim = madMela::params_r_.mdl_cehim;
  double& mdl_cuhim = madMela::params_r_.mdl_cuhim;
  double& mdl_cdhim = madMela::params_r_.mdl_cdhim;
  double& mdl_cledqim = madMela::params_r_.mdl_cledqim;
  double& mdl_cquqd1im = madMela::params_r_.mdl_cquqd1im;
  double& mdl_cquqd8im = madMela::params_r_.mdl_cquqd8im;
  double& mdl_cquqd11im = madMela::params_r_.mdl_cquqd11im;
  double& mdl_cquqd81im = madMela::params_r_.mdl_cquqd81im;
  double& mdl_clequ1im = madMela::params_r_.mdl_clequ1im;
  double& mdl_clequ3im = madMela::params_r_.mdl_clequ3im;

  //madMELA CKM Values
  double& mdl_ckmlambda = madMela::params_r_.mdl_ckmlambda;
  double& mdl_ckma = madMela::params_r_.mdl_ckma;
  double& mdl_ckmrho = madMela::params_r_.mdl_ckmrho;
  double& mdl_ckmeta = madMela::params_r_.mdl_ckmeta;


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

