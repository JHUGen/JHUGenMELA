#ifndef TCOUPLINGSBASE_HH
#define TCOUPLINGSBASE_HH



/**
 * @file TCouplingsBase.hh
 * This file defines the indices for different couplings in MELA
 * 
*/

//---------------------------------
// Coupling array sizes
//---------------------------------
namespace{
  /**
   * @brief This is the enumeration for couplings between the Higgs and the kappa formulation of quarks
   * @sa Mela::selfDHqqcoupl
   * @remark Table listed as enumeration name on the left
   * @remark the corresponding pyMela/colloquial coupling name, as well as the value of the enum, on the right
  */
  enum CouplingIndex_HQQ{
    gHIGGS_KAPPA, /*!< `kappa` (Value=0) */
    gHIGGS_KAPPA_TILDE, /*!< `kappa` (Value=1) */

    SIZE_HQQ /*!< The size of the array **(Value=2)** */
  };

  /**
   * @brief This is the enumerator for the couplings between the Higgs and gluons
   * @sa Mela::selfDHggcoupl
   * @remark Table listed as enumeration name on the left
   * @remark the corresponding pyMela/colloquial coupling name, as well as the value of the enum, on the right
  */
  enum CouplingIndex_HGG{
    gHIGGS_GG_2, /*!< `ghg2` (Value=0) */
    gHIGGS_GG_3, /*!< `ghg3` (Value=1) */
    gHIGGS_GG_4, /*!< `ghg4` (Value=2) */

    SIZE_HGG /*!<  The size of the array **(Value=3) */
  };

  /**
   * @brief This is the enumerator for the couplings between the Higgs and the vector bosons (Z/Z' & W/W')
   * @sa Mela::selfDHzzcoupl, Mela::selfDHzzpcoupl, Mela::selfDHzpzpcoupl
   * @sa Mela::selfDHwwcoupl, Mela::selfDHwwpcoupl, Mela::selfDHwpwpcoupl
   * @attention Coupling names are given as though they are couplings between the z and H
   * @attention See the full coupling array list in Mela.h for a thorough definition of what couplings use these values
   * @remark Table listed as enumeration name on the left
   * @remark the corresponding pyMela/colloquial coupling name, as well as the value of the enum, on the right
  */
  enum CouplingIndex_HVV{
    gHIGGS_VV_1, /*!< `gh{z,w}1` (Value=0) */
    gHIGGS_VV_2, /*!< `gh{z,w}2` (Value=1) */
    gHIGGS_VV_3, /*!< `gh{z,w}3` (Value=2) */
    gHIGGS_VV_4, /*!< `gh{z,w}4` (Value=3) */

    gHIGGS_ZA_2, /*!< `ghzgs2` (Value=4) */
    gHIGGS_ZA_3, /*!< `ghzgs3` (Value=5) */
    gHIGGS_ZA_4, /*!< `ghzgs4` (Value=6) */

    gHIGGS_AA_2, /*!< `ghgsgs2` (Value=7) */
    gHIGGS_AA_3, /*!< `ghgsgs3` (Value=8) */
    gHIGGS_AA_4, /*!< `ghgsgs4` (Value=9) */

    gHIGGS_VV_1_PRIME,  /*!< `ghz1_prime` (Value=10) */
    gHIGGS_VV_1_PRIME2, /*!< `ghz1_prime2` (Value=11) */
    gHIGGS_VV_1_PRIME3, /*!< `ghz1_prime3` (Value=12) */
    gHIGGS_VV_1_PRIME4, /*!< `ghz1_prime4` (Value=13) */
    gHIGGS_VV_1_PRIME5, /*!< `ghz1_prime5` (Value=14) */

    gHIGGS_VV_2_PRIME,  /*!< `ghz2_prime` (Value=15) */
    gHIGGS_VV_2_PRIME2, /*!< `ghz2_prime2` (Value=16) */
    gHIGGS_VV_2_PRIME3, /*!< `ghz2_prime3` (Value=17) */
    gHIGGS_VV_2_PRIME4, /*!< `ghz2_prime4` (Value=18) */
    gHIGGS_VV_2_PRIME5, /*!< `ghz2_prime5` (Value=19) */

    gHIGGS_VV_3_PRIME,  /*!< `ghz3_prime` (Value=20) */
    gHIGGS_VV_3_PRIME2, /*!< `ghz3_prime2` (Value=21) */
    gHIGGS_VV_3_PRIME3, /*!< `ghz3_prime3` (Value=22) */
    gHIGGS_VV_3_PRIME4, /*!< `ghz3_prime4` (Value=23) */
    gHIGGS_VV_3_PRIME5, /*!< `ghz3_prime5` (Value=24) */

    gHIGGS_VV_4_PRIME,  /*!< `ghz4_prime`  (Value=25) */
    gHIGGS_VV_4_PRIME2, /*!< `ghz4_prime2` (Value=26) */
    gHIGGS_VV_4_PRIME3, /*!< `ghz4_prime3` (Value=27) */
    gHIGGS_VV_4_PRIME4, /*!< `ghz4_prime4` (Value=28) */
    gHIGGS_VV_4_PRIME5, /*!< `ghz4_prime5` (Value=29) */

    gHIGGS_ZA_1_PRIME2, /*!< `ghzgs1_prime2` (Value=30) */

    gHIGGS_VV_1_PRIME6, /*!< `ghz1_prime6` (Value=31) */
    gHIGGS_VV_1_PRIME7, /*!< `ghz1_prime6` (Value=32) */
    gHIGGS_VV_2_PRIME6, /*!< `ghz2_prime6` (Value=33) */
    gHIGGS_VV_2_PRIME7, /*!< `ghz2_prime6` (Value=34) */
    gHIGGS_VV_3_PRIME6, /*!< `ghz3_prime6` (Value=35) */
    gHIGGS_VV_3_PRIME7, /*!< `ghz3_prime6` (Value=36) */
    gHIGGS_VV_4_PRIME6, /*!< `ghz4_prime6` (Value=37) */
    gHIGGS_VV_4_PRIME7, /*!< `ghz4_prime6` (Value=38) */

    SIZE_HVV /*!<  The size of the array   (Value=39) */
  };

  /**
   * @brief This is the enumeration for couplings between the Higgs and the lambda Z/W couplings
   * @sa Mela::selfDHzzLambda_qsq, selfDHwwLambda_qsq
   * @remark Table listed as enumeration name on the left
   * @remark the corresponding pyMela/colloquial coupling name, as well as the value of the enum, on the right
  */
  enum CouplingIndex_LAMBDAQSQ{
    LambdaHIGGS_QSQ_VV_1 = 0,
    LambdaHIGGS_QSQ_VV_2 = 1,
    LambdaHIGGS_QSQ_VV_3 = 2,
    LambdaHIGGS_QSQ_VV_4 = 3,

    SIZE_HVV_LAMBDAQSQ = 4
  };
  enum CouplingIndex_HVV_CQSQ{
    cLambdaHIGGS_VV_QSQ1 = 0,
    cLambdaHIGGS_VV_QSQ2 = 1,
    cLambdaHIGGS_VV_QSQ12 = 2,

    SIZE_HVV_CQSQ = 3
  };
  enum CouplingIndex_Vpff{
    gHIGGS_Vp_El_left,
    gHIGGS_Vp_El_right,
    gHIGGS_Vp_Mu_left,
    gHIGGS_Vp_Mu_right,
    gHIGGS_Vp_Ta_left,
    gHIGGS_Vp_Ta_right,
    gHIGGS_Vp_NuE_left,
    gHIGGS_Vp_NuE_right,

    gHIGGS_Vp_Dn_left,
    gHIGGS_Vp_Dn_right,
    gHIGGS_Vp_Up_left,
    gHIGGS_Vp_Up_right,
    gHIGGS_Vp_Str_left,
    gHIGGS_Vp_Str_right,
    gHIGGS_Vp_Chm_left,
    gHIGGS_Vp_Chm_right,
    gHIGGS_Vp_Bot_left,
    gHIGGS_Vp_Bot_right,
    gHIGGS_Vp_Top_left,
    gHIGGS_Vp_Top_right,

    SIZE_Vpff
  };
  enum CouplingIndex_ZQQ{
    gZPRIME_QQ_LEFT,
    gZPRIME_QQ_RIGHT,

    SIZE_ZQQ
  };
  enum CouplingIndex_ZVV{
    gZPRIME_VV_1,
    gZPRIME_VV_2,

    SIZE_ZVV
  };
  enum CouplingIndex_GQQ{
    gGRAVITON_QQ_LEFT,
    gGRAVITON_QQ_RIGHT,

    SIZE_GQQ
  };
  enum CouplingIndex_GGG{
    gGRAVITON_GG_1,
    gGRAVITON_GG_2,
    gGRAVITON_GG_3,
    gGRAVITON_GG_4,
    gGRAVITON_GG_5,

    SIZE_GGG
  };
  enum CouplingIndex_GVV{
    gGRAVITON_VV_1,
    gGRAVITON_VV_2,
    gGRAVITON_VV_3,
    gGRAVITON_VV_4,
    gGRAVITON_VV_5,
    gGRAVITON_VV_6,
    gGRAVITON_VV_7,
    gGRAVITON_VV_8,
    gGRAVITON_VV_9,
    gGRAVITON_VV_10,

    gGRAVITON_ZA_1,
    gGRAVITON_ZA_2,
    gGRAVITON_ZA_3,
    gGRAVITON_ZA_4,
    gGRAVITON_ZA_8,

    gGRAVITON_AA_1,
    gGRAVITON_AA_2,
    gGRAVITON_AA_3,
    gGRAVITON_AA_4,
    gGRAVITON_AA_8,

    SIZE_GVV
  };
  enum CouplingIndex_ATQGC{
    gATQGC_dVA,
    gATQGC_dPA,
    gATQGC_dMA,
    gATQGC_dFourA,

    gATQGC_dVZ,
    gATQGC_dPZ,
    gATQGC_dMZ,
    gATQGC_dFourZ,

    gATQGC_dAAWpWm,
    gATQGC_dZAWpWm,
    gATQGC_dZZWpWm,

    SIZE_ATQGC
  };
  enum CouplingIndex_AZff{
    gAZff_ZllRH,
    gAZff_ZllLH,
    gAZff_ZuuRH,
    gAZff_ZuuLH,
    gAZff_ZddRH,
    gAZff_ZddLH,
    gAZff_ZnunuRH,
    gAZff_ZnunuLH,

    gAZff_uZRH,
    gAZff_uZLH,
    gAZff_dZRH,
    gAZff_dZLH,

    SIZE_AZff
  };
  enum CouplingIndex_SMEFT{
    gMDL_ch,
    gMDL_chbox,
    gMDL_chdd,
    gMDL_chg,
    gMDL_chw,
    gMDL_chb,
    gMDL_chwb,
    gMDL_cehre,
    gMDL_cuhre,
    gMDL_cdhre,
    gMDL_cewre,
    gMDL_cebre,
    gMDL_cugre,
    gMDL_cuwre,
    gMDL_cubre,
    gMDL_cdgre,
    gMDL_cdwre,
    gMDL_cdbre,
    gMDL_chl1,
    gMDL_chl3,
    gMDL_che,
    gMDL_chq1,
    gMDL_chq3,
    gMDL_chu,
    gMDL_chd,
    gMDL_chudre,
    gMDL_cll,
    gMDL_cll1,
    gMDL_cqq1,
    gMDL_cqq11,
    gMDL_cqq3,
    gMDL_cqq31,
    gMDL_clq1,
    gMDL_clq3,
    gMDL_cee,
    gMDL_cuu,
    gMDL_cuu1,
    gMDL_cdd,
    gMDL_cdd1,
    gMDL_ceu,
    gMDL_ced,
    gMDL_cud1,
    gMDL_cud8,
    gMDL_cle,
    gMDL_clu,
    gMDL_cld,
    gMDL_cqe,
    gMDL_cqu1,
    gMDL_cqu8,
    gMDL_cqd1,
    gMDL_cqd8,
    gMDL_cledqre,
    gMDL_cquqd1re,
    gMDL_cquqd11re,
    gMDL_cquqd8re,
    gMDL_cquqd81re,
    gMDL_clequ1re,
    gMDL_clequ3re,
    gMDL_cgtil,
    gMDL_cwtil,
    gMDL_chgtil,
    gMDL_chwtil,
    gMDL_chbtil,
    gMDL_chwbtil,
    gMDL_cewim,
    gMDL_cebim,
    gMDL_cugim,
    gMDL_cuwim,
    gMDL_cubim,
    gMDL_cdgim,
    gMDL_cdwim,
    gMDL_cdbim,
    gMDL_chudim,
    gMDL_cehim,
    gMDL_cuhim,
    gMDL_cdhim,
    gMDL_cledqim,
    gMDL_cquqd1im,
    gMDL_cquqd8im,
    gMDL_cquqd11im,
    gMDL_cquqd81im,
    gMDL_clequ1im,
    gMDL_clequ3im,

    SIZE_SMEFT
  };
}


#endif
