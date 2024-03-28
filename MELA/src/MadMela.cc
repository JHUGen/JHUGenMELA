#include <iostream>
#include <stdexcept>
#include "MadMela.h"
#include "MELAStreamHelpers.hh"
#include "TLorentzVector.h"
#include "TUtil.hh"
#include "TVar.hh"
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;

madMela::madMelaInput* madMela::madMelaCandidate = nullptr;
map<int, double*> madMela::mass_map;
map<int, double*> madMela::width_map;
map<pair<int, int>, pair<complex<double>*, complex<double>*> > madMela::CKM_map;
TVar::VerbosityLevel madMela::myVerbosity_ = TVar::SILENT;

void madMela::setDefaultMadgraphValues(){
    //Setting the CKM Matrix for MADGRAPH
    madMela::params_c_.mdl_conjg__ckm1x1 =  complex<double>( 0.97434887500000000     , -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm1x2 =  complex<double>( 0.22650000000000001     , -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm1x3 =  complex<double>(  1.2943473026287501E-003,  3.2771772130387503E-003);
    madMela::params_c_.mdl_conjg__ckm2x1 =  complex<double>(-0.22650000000000001     , -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm2x2 =  complex<double>( 0.97434887500000000     , -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm2x3 =  complex<double>(  4.0528777500000002E-002, -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm3x1 =  complex<double>(  7.8854208011212516E-003,  3.2771772130387503E-003);
    madMela::params_c_.mdl_conjg__ckm3x2 =  complex<double>( -4.0528777500000002E-002, -0.0000000000000000     );
    madMela::params_c_.mdl_conjg__ckm3x3 =  complex<double>(  1.0000000000000000     , -0.0000000000000000     );

    madMela::params_c_.mdl_ckm1x1 =  complex<double>( 0.97434887500000000     ,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm1x2 =  complex<double>( 0.22650000000000001     ,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm1x3 =  complex<double>(  1.2943473026287501E-003, -3.2771772130387503E-003);
    madMela::params_c_.mdl_ckm2x1 =  complex<double>(-0.22650000000000001     ,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm2x2 =  complex<double>( 0.97434887500000000     ,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm2x3 =  complex<double>(  4.0528777500000002E-002,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm3x1 =  complex<double>(  7.8854208011212516E-003, -3.2771772130387503E-003);
    madMela::params_c_.mdl_ckm3x2 =  complex<double>( -4.0528777500000002E-002,  0.0000000000000000     );
    madMela::params_c_.mdl_ckm3x3 =  complex<double>(  1.0000000000000000     ,  0.0000000000000000     );

    madMela::params_r_.mdl_ckmlambda =   0.22650000000000001;
    madMela::params_r_.mdl_ckmlambda__exp__2 =    5.1302250000000001E-002;
    madMela::params_r_.mdl_ckmlambda__exp__3 =    1.1619959625000001E-002;
    madMela::params_r_.mdl_ckma =   0.79000000000000004;
    madMela::params_r_.mdl_ckmrho =   0.14099999999999999;
    madMela::params_r_.mdl_ckmeta =   0.35699999999999998;

    //Nominal Masses
    madMela::mad_masses_.mdl_ms = 0.093;
    madMela::mad_masses_.mdl_mu = 0.00216;
    madMela::mad_masses_.mdl_mt = madMela::mad_masses_.mdl_mt1 = 172.9;

    madMela::mad_masses_.mdl_mw = madMela::mad_masses_.mdl_mw1 = madMela::params_r_.mdl_mwsm =    79.831336335943078;
    madMela::params_r_.mdl_mwsm__exp__2 =    6373.0422611824652;
    madMela::params_r_.mdl_mwsm__exp__4 =    40615667.662817709;
    madMela::params_r_.mdl_mwsm__exp__6 =    258845366481.27933;


    madMela::mad_masses_.mdl_mh = madMela::mad_masses_.mdl_mh1 = 125;
    madMela::mad_masses_.mdl_mta = 1.777;
    madMela::mad_masses_.mdl_mc = 1.27;
    madMela::mad_masses_.mdl_mb = 4.18;
    madMela::mad_masses_.mdl_me = 0.000511;
    madMela::mad_masses_.mdl_md = 0.00467;
    madMela::mad_masses_.mdl_mmu = 0.10566;

    madMela::mad_masses_.mdl_mz = 91.1876;
    madMela::mad_masses_.mdl_mz1 = 91.1876;
    madMela::params_r_.mdl_mz__exp__2 =    8315.1783937600012;
    madMela::params_r_.mdl_mz__exp__3 =    758241.16129882948;
    madMela::params_r_.mdl_mz__exp__4 =    69142191.720053151;
    madMela::params_r_.mdl_mz__exp__6 =    574929658687.79749;

    madMela::params_r_.mdl_mt__exp__2 =    29846.017599999996;

    madMela::params_r_.mdl_mh__exp__4 =    244844509.73956564;
    madMela::params_r_.mdl_mh__exp__6 =    3831206449390.3818;

    madMela::params_r_.mdl_mb__exp__2 =    17.472399999999997;

    //Nominal widths
    madMela::widths_.mdl_wz = madMela::widths_.mdl_wz1 =   2.4952000000000001;
    madMela::widths_.mdl_ww = madMela::widths_.mdl_ww1 =   2.0850000000000000;
    madMela::widths_.mdl_wt = madMela::widths_.mdl_wt1 =   1.3300000000000001;
    madMela::widths_.mdl_wh = madMela::widths_.mdl_wh1 =   4.0699999999999998E-003;

    //Yukawa Couplings
    madMela::params_r_.mdl_ymdo =    4.6699999999999997E-003;
    madMela::params_r_.mdl_ymup =    2.1600000000000000E-003;
    madMela::params_r_.mdl_yms =    9.2999999999999999E-002;
    madMela::params_r_.mdl_ymc =    1.2700000000000000;
    madMela::params_r_.mdl_ymb =    4.1799999999999997;
    madMela::params_r_.mdl_ymt =    172.75999999999999;
    madMela::params_r_.mdl_yme =    5.1099999999999995E-004;
    madMela::params_r_.mdl_ymm =   0.10566000000000000;
    madMela::params_r_.mdl_ymtau =    1.7769999999999999;

    madMela::params_r_.mdl_yb =    2.4008698866559607E-002;
    madMela::params_r_.mdl_yc =    7.2945089857728964E-003;
    madMela::params_r_.mdl_ydo =    2.6823115719338128E-005;
    madMela::params_r_.mdl_ye =    2.9350347178975982E-006;
    madMela::params_r_.mdl_ym =    6.0688017278485377E-004;
    madMela::params_r_.mdl_ys =    5.3416483124163721E-004;
    madMela::params_r_.mdl_yt =   0.99228297037962632;
    madMela::params_r_.mdl_ytau =    1.0206568872219240E-002;
    madMela::params_r_.mdl_yup =    1.2406408983676737E-005;

    //Stored values that are derived
    madMela::params_c_.mdl_complexi =  complex<double>(  0.0000000000000000     ,  1.0000000000000000     );
    madMela::params_r_.mdl_nb__10__exp___m_40 =    9.9999999999999993E-041;
    madMela::params_r_.mdl_lambdasmeft__exp__2 =    1000000.0000000000;
    madMela::params_r_.mdl_ee__exp__2 =    9.8213135554166664E-002;

    madMela::params_r_.mdl_propcorr =    0.0000000000000000;
    madMela::params_r_.mdl_propcorr__exp__2 =    0.0000000000000000;
    madMela::params_r_.mdl_propcorr__exp__3 =    0.0000000000000000;
    madMela::params_r_.mdl_propcorr__exp__4 =    0.0000000000000000;
    madMela::params_r_.mdl_ee__exp__3 =    3.0778990021400259E-002;


    //A variety of standard values
    madMela::params_r_.mdl_lambdasmeft = 1000;
    madMela::params_r_.mdl_aew = 7.8155529999999994E-003;
    madMela::params_r_.mdl_ee =   0.31338975023788934;

    madMela::params_r_.mdl_gf = 1.1663789999999999E-005;
    madMela::params_r_.mdl_sqrt__gf =    3.4152291284773265E-003;

    madMela::params_r_.as = 0.11790000000000000;

    madMela::params_r_.mdl_vev =    246.21961912951551;
    madMela::params_r_.mdl_vevhat =    246.21961912951551;
    madMela::params_r_.mdl_vevhat__exp__2 =    60624.100844283676;
    madMela::params_r_.mdl_vevhat__exp__3 =    14926843.019948866;
    madMela::params_r_.mdl_vevt =    246.21961912951551;

    madMela::params_r_.mdl_lam =   0.12905352724481209;

    madMela::params_r_.mdl_cth =   0.87546263237483035;
    madMela::params_r_.mdl_cth__exp__2 =   0.76643482068466739;
    madMela::params_r_.mdl_cth__exp__3 =   0.67098504566033002;
    madMela::params_r_.mdl_cth__exp__4 =   0.58742233435793823;
    madMela::params_r_.mdl_cth__exp__5 =   0.51426630315276844;

    madMela::params_r_.mdl_sth2 =   0.23356517931533255;
    madMela::params_r_.mdl_sqrt__sth2 =   0.48328581534670823;

    madMela::params_r_.mdl_sth =   0.48328581534670823;
    madMela::params_r_.mdl_sth__exp__2 =   0.23356517931533255;
    madMela::params_r_.mdl_sth__exp__3 =   0.11287873812201060;
    madMela::params_r_.mdl_sth__exp__4 =    5.4552692988603449E-002;
    madMela::params_r_.mdl_sth__exp__6 =    1.2741609520017450E-002;
    madMela::params_r_.mdl_sth__exp__5 =    2.6364542710355869E-002;

    //A series of shifts
    madMela::params_r_.mdl_dgf =    0.0000000000000000;
    madMela::params_r_.mdl_dkh =    0.0000000000000000;
    madMela::params_r_.mdl_dmz2 =   0.0000000000000000;
    madMela::params_r_.mdl_dmh2 =   0.0000000000000000;
    madMela::params_r_.mdl_dg1 =    0.0000000000000000;
    madMela::params_r_.mdl_dgw =    0.0000000000000000;
    madMela::params_r_.mdl_dwz =    0.0000000000000000;
    madMela::params_r_.mdl_dmw =    0.0000000000000000;
    madMela::params_r_.mdl_dwt =    0.0000000000000000;
    madMela::params_r_.mdl_dww =    0.0000000000000000;
    madMela::params_r_.mdl_dwhc =   0.0000000000000000;
    madMela::params_r_.mdl_dwhb =   0.0000000000000000;
    madMela::params_r_.mdl_dwhta =  0.0000000000000000;

    madMela::params_r_.mdl_g1 =   0.35797044745104684;
    madMela::params_r_.mdl_gw =   0.64845633843621953;

    madMela::params_r_.mdl_barlam =   0.12905352724481209;
    madMela::params_r_.mdl_ghaa =   -2.0035938810065040E-003;
    madMela::params_r_.mdl_ghza =   -6.8918046831024345E-003;
    madMela::params_r_.mdl_g1sh =   0.35797044745104684;
    madMela::params_r_.mdl_gwsh =   0.64845633843621953;


    madMela::params_r_.mdl_sqrt__as =   0.34336569426778790;

    madMela::params_r_.mdl_ghgg1 =    3.1273946317557436E-003;
    madMela::params_r_.mdl_ghgg2 =   -3.6486270703817007E-004;
    madMela::params_r_.mdl_ghgg3 =    7.6133275785037098E-004;
    madMela::params_r_.mdl_ghgg4 =    1.0424648772519146E-004;
    madMela::params_r_.mdl_ghgg5 =    1.8764367790534461E-003;

    madMela::params_r_.mdl_g__exp__2 =    1.4815750954329465;
    madMela::params_c_.mdl_g__exp__3 =  complex<double>(  1.8033727530176447     ,  0.0000000000000000     );
    madMela::params_r_.mdl_dwh =    0.0000000000000000;

    //Complex valued couplings
    madMela::params_c_.mdl_ceh = complex<double>(0, 0);
    madMela::params_c_.mdl_cuh = complex<double>(0, 0);
    madMela::params_c_.mdl_cdh = complex<double>(0, 0);
    madMela::params_c_.mdl_cew = complex<double>(0, 0);
    madMela::params_c_.mdl_ceb = complex<double>(0, 0);
    madMela::params_c_.mdl_cug = complex<double>(0, 0);
    madMela::params_c_.mdl_cuw = complex<double>(0, 0);
    madMela::params_c_.mdl_cub = complex<double>(0, 0);
    madMela::params_c_.mdl_cdg = complex<double>(0, 0);
    madMela::params_c_.mdl_cdw = complex<double>(0, 0);
    madMela::params_c_.mdl_cdb = complex<double>(0, 0);
    madMela::params_c_.mdl_chud = complex<double>(0, 0);
    madMela::params_c_.mdl_cledq = complex<double>(0, 0);
    madMela::params_c_.mdl_cquqd1 = complex<double>(0, 0);
    madMela::params_c_.mdl_cquqd11 = complex<double>(0, 0);
    madMela::params_c_.mdl_cquqd8 = complex<double>(0, 0);
    madMela::params_c_.mdl_cquqd81 = complex<double>(0, 0);
    madMela::params_c_.mdl_clequ1 = complex<double>(0, 0);
    madMela::params_c_.mdl_clequ3 = complex<double>(0, 0);

    //Real valued couplings
    madMela::params_r_.mdl_ch = 0.;
    madMela::params_r_.mdl_chbox = 0.;
    madMela::params_r_.mdl_chdd = 0.;
    madMela::params_r_.mdl_chg = 0.;
    madMela::params_r_.mdl_chw = 0.;
    madMela::params_r_.mdl_chb = 0.;
    madMela::params_r_.mdl_chwb = 0.;
    madMela::params_r_.mdl_cehre = 0.;
    madMela::params_r_.mdl_cuhre = 0.;
    madMela::params_r_.mdl_cdhre = 0.;
    madMela::params_r_.mdl_cewre = 0.;
    madMela::params_r_.mdl_cebre = 0.;
    madMela::params_r_.mdl_cugre = 0.;
    madMela::params_r_.mdl_cuwre = 0.;
    madMela::params_r_.mdl_cubre = 0.;
    madMela::params_r_.mdl_cdgre = 0.;
    madMela::params_r_.mdl_cdwre = 0.;
    madMela::params_r_.mdl_cdbre = 0.;
    madMela::params_r_.mdl_chl1 = 0.;
    madMela::params_r_.mdl_chl3 = 0.;
    madMela::params_r_.mdl_che = 0.;
    madMela::params_r_.mdl_chq1 = 0.;
    madMela::params_r_.mdl_chq3 = 0.;
    madMela::params_r_.mdl_chu = 0.;
    madMela::params_r_.mdl_chd = 0.;
    madMela::params_r_.mdl_chudre = 0.;
    madMela::params_r_.mdl_cll = 0.;
    madMela::params_r_.mdl_cll1 = 0.;
    madMela::params_r_.mdl_cqq1 = 0.;
    madMela::params_r_.mdl_cqq11 = 0.;
    madMela::params_r_.mdl_cqq3 = 0.;
    madMela::params_r_.mdl_cqq31 = 0.;
    madMela::params_r_.mdl_clq1 = 0.;
    madMela::params_r_.mdl_clq3 = 0.;
    madMela::params_r_.mdl_cee = 0.;
    madMela::params_r_.mdl_cuu = 0.;
    madMela::params_r_.mdl_cuu1 = 0.;
    madMela::params_r_.mdl_cdd = 0.;
    madMela::params_r_.mdl_cdd1 = 0.;
    madMela::params_r_.mdl_ceu = 0.;
    madMela::params_r_.mdl_ced = 0.;
    madMela::params_r_.mdl_cud1 = 0.;
    madMela::params_r_.mdl_cud8 = 0.;
    madMela::params_r_.mdl_cle = 0.;
    madMela::params_r_.mdl_clu = 0.;
    madMela::params_r_.mdl_cld = 0.;
    madMela::params_r_.mdl_cqe = 0.;
    madMela::params_r_.mdl_cqu1 = 0.;
    madMela::params_r_.mdl_cqu8 = 0.;
    madMela::params_r_.mdl_cqd1 = 0.;
    madMela::params_r_.mdl_cqd8 = 0.;
    madMela::params_r_.mdl_cledqre = 0.;
    madMela::params_r_.mdl_cquqd1re = 0.;
    madMela::params_r_.mdl_cquqd11re = 0.;
    madMela::params_r_.mdl_cquqd8re = 0.;
    madMela::params_r_.mdl_cquqd81re = 0.;
    madMela::params_r_.mdl_clequ1re = 0.;
    madMela::params_r_.mdl_clequ3re = 0.;
    madMela::params_r_.mdl_cgtil = 0.;
    madMela::params_r_.mdl_cwtil = 0.;
    madMela::params_r_.mdl_chgtil = 0.;
    madMela::params_r_.mdl_chwtil = 0.;
    madMela::params_r_.mdl_chbtil = 0.;
    madMela::params_r_.mdl_chwbtil = 0.;
    madMela::params_r_.mdl_cewim = 0.;
    madMela::params_r_.mdl_cebim = 0.;
    madMela::params_r_.mdl_cugim = 0.;
    madMela::params_r_.mdl_cuwim = 0.;
    madMela::params_r_.mdl_cubim = 0.;
    madMela::params_r_.mdl_cdgim = 0.;
    madMela::params_r_.mdl_cdwim = 0.;
    madMela::params_r_.mdl_cdbim = 0.;
    madMela::params_r_.mdl_chudim = 0.;
    madMela::params_r_.mdl_cehim = 0.;
    madMela::params_r_.mdl_cuhim = 0.;
    madMela::params_r_.mdl_cdhim = 0.;
    madMela::params_r_.mdl_cledqim = 0.;
    madMela::params_r_.mdl_cquqd1im = 0.;
    madMela::params_r_.mdl_cquqd8im = 0.;
    madMela::params_r_.mdl_cquqd11im = 0.;
    madMela::params_r_.mdl_cquqd81im = 0.;
    madMela::params_r_.mdl_clequ1im = 0.;
    madMela::params_r_.mdl_clequ3im = 0.;

    delete madMela::madMelaCandidate;
    madMela::madMelaCandidate = nullptr;
}

void madMela::initialize_madMELA(){
    madMela::setDefaultMadgraphValues();
    madMela::mass_map[1] =  &madMela::mad_masses_.mdl_md;
    madMela::mass_map[2] =  &madMela::mad_masses_.mdl_mu;
    madMela::mass_map[3] =  &madMela::mad_masses_.mdl_ms;
    madMela::mass_map[4] =  &madMela::mad_masses_.mdl_mc;
    madMela::mass_map[5] =  &madMela::mad_masses_.mdl_mb;
    madMela::mass_map[6] =  &madMela::mad_masses_.mdl_mt;

    madMela::mass_map[11] = &madMela::mad_masses_.mdl_me;
    madMela::mass_map[13] = &madMela::mad_masses_.mdl_mmu;
    madMela::mass_map[15] = &madMela::mad_masses_.mdl_mta;

    madMela::mass_map[23] = &madMela::mad_masses_.mdl_mz;
    madMela::mass_map[24] = &madMela::mad_masses_.mdl_mw;
    madMela::mass_map[25] = &madMela::mad_masses_.mdl_mh;

    madMela::width_map[6] = &madMela::widths_.mdl_wt;
    madMela::width_map[23] = &madMela::widths_.mdl_wz;
    madMela::width_map[24] = &madMela::widths_.mdl_ww;
    madMela::width_map[25] = &madMela::widths_.mdl_wh;


    madMela::CKM_map[pair<int, int>(1,1)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm1x1, &madMela::params_c_.mdl_conjg__ckm1x1);
    madMela::CKM_map[pair<int, int>(1,2)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm1x2, &madMela::params_c_.mdl_conjg__ckm1x2);
    madMela::CKM_map[pair<int, int>(1,3)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm1x3, &madMela::params_c_.mdl_conjg__ckm1x3);
    madMela::CKM_map[pair<int, int>(2,1)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm2x1, &madMela::params_c_.mdl_conjg__ckm2x1);
    madMela::CKM_map[pair<int, int>(2,2)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm2x2, &madMela::params_c_.mdl_conjg__ckm2x2);
    madMela::CKM_map[pair<int, int>(2,3)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm2x3, &madMela::params_c_.mdl_conjg__ckm2x3);
    madMela::CKM_map[pair<int, int>(3,1)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm3x1, &madMela::params_c_.mdl_conjg__ckm3x1);
    madMela::CKM_map[pair<int, int>(3,2)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm3x2, &madMela::params_c_.mdl_conjg__ckm3x2);
    madMela::CKM_map[pair<int, int>(3,3)] = pair<complex<double>*, complex<double>*>(&madMela::params_c_.mdl_ckm3x3, &madMela::params_c_.mdl_conjg__ckm3x3);

}

void madMela::setInputEvent(SimpleParticleCollection_t* pDaughters, SimpleParticleCollection_t* pAssociated, SimpleParticleCollection_t* pMothers, bool isGen){
    if (!(pDaughters) || (pDaughters->size() == 0)){ throw invalid_argument("madMela::setInputEvent: No daughters!");}
    else if (pDaughters->size()>4){
        MELAerr << "madMela::setInputEvent: Daughter size " << pDaughters->size() << endl;
        throw invalid_argument("madMela::setInputEvent: Daughter size >4 is not supported!");
    }
    else if (!(pMothers) || pMothers->size()!=2){
        MELAerr << "TUtil::ConvertVectorFormat: Mothers momentum size (" << pMothers->size() << ") has to have had been 2!" << endl;
        throw invalid_argument("madMela::setInputEvent: Mother size != 4 is not supported!");
    } else if(!(pAssociated)){
        if(myVerbosity_>=TVar::DEBUG) MELAout << "No associated particles found. Instantiating empty vector." << endl;
        pAssociated = new SimpleParticleCollection_t();
    }

    const int nPDG = pDaughters->size() + pAssociated->size() + pMothers->size();
    vector<int> pdgs(nPDG);
    vector<vector<double>> p(10, vector<double>(4));

    int i = 0;
    for(SimpleParticle_t particle : *pMothers){
        pdgs[i] = particle.first;
        TLorentzVector motherVec = particle.second;
        p[i][0] = motherVec.E();
        p[i][1] = motherVec.Px();
        p[i][2] = motherVec.Py();
        p[i][3] = motherVec.Pz();
        i++;
    }

    bool previously_swapped = false;
    for(SimpleParticle_t particle : *pDaughters){
        bool swap_spaces = false;
        if((i == 2 || i == 4) && particle.first > 0){
            i++;
            swap_spaces = true;
            if(myVerbosity_ >= TVar::DEBUG) MELAout << "Swapping daughters at index " << i << " and " << i-1 << endl;
        }

        pdgs[i] = particle.first;
        TLorentzVector motherVec = particle.second;
        p[i][0] = motherVec.E();
        p[i][1] = motherVec.Px();
        p[i][2] = motherVec.Py();
        p[i][3] = motherVec.Pz();

        if(previously_swapped){
            i += 2;
            previously_swapped = false;
        }
        else if(swap_spaces){
            i--;
            previously_swapped = true;
        } else{
            i++;
        }
    }

    for(SimpleParticle_t particle : *pAssociated){
        pdgs[i] = particle.first;
        TLorentzVector motherVec = particle.second;
        p[i][0] = motherVec.E();
        p[i][1] = motherVec.Px();
        p[i][2] = motherVec.Py();
        p[i][3] = motherVec.Pz();
        i++;
    }

    madMela::madMelaCandidate = new madMela::madMelaInput(pdgs, p);
}

void madMela::computeP(double& prob, int nhel){
    //Calculate the effects of the set couplings
    madMela::coup_();

    vector<int> pdgs = madMela::madMelaCandidate->first;
    vector<vector<double>> p = madMela::madMelaCandidate->second;

    const int nPDG = pdgs.size();
    int pdgs_for_fortran[nPDG];
    double p_for_fortran[nPDG][4];

    copy(pdgs.begin(), pdgs.end(), pdgs_for_fortran);
    for(int n = 0; n < nPDG; n++){
        copy(p[n].begin(), p[n].end(), p_for_fortran[n]);
    }
    double scale2 = 0;

    if(madMela::myVerbosity_ >= TVar::DEBUG_VERBOSE){
        MELAout << "Input vectors to MADGRAPH function in order as id, px, py, pz, E:" << endl;
        for(int i = 0; i < nPDG; i++){
            MELAout << "id of " << pdgs[i] << " & vector of " << p[i][1] << ", " << p[i][2] << ", " << p[i][3] << ", " << p[i][0] << endl;
        }
    }
    smatrixhel_(pdgs_for_fortran, nPDG, p_for_fortran, madMela::params_r_.as, scale2, nhel, prob);
    if(madMela::myVerbosity_ >= TVar::DEBUG) MELAout << " smatrixhel returns prob of " << prob << endl;
    madMela::setDefaultMadgraphValues();//reset couplings
}

void madMela::set_mHiggs(double myHiggsMass, int index){
    if((index < 0) || (index > 1)){
        throw out_of_range("Higgs mass index can only be 1 or 0!");
    }
    else if(index == 0){
        madMela::mad_masses_.mdl_mh = myHiggsMass;
    } else{
        madMela::mad_masses_.mdl_mh1 = myHiggsMass;
    }
}

void madMela::set_wHiggs(double myHiggsWidth, int index){
    if((index < 0) || (index > 1)){
        throw out_of_range("Higgs mass index can only be 1 or 0!");
    }
    else if(index == 0){
        madMela::widths_.mdl_wh = myHiggsWidth;
    } else{
        madMela::widths_.mdl_wh1 = myHiggsWidth;
    }
}




// void madMela::setCKMMatEl(int i, int j, complex<double> matEl){
//     if( (i < 1) || (i > 3) || (j < 1) || (j > 3) ){
//         throw out_of_range("1 < {i,j} < 3!");
//     }

//     pair<complex<double>*, complex<double>*> matElAndConjugate= madMela::CKM_map(<int, int>(i, j));
//     if(myVerbosity_ >= TVar::DEBUG) MELAout << "Setting CKM matrix element " << i << "," << j << " to " << matEl.first.real << " + " << matEl.first.imag << "i" << endl;

//     *matElAndConjugate.first = matEl;
//     *matElAndConjugate.second = conj(matEl);
// }

// const pair<const complex<double>&, const complex<double>&> madMela::getCKMMatEl(int i, int j) const{
//     pair<complex<double>*, complex<double>*> matElAndConjugate= madMela::CKM_map(<int, int>(i, j));

//     const complex<double>& matEl = as_const(*matElAndConjugate.first);
//     const complex<double>& Conjugate = as_const(*matElAndConjugate.second);
//     return const pair<const complex<double>&, const complex<double>&>(matEl, Conjugate);
// }
