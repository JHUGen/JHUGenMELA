#ifndef madMELA
#define madMELA

#include <array>
#include <map>
#include <functional>
#include "TVar.hh"
static_assert(sizeof(std::array<double,2>) == sizeof(double)*2);

namespace madMela{
    extern "C"{
        extern struct{
            double mdl_mz,mdl_mw,mdl_mt,mdl_mb,mdl_mh,mdl_mw1,
            mdl_mz1,mdl_me,mdl_mmu,mdl_mta,mdl_mh1,mdl_ms,mdl_md,mdl_mu,
            mdl_mc,mdl_mt1;
        }mad_masses_;

        extern struct{
            double mdl_wz1,mdl_wz,mdl_wh1,mdl_wt,mdl_wh,mdl_wt1,
            mdl_ww,mdl_ww1;
        }widths_;

        extern struct{
            double as,mdl_aew,mdl_barlam,mdl_cdbim,
            mdl_cdbre,mdl_cdd,mdl_cdd1,mdl_cdgim,
            mdl_cdgre,mdl_cdhim,mdl_cdhre,mdl_cdwim,
            mdl_cdwre,mdl_cebim,mdl_cebre,mdl_ced,
            mdl_cee,mdl_cehim,mdl_cehre,mdl_ceu,
            mdl_cewim,mdl_cewre,mdl_cg,mdl_cgtil,
            mdl_ch,mdl_chb,mdl_chbox,mdl_chbtil,
            mdl_chd,mdl_chdd,mdl_che,mdl_chg,
            mdl_chgtil,mdl_chl1,mdl_chl3,mdl_chq1,
            mdl_chq3,mdl_chu,mdl_chudim,mdl_chudre,
            mdl_chw,mdl_chwb,mdl_chwbtil,mdl_chwtil,
            mdl_ckma,mdl_ckmeta,mdl_ckmlambda,mdl_ckmlambda__exp__2,
            mdl_ckmlambda__exp__3,mdl_ckmrho,mdl_cld,mdl_cle,
            mdl_cledqim,mdl_cledqre,mdl_clequ1im,mdl_clequ1re,
            mdl_clequ3im,mdl_clequ3re,mdl_cll,mdl_cll1,
            mdl_clq1,mdl_clq3,mdl_clu,mdl_cqd1,
            mdl_cqd8,mdl_cqe,mdl_cqq1,mdl_cqq11,
            mdl_cqq3,mdl_cqq31,mdl_cqu1,mdl_cqu8,
            mdl_cquqd11im,mdl_cquqd11re,mdl_cquqd1im,mdl_cquqd1re,
            mdl_cquqd81im,mdl_cquqd81re,mdl_cquqd8im,mdl_cquqd8re,
            mdl_cth,mdl_cth__exp__2,mdl_cth__exp__3,mdl_cth__exp__4,
            mdl_cth__exp__5,mdl_cubim,mdl_cubre,mdl_cud1,
            mdl_cud8,mdl_cugim,mdl_cugre,mdl_cuhim,
            mdl_cuhre,mdl_cuu,mdl_cuu1,mdl_cuwim,
            mdl_cuwre,mdl_cw,mdl_cwtil,mdl_dg1,
            mdl_dgf,mdl_dgw,mdl_dkh,mdl_dmh2,
            mdl_dmw,mdl_dmz2,mdl_dwh,mdl_dwhb,
            mdl_dwhc,mdl_dwhta,mdl_dwt,mdl_dww,
            mdl_dwz,mdl_ee,mdl_ee__exp__2,mdl_ee__exp__3,
            mdl_g1,mdl_g1sh,mdl_gf,mdl_ghaa,
            mdl_ghgg1,mdl_ghgg2,mdl_ghgg3,mdl_ghgg4,
            mdl_ghgg5,mdl_ghza,mdl_gw,mdl_gwsh,
            mdl_g__exp__2,mdl_lam,mdl_lambdasmeft,mdl_lambdasmeft__exp__2,
            mdl_linearpropcorrections,mdl_mb__exp__2,mdl_mh__exp__2,mdl_mh__exp__4,
            mdl_mh__exp__6,mdl_mt__exp__2,mdl_mwsm,mdl_mwsm__exp__2,
            mdl_mwsm__exp__4,mdl_mwsm__exp__6,mdl_mz__exp__2,mdl_mz__exp__3,
            mdl_mz__exp__4,mdl_mz__exp__6,mdl_nb__10__exp___m_40,mdl_nb__2__exp__0_25,
            mdl_propcorr,mdl_propcorr__exp__2,mdl_propcorr__exp__3,mdl_propcorr__exp__4,
            mdl_sqrt__2,mdl_sqrt__aew,mdl_sqrt__as,mdl_sqrt__gf,
            mdl_sqrt__sth2,mdl_sth,mdl_sth2,mdl_sth__exp__2,
            mdl_sth__exp__3,mdl_sth__exp__4,mdl_sth__exp__5,mdl_sth__exp__6,
            mdl_vev,mdl_vevhat,mdl_vevhat__exp__2,mdl_vevhat__exp__3,
            mdl_vevt,mdl_yb,mdl_yc,mdl_ydo,
            mdl_ye,mdl_ym,mdl_ymb,mdl_ymc,
            mdl_ymdo,mdl_yme,mdl_ymm,mdl_yms,
            mdl_ymt,mdl_ymtau,mdl_ymup,mdl_ys,
            mdl_yt,mdl_ytau,mdl_yup;
        }params_r_;

        extern struct{
            std::array<double,2> mdl_cdb,mdl_cdg,mdl_cdh,mdl_cdw,
            mdl_ceb,mdl_ceh,mdl_cew,mdl_chud,
            mdl_ckm1x1,mdl_ckm1x2,mdl_ckm1x3,mdl_ckm2x1,
            mdl_ckm2x2,mdl_ckm2x3,mdl_ckm3x1,mdl_ckm3x2,
            mdl_ckm3x3,mdl_cledq,mdl_clequ1,mdl_clequ3,
            mdl_complexi,mdl_conjg__ckm1x1,mdl_conjg__ckm1x2,mdl_conjg__ckm1x3,
            mdl_conjg__ckm2x1,mdl_conjg__ckm2x2,mdl_conjg__ckm2x3,mdl_conjg__ckm3x1,
            mdl_conjg__ckm3x2,mdl_conjg__ckm3x3,mdl_cquqd1,mdl_cquqd11,
            mdl_cquqd8,mdl_cquqd81,mdl_cub,mdl_cug,
            mdl_cuh,mdl_cuw,mdl_g__exp__3;
        }params_c_;

        extern struct{
            std::array<double,2> gc_1,gc_2,gc_4,gc_12,gc_13,gc_14,gc_24,
            gc_25,gc_26,gc_27,gc_28,gc_29,gc_30,gc_40,gc_43,
            gc_135,gc_136,gc_195,gc_196,gc_197,gc_225,gc_229,gc_234,
            gc_238,gc_261,gc_263,gc_264,gc_269,gc_270,gc_271,gc_272,
            gc_285,gc_305,gc_306,gc_307,gc_309,gc_310,gc_311,gc_312,
            gc_313,gc_319,gc_320,gc_321,gc_322,gc_323,gc_324,gc_325,
            gc_326,gc_327,gc_328,gc_329,gc_330,gc_350,gc_352,gc_353,
            gc_354,gc_355,gc_356,gc_357,gc_361,gc_363,gc_364,gc_383,
            gc_384,gc_385,gc_386,gc_387,gc_388,gc_390,gc_391,gc_392,
            gc_393,gc_394,gc_395,gc_396,gc_397,gc_398,gc_429,gc_430,
            gc_431,gc_432,gc_433,gc_435,gc_436,gc_437,gc_438,gc_464,
            gc_465,gc_466,gc_467,gc_526,gc_527,gc_528,gc_535,gc_537,
            gc_549,gc_550,gc_551,gc_552,gc_553,gc_554,gc_557,gc_558,
            gc_561,gc_563,gc_574,gc_575,gc_576,gc_577,gc_578,gc_579,
            gc_580,gc_581,gc_582,gc_583,gc_602,gc_610,gc_611,gc_612,
            gc_613,gc_625,gc_626,gc_627,gc_628,gc_631,gc_632,gc_636,
            gc_637,gc_638,gc_639,gc_650,gc_651,gc_652,gc_653,gc_654,
            gc_655,gc_656,gc_657,gc_658,gc_659,gc_678,gc_684,gc_685,
            gc_686,gc_687,gc_697,gc_698,gc_699,gc_700,gc_701,gc_702,
            gc_706,gc_707,gc_708,gc_709,gc_718,gc_719,gc_720,gc_721,
            gc_722,gc_723,gc_724,gc_725,gc_726,gc_727,gc_737,gc_738,
            gc_742,gc_746,gc_747,gc_748,gc_749,gc_756,gc_757,gc_764,
            gc_765,gc_769,gc_773,gc_779,gc_780,gc_781,gc_782,gc_792,
            gc_793,gc_794,gc_795,gc_796,gc_797,gc_801,gc_802,gc_803,
            gc_804,gc_813,gc_814,gc_815,gc_816,gc_817,gc_818,gc_819,
            gc_820,gc_821,gc_822,gc_832,gc_833,gc_837,gc_841,gc_842,
            gc_843,gc_844,gc_851,gc_852,gc_859,gc_860,gc_864,gc_868,
            gc_876,gc_877,gc_878,gc_879,gc_891,gc_892,gc_893,gc_894,
            gc_897,gc_898,gc_902,gc_903,gc_904,gc_905,gc_916,gc_917,
            gc_918,gc_919,gc_920,gc_921,gc_922,gc_923,gc_924,gc_925,
            gc_944,gc_945,gc_949,gc_953,gc_954,gc_958,gc_1261,gc_1262,
            gc_1263,gc_1270,gc_1272,gc_1284,gc_1285,gc_1286,gc_1287,gc_1288,
            gc_1289,gc_1292,gc_1293,gc_1296,gc_1298,gc_1309,gc_1310,gc_1311,
            gc_1312,gc_1313,gc_1314,gc_1315,gc_1316,gc_1317,gc_1318,gc_1355,
            gc_1356,gc_1360,gc_1364,gc_1365,gc_1366,gc_1367,gc_1374,gc_1375,
            gc_1382,gc_1383,gc_1387,gc_1391,gc_1392,gc_1393,gc_1394,gc_1401,
            gc_1402;
        }couplings_; //you shouldn't have to edit these, but they are here for safekeeping!

        void ggFSIG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void ggFSIG_update_all_coup_();
        void ggFSIG_get_pdg_order_(int pdg[6][3], int allproc[3]);

        void qq4lSIG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void qq4lSIG_update_all_coup_();
        void qq4lSIG_get_pdg_order_(int pdg[6][12], int allproc[12]);

        void qq4lBKG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void qq4lBKG_update_all_coup_();
        void qq4lBKG_get_pdg_order_(int pdg[6][12], int allproc[12]);

        void qq4lBSI_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void qq4lBSI_update_all_coup_();
        void qq4lBSI_get_pdg_order_(int pdg[6][12], int allproc[12]);

    }

    struct MG_process_double {
        std::function<void(int*, int&, const int&, double*, double&, double&, int&, double&)> computeFunc;
        std::function<void()> updateFunc;
    };

    void setDefaultMadgraphValues();
    void initialize_madMELA();
    void update_all_coup(const TVar::Process& process, const TVar::Production& production);
    void smatrixhel(
        const TVar::Process& process, const TVar::Production& production,
        int pdgs[], int& procid, const int& npdg, double* p, double& alphas, 
        double& scale2, int& nhel, double& ans
        );
    extern std::map<std::pair<TVar::Process, TVar::Production>, MG_process_double>* updateMap;
}

#endif
