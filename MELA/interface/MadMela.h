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

    typedef std::pair<std::function<void(int*, int&, const int&, double*, double&, double&, int&, double&)>, std::function<void()>> MG_process_double;

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
