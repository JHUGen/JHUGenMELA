#ifndef madMELA
#define madMELA

#include <complex>
#include <map>
#include <vector>
#include "TVar.hh"
using namespace std;

namespace madMela{

    typedef pair<vector<int>, vector<vector<double>>> madMelaInput;

    extern "C"{
        extern struct{
            double mdl_mz,mdl_mw,mdl_mt,mdl_mb,mdl_mh,mdl_mw1,mdl_mz1,mdl_me,mdl_mmu,mdl_mta,mdl_mh1,mdl_ms,mdl_md,mdl_mu,mdl_mc,mdl_mt1;
        }mad_masses_;

        extern struct{
            double mdl_wz1,mdl_wz,mdl_wh1,mdl_wt,mdl_wh,mdl_wt1,mdl_ww,mdl_ww1;
        }widths_;

        extern struct{
            double mdl_sqrt__as,mdl_ghgg2,mdl_ghgg4,mdl_ghgg5,mdl_g__exp__2,mdl_ghgg1,mdl_ghgg3,mdl_dwh,mdl_nb__2__exp__0_25,mdl_mh__exp__2,mdl_sqrt__2,mdl_sqrt__aew,mdl_ee,mdl_mz__exp__2,mdl_ckmlambda__exp__2,mdl_ckmlambda__exp__3,
            mdl_nb__10__exp___m_40,mdl_propcorr,mdl_lambdasmeft__exp__2,mdl_ee__exp__2,mdl_mt__exp__2,mdl_mh__exp__6,mdl_mh__exp__4,mdl_mz__exp__4,mdl_mz__exp__6,mdl_mb__exp__2,mdl_mz__exp__3,mdl_propcorr__exp__2,mdl_propcorr__exp__3,
            mdl_propcorr__exp__4,mdl_ee__exp__3,mdl_sqrt__gf,mdl_vevhat,mdl_lam,mdl_sth2,mdl_cth,mdl_sqrt__sth2,mdl_sth,mdl_yb,mdl_yc,mdl_ydo,mdl_ye,mdl_ym,mdl_ys,mdl_yt,mdl_ytau,mdl_yup,mdl_vevhat__exp__2,mdl_dgf,mdl_dkh,mdl_mwsm,
            mdl_vevt,mdl_g1,mdl_gw,mdl_dmz2,mdl_dmh2,mdl_barlam,mdl_vev,mdl_mwsm__exp__6,mdl_mwsm__exp__4,mdl_mwsm__exp__2,mdl_ghaa,mdl_cth__exp__2,mdl_sth__exp__2,mdl_ghza,mdl_cth__exp__3,mdl_dg1,mdl_sth__exp__3,mdl_dgw,
            mdl_sth__exp__4,mdl_sth__exp__6,mdl_sth__exp__5,mdl_dwz,mdl_g1sh,mdl_gwsh,mdl_dmw,mdl_dwt,mdl_dww,mdl_dwhc,mdl_dwhb,mdl_dwhta,mdl_vevhat__exp__3,mdl_cth__exp__4,mdl_cth__exp__5,mdl_ckmlambda,mdl_ckma,mdl_ckmrho,mdl_ckmeta,
            mdl_cg,mdl_cw,mdl_ch,mdl_chbox,mdl_chdd,mdl_chg,mdl_chw,mdl_chb,mdl_chwb,mdl_cehre,mdl_cuhre,mdl_cdhre,mdl_cewre,mdl_cebre,mdl_cugre,mdl_cuwre,mdl_cubre,mdl_cdgre,mdl_cdwre,mdl_cdbre,mdl_chl1,mdl_chl3,mdl_che,mdl_chq1,
            mdl_chq3,mdl_chu,mdl_chd,mdl_chudre,mdl_cll,mdl_cll1,mdl_cqq1,mdl_cqq11,mdl_cqq3,mdl_cqq31,mdl_clq1,mdl_clq3,mdl_cee,mdl_cuu,mdl_cuu1,mdl_cdd,mdl_cdd1,mdl_ceu,mdl_ced,mdl_cud1,mdl_cud8,mdl_cle,mdl_clu,mdl_cld,mdl_cqe,
            mdl_cqu1,mdl_cqu8,mdl_cqd1,mdl_cqd8,mdl_cledqre,mdl_cquqd1re,mdl_cquqd11re,mdl_cquqd8re,mdl_cquqd81re,mdl_clequ1re,mdl_clequ3re,mdl_cgtil,mdl_cwtil,mdl_chgtil,mdl_chwtil,mdl_chbtil,mdl_chwbtil,mdl_cewim,mdl_cebim,
            mdl_cugim,mdl_cuwim,mdl_cubim,mdl_cdgim,mdl_cdwim,mdl_cdbim,mdl_chudim,mdl_cehim,mdl_cuhim,mdl_cdhim,mdl_cledqim,mdl_cquqd1im,mdl_cquqd8im,mdl_cquqd11im,mdl_cquqd81im,mdl_clequ1im,mdl_clequ3im,mdl_lambdasmeft,mdl_aew,
            mdl_gf,as,mdl_linearpropcorrections,mdl_ymdo,mdl_ymup,mdl_yms,mdl_ymc,mdl_ymb,mdl_ymt,mdl_yme,mdl_ymm,mdl_ymtau;
        }params_r_;

        extern struct{
            complex<double> mdl_g__exp__3,mdl_complexi,mdl_ceh,mdl_cuh,mdl_cdh,mdl_cew,mdl_ceb,mdl_cug,mdl_cuw,mdl_cub,mdl_cdg,mdl_cdw,mdl_cdb,mdl_chud,mdl_cledq,mdl_cquqd1,mdl_cquqd11,mdl_cquqd8,mdl_cquqd81,mdl_clequ1,mdl_clequ3,
            mdl_ckm1x1,mdl_ckm1x2,mdl_ckm1x3,mdl_ckm2x1,mdl_ckm2x2,mdl_ckm2x3,mdl_ckm3x1,mdl_ckm3x2,mdl_ckm3x3,mdl_conjg__ckm1x1,mdl_conjg__ckm1x2,mdl_conjg__ckm1x3,mdl_conjg__ckm2x1,mdl_conjg__ckm2x2,mdl_conjg__ckm2x3,
            mdl_conjg__ckm3x1,mdl_conjg__ckm3x2,mdl_conjg__ckm3x3;
        }params_c_;

        extern struct{
            complex<double> gc_225, gc_229, gc_234, gc_238, gc_4, gc_136,
            gc_197, gc_269, gc_270, gc_271, gc_272, gc_285, gc_319, gc_320,
            gc_321, gc_322, gc_323, gc_324, gc_325, gc_326, gc_327,
            gc_328, gc_329, gc_330, gc_352, gc_353, gc_364, gc_384, gc_387,
            gc_392, gc_395, gc_398, gc_435, gc_436, gc_437, gc_438,
            gc_464, gc_465, gc_466, gc_467, gc_708, gc_709, gc_718, gc_719,
            gc_803, gc_804, gc_813, gc_814;
        }couplings_; //you shouldn't have to edit these, but they are here for safekeeping!

        void smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void update_all_coup_();
    }
    void setDefaultMadgraphValues();
    void initialize_madMELA();

    extern TVar::VerbosityLevel myVerbosity_;
    const pair<const complex<double>&, const complex<double>&> getCKMMatEl(int i, int j);

    extern madMelaInput* madMelaCandidate;
    void setInputEvent(SimpleParticleCollection_t* pDaughters, SimpleParticleCollection_t* pAssociated, SimpleParticleCollection_t* pMothers, bool isGen);

    void computeP(double& prob, int nhel=-1);

    void set_mHiggs(double myHiggsMass, int index);
    void set_wHiggs(double myHiggsWidth, int index);
}

#endif
