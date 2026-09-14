/**
 * @file MadMela.h
 * @brief This is the code that interfaces with the generated Madgraph
 * matrix elements in fortran for specified Madgraph processes.
 * 
 * MadMELA was made possible primarily through the efforts of
 * Jeffrey Davis and Mohit Srivastav, as well as the rest of
 * the JHUGenMELA development team. Code written
 * by Mohit Srivastav.
 * 
 * @author Mohit Srivastav
*/

#ifndef madMELA
#define madMELA

#include <array>
#include <map>
#include <vector>
#include <functional>
#include "TVar.hh"
static_assert(sizeof(std::array<double,2>) == sizeof(double)*2);

namespace madMela{
    extern "C"{
        extern struct{
            double mdl_mt1,mdl_mh,mdl_mt,mdl_mw1,mdl_mz1,mdl_mz,
            mdl_mw,mdl_mh1;
        }mad_masses_;

        extern struct{
            double mdl_wz1,mdl_wz,mdl_wt1,mdl_wh1,mdl_ww1,mdl_wh,
            mdl_ww,mdl_wt;
        }widths_;

        extern struct{
            double as,mdl_aew,mdl_barlam,mdl_ch,
            mdl_chb,mdl_chbox,mdl_chbtil,mdl_chdd,
            mdl_chg,mdl_chj3,mdl_chl3,mdl_chw,
            mdl_chwb,mdl_chwbtil,mdl_chwtil,mdl_conjg__cthh,
            mdl_cth,mdl_cth__exp__2,mdl_cth__exp__3,mdl_cw,
            mdl_cwtil,mdl_dg1,mdl_dgf,mdl_dgw,
            mdl_dkh,mdl_dmh2,mdl_dwh,mdl_dwhb,
            mdl_dwhc,mdl_dwhta,mdl_dww,mdl_ee,
            mdl_ee__exp__2,mdl_ee__exp__3,mdl_g1,mdl_g1sh,
            mdl_gf,mdl_ghaa,mdl_ghgg1,mdl_ghgg2,
            mdl_ghgg3,mdl_ghgg4,mdl_ghgg5,mdl_ghza,
            mdl_gw,mdl_gwsh,mdl_g__exp__2,mdl_lam,
            mdl_lambdasmeft,mdl_lambdasmeft__exp__2,mdl_mh__exp__2,mdl_mh__exp__4,
            mdl_mh__exp__6,mdl_mt__exp__2,mdl_mwsm,mdl_mwsm__exp__2,
            mdl_mwsm__exp__4,mdl_mwsm__exp__6,mdl_mw__exp__2,mdl_mz__exp__2,
            mdl_mz__exp__3,mdl_mz__exp__4,mdl_mz__exp__6,mdl_nb__10__exp___m_40,
            mdl_nb__2__exp__0_25,mdl_sqrt__2,mdl_sqrt__aew,mdl_sqrt__as,
            mdl_sqrt__gf,mdl_sqrt__sth2,mdl_sth,mdl_sth2,
            mdl_sth__exp__2,mdl_sth__exp__3,mdl_sth__exp__4,mdl_sth__exp__5,
            mdl_sth__exp__6,mdl_vev,mdl_vevhat,mdl_vevhat__exp__2,
            mdl_vevhat__exp__3,mdl_vevt,mdl_ymt,mdl_yt;
        }params_r_;

        extern struct{
            std::array<double,2> mdl_complexi,mdl_g__exp__3,mdl_yt0;
        }params_c_;

        // void ggFSIG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        // void ggFSIG_update_all_coup_();
        // void ggFSIG_get_pdg_order_(int pdg[6][3], int allproc[3]);

        // void qq4lSIG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        // void qq4lSIG_update_all_coup_();
        // void qq4lSIG_get_pdg_order_(int pdg[6][30], int allproc[30]);

        // void qq4lBKG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        // void qq4lBKG_update_all_coup_();
        // void qq4lBKG_get_pdg_order_(int pdg[6][30], int allproc[30]);

        // void qq4lBSI_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        // void qq4lBSI_update_all_coup_();
        // void qq4lBSI_get_pdg_order_(int pdg[6][30], int allproc[30]);

        void VBSSIG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void VBSSIG_update_all_coup_();
        void VBSSIG_get_pdg_order_(int pdg[8][261], int allproc[261]);

        void VBSBSI_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void VBSBSI_update_all_coup_();
        void VBSBSI_get_pdg_order_(int pdg[8][261], int allproc[261]);

        void VBSBKG_smatrixhel_(int pdgs[], int& procid, const int& npdg, double* p, double& alphas, double& scale2, int& nhel, double& ans);
        void VBSBKG_update_all_coup_();
        void VBSBKG_get_pdg_order_(int pdg[8][261], int allproc[261]);
    }

    typedef std::pair<std::function<void(int*, int&, const int&, double*, double&, double&, int&, double&)>, std::function<void()>> MG_process_double;
    /**
     * @brief Sets the Madgraph values to their defaults.
     * This function is called after every call to a compute function
     */
    void setDefaultMadgraphValues();

    /**
     * @brief Initializes all of the values in Madgraph for proper usage.
     * These include values like alpha, the fermi constant, 
     * quark masses, and force carrier boson masses (
     * including the Higgs)
     */
    void initialize_madMELA();

    /**
     * @brief This is the function that redirects the 
     * SMEFTSim calculations to the correct process and production mode
     * 
     * @param[in] process The process being used (which correlate to
     * signal, background, signal+background+interference)
     * @param[in] production The production mode
     */
    void update_all_coup(const TVar::Process& process, const TVar::Production& production);

    /**
     * @brief This is the function that redirects the SMEFTsim
     * matrix element calculation to the correct smatrix
     * function call
     * 
     * @param[in] process The process being used (which correlate to
     * signal, background, signal+background+interference)
     * @param[in] production The production mode
     * @param[in] pdgs The particle ID of the particles being used
     * @param[in] procid The process ID of the process (This 
     * is a useless variable, it is never utilized
     * and only exists as a relic of the autogenerated
     * Madgraph function)
     * @param[in] npdg The number of particles that are being passed.
     * In the current use-case, this is always 6 (2 mother 
     * particles and 4 final-state leptons)
     * @param[in] p The 4-vectors to be passed
     * in Fortran column-major format.
     * These vectors are passed as `<E, px, py, pz>`.
     * @param[in] alphas The value of alpha strong
     * @param[in] scale2 Another irrelevant variable.
     * @param[in] nhel This should be -1 to sum over all possible,
     * helicities. However, you can set this to be a specific number
     * that corresponds to a specific index in the subprocess
     * @param[out] ans This is where the matrix element goes
     */
    void smatrixhel(
        const TVar::Process& process, const TVar::Production& production,
        int pdgs[], int& procid, const int& npdg, double* p, double& alphas, 
        double& scale2, int& nhel, double& ans
        );
    
    std::vector<std::vector<int>> get_pdg_order(
        const TVar::Process& process, const TVar::Production& production
        );
    
    extern std::map<std::pair<TVar::Process, TVar::Production>, MG_process_double>* updateMap;
}

#endif
