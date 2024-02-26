#include "Mela.h"
#include "TVar.hh"
#include "TCouplingsBase.hh"
#include "TMCFM.hh"
#include "TLorentzVector.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
namespace py = pybind11;
using namespace std;

namespace pymela{
    enum HQQ_indices{
        gHIGGS_KAPPA,
        gHIGGS_KAPPA_TILDE,

        SIZE_HQQ
    };

    enum HVV_indices{
        gHIGGS_VV_1,
        gHIGGS_VV_2,
        gHIGGS_VV_3,
        gHIGGS_VV_4,

        gHIGGS_ZA_2,
        gHIGGS_ZA_3,
        gHIGGS_ZA_4,

        gHIGGS_AA_2,
        gHIGGS_AA_3,
        gHIGGS_AA_4,

        gHIGGS_VV_1_PRIME,
        gHIGGS_VV_1_PRIME2,
        gHIGGS_VV_1_PRIME3,
        gHIGGS_VV_1_PRIME4,
        gHIGGS_VV_1_PRIME5,

        gHIGGS_VV_2_PRIME,
        gHIGGS_VV_2_PRIME2,
        gHIGGS_VV_2_PRIME3,
        gHIGGS_VV_2_PRIME4,
        gHIGGS_VV_2_PRIME5,

        gHIGGS_VV_3_PRIME,
        gHIGGS_VV_3_PRIME2,
        gHIGGS_VV_3_PRIME3,
        gHIGGS_VV_3_PRIME4,
        gHIGGS_VV_3_PRIME5,

        gHIGGS_VV_4_PRIME,
        gHIGGS_VV_4_PRIME2,
        gHIGGS_VV_4_PRIME3,
        gHIGGS_VV_4_PRIME4,
        gHIGGS_VV_4_PRIME5,

        gHIGGS_ZA_1_PRIME2,

        gHIGGS_VV_1_PRIME6,
        gHIGGS_VV_1_PRIME7,
        gHIGGS_VV_2_PRIME6,
        gHIGGS_VV_2_PRIME7,
        gHIGGS_VV_3_PRIME6,
        gHIGGS_VV_3_PRIME7,
        gHIGGS_VV_4_PRIME6,
        gHIGGS_VV_4_PRIME7,

        SIZE_HVV
    };

    enum HVV_LAMBDAQSQ_indices{
        LambdaHIGGS_QSQ_VV_1,
        LambdaHIGGS_QSQ_VV_2,
        LambdaHIGGS_QSQ_VV_3,
        LambdaHIGGS_QSQ_VV_4,

        SIZE_HVV_LAMBDAQSQ
    };

    enum HVV_CQSQ_indices{
        cLambdaHIGGS_VV_QSQ1,
        cLambdaHIGGS_VV_QSQ2,
        cLambdaHIGGS_VV_QSQ12,

        SIZE_HVV_CQSQ
    };

    enum VPff_indices{
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

    enum ZQQ_indices{
        gZPRIME_QQ_LEFT,
        gZPRIME_QQ_RIGHT,

        SIZE_ZQQ
    };

    enum ZVV_indices{
        gZPRIME_VV_1,
        gZPRIME_VV_2,

        SIZE_ZVV
    };

    enum GQQ_indices{
        gGRAVITON_QQ_LEFT,
        gGRAVITON_QQ_RIGHT,

        SIZE_GQQ
    };

    enum GGG_indices{
        gGRAVITON_GG_1,
        gGRAVITON_GG_2,
        gGRAVITON_GG_3,
        gGRAVITON_GG_4,
        gGRAVITON_GG_5,

        SIZE_GGG
    };

    enum GVV_indices{
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

    enum ATQGC_indices{
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

    enum AZff_indices{
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
}
/**
 These are functions that are pass by reference!
 They are turned into returnable functions
*/
float computeP(Mela& mela, bool useConstant=true){
    float result;
    mela.computeP(result, useConstant);
    return result;
}

float computeProdP(Mela& mela, bool useConstant=true) {
    float result;
    mela.computeProdP(result, useConstant);
    return result;
}

float computeProdDecP(Mela& mela, bool useConstant=true) {
    float result;
    mela.computeProdDecP(result, useConstant);
    return result;
}

float computeD_CP(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
    float result;
    mela.computeD_CP(myME, myType, result);
    return result;
}

float computeProdP_VH(Mela& mela, bool includeHiggsDecay, bool useConstant) {
    float result;
    mela.computeProdP_VH(result, includeHiggsDecay, useConstant);
    return result;
}

float computeProdP_ttH(Mela& mela, int topProcess, int topDecay, bool useConstant) {
    float result;
    mela.computeProdP_ttH(result, topProcess, topDecay, useConstant);
    return result;
}
float compute4FermionWeight(Mela& mela) {
    float result;
    mela.compute4FermionWeight(result);
    return result;
}
float getXPropagator(Mela& mela, TVar::ResonancePropagatorScheme scheme) {
    float result;
    mela.getXPropagator(scheme, result);
    return result;
}
float computePM4l(Mela& mela, TVar::SuperMelaSyst syst) {
    float result;
    mela.computePM4l(syst, result);
    return result;
}
float computeD_gg(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
    float result;
    mela.computeD_gg(myME, myType, result);
    return result;
}
float getConstant(Mela& mela) {
    float result;
    mela.getConstant(result);
    return result;
}
float computeDijetConvBW(Mela& mela, bool useTrueBW) {
    float result;
    mela.computeDijetConvBW(result, useTrueBW);
    return result;
}

py::array getWeightedMEArray(Mela& mela) {
    double resultarray[nmsq][nmsq] = {{0}};
    MelaIO* melaio = mela.getIORecord();
    melaio->getWeightedMEArray(resultarray);
    vector<vector<double>> result;
    for (const auto& lst : resultarray) {
        vector<double> vctr;
        for (const auto& number : lst) {
            vctr.push_back(number);
        }
        result.push_back(vctr);
    }
    py::array py_result = py::cast(result);
    return py_result;
}
py::array getUnweightedMEArray(Mela& mela) {
    double resultarray[nmsq][nmsq];
    MelaIO* melaio = mela.getIORecord();
    melaio->getUnweightedMEArray(resultarray);
    vector<vector<double>> result;
    for (const auto& lst : resultarray) {
        vector<double> vctr;
        for (const auto& number : lst) {
            vctr.push_back(number);
        }
        result.push_back(vctr);
    }
    py::array py_result = py::cast(result);
    return py_result;
}
py::tuple getPartonWeights(Mela& mela) {
    std::pair<vector<double>, vector<double>> result;
    result.first.resize(nmsq);
    result.second.resize(nmsq);
    MelaIO* melaio = mela.getIORecord();
    melaio->getPartonWeights(result.first.data(), result.second.data());
    py::tuple py_result = py::make_tuple(py::cast(result.first), py::cast(result.second));
    return py_result;
}

float getPAux(Mela& mela) {
    float result;
    mela.getPAux(result);
    return result;
}
array<float, 8> computeDecayAngles(Mela& mela) {
    array<float, 8> result;
    mela.computeDecayAngles(
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
        result[6],
        result[7]
    );
    return result;
}
array<float, 7> computeVBFAngles(Mela& mela) {
    array<float, 7> result;
    mela.computeVBFAngles(
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
        result[6]
    );
    return result;
}
array<float, 9> computeVBFAngles_ComplexBoost(Mela& mela) {
    array<float, 9> result;
    mela.computeVBFAngles_ComplexBoost(
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
        result[6],
        result[7],
        result[8]
    );
    return result;
}
array<float, 7> computeVHAngles(Mela& mela, TVar::Production prod) {
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, prod);
    array<float, 7> result;
    mela.computeVHAngles(
        result[0],
        result[1],
        result[2],
        result[3],
        result[4],
        result[5],
        result[6]
    );
    return result;
}


SimpleParticle_t particle_initializer(int id, float x, float y, float z, float e, bool ptEtaPhi=false){
    TLorentzVector vec = TLorentzVector();
    if(ptEtaPhi){
        vec.SetPtEtaPhiM(x, y, z, e);
    } else{
        vec.SetPxPyPzE(x, y, z, e);
    }
    return SimpleParticle_t(id, vec);
}

SimpleParticleCollection_t collection_initializer_from_column(std::vector<int> ids, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> e, bool ptEtaPhi=false){
    SimpleParticleCollection_t collection = SimpleParticleCollection_t();

    unsigned int i;
    for(i = 0; i < ids.size(); i++){
        collection.push_back(
            particle_initializer(ids[i], x[i], y[i], z[i], e[i], ptEtaPhi)
        );
    }
    
    return collection;
}

SimpleParticleCollection_t collection_initializer(py::list listOfParticles){

    SimpleParticleCollection_t collection = SimpleParticleCollection_t();
    for(py::handle P : listOfParticles){
        SimpleParticle_t particle = P.cast<SimpleParticle_t>();
        collection.push_back(particle);
    }
    return collection;
}

PYBIND11_MAKE_OPAQUE(SimpleParticle_t)
PYBIND11_MAKE_OPAQUE(SimpleParticleCollection_t)
PYBIND11_MODULE(Mela, m) {
    py::class_<Mela>(m, "Mela")
        .def(py::init<double, double, TVar::VerbosityLevel>())
        .def("setProcess", &Mela::setProcess)
        .def("setVerbosity", &Mela::setVerbosity)
        .def("setInputEvent", [](Mela &D, SimpleParticleCollection_t* pDaughters, SimpleParticleCollection_t* pAssociated, SimpleParticleCollection_t* pMothers, bool isGen){
            D.setInputEvent(pDaughters, pAssociated, pMothers, isGen);
        })
        // .def("setInputEvent_noMothers", [](Mela &D, SimpleParticleCollection_t* pDaughters, SimpleParticleCollection_t* pAssociated, bool isGen){
        //     D.setInputEvent(pDaughters, pAssociated, NULL, isGen);
        // })
        // .def("setInputEvent_noAssociated", [](Mela &D, SimpleParticleCollection_t* pDaughters, SimpleParticleCollection_t* pMothers, bool isGen){
        //     D.setInputEvent(pDaughters, NULL, pMothers, isGen);
        // })
        // .def("setInputEvent_onlyDaughters", [](py::object &obj, SimpleParticleCollection_t* pDaughters, bool isGen=false){
        //     Mela &D = obj.cast<Mela&>();
        //     D.setInputEvent(pDaughters, NULL, NULL, isGen);
        // })
        
        .def("resetInputEvent", &Mela::resetInputEvent)
        .def("resetMass", &Mela::resetMass)
        .def("resetWidth", &Mela::resetWidth)
        .def("resetQuarkMasses", &Mela::resetQuarkMasses)
        .def("resetMCFM_EWKParameters", &Mela::resetMCFM_EWKParameters)


        .def("getPrimaryMass", &Mela::getPrimaryMass)
        .def("getPrimaryWidth", &Mela::getPrimaryWidth)
        .def("getHiggsWidthAtPoleMass", &Mela::getHiggsWidthAtPoleMass)

        .def("getIORecord", &Mela::getIORecord)
        .def("getWeightedMEArray", &getWeightedMEArray)
        .def("getUnweightedMEArray", &getUnweightedMEArray)
        .def("getPartonWeights", &getPartonWeights)
        .def("getPAux", &getPAux)
        
        .def("computeP", &computeP)
        .def("computeProdP", &computeProdP)
        .def("computeProdDecP", &computeProdDecP)
        .def("compute4FermionWeight", &compute4FermionWeight)
        .def("getXPropagator", &getXPropagator)
        .def("computePM4l", &computePM4l)
        .def("computeD_gg", &computeD_gg)
        .def("computeProdP_VH", &computeProdP_VH)
        .def("computeProdP_ttH", &computeProdP_ttH)
        .def("computeDijetConvBW", &computeDijetConvBW)
        .def("computeD_CP", &computeD_CP)
    
        .def("getConstant", &getConstant)

        .def("computeDecayAngles", &computeDecayAngles)
        .def("computeVBFAngles", &computeVBFAngles)
        .def("computeVBFAngles_ComplexBoost", &computeVBFAngles_ComplexBoost)
        .def("computeVHAngles", &computeVHAngles)
       
        .def_readwrite("differentiate_HWW_HZZ", &Mela::differentiate_HWW_HZZ)

        .def("selfDHggcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HGG, 2}, (const double*) &D.selfDHggcoupl, obj);
        })


        .def("selfDHg4g4coupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HGG, 2}, (const double*) &D.selfDHg4g4coupl, obj);
        })


        .def("selfDHqqcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HQQ, 2}, (const double*) &D.selfDHqqcoupl, obj);
        })


        .def("selfDHbbcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HQQ, 2}, (const double*) &D.selfDHbbcoupl, obj);
        })


        .def("selfDHttcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HQQ, 2}, (const double*) &D.selfDHttcoupl, obj);
        })


        .def("selfDHb4b4coupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HQQ, 2}, (const double*) &D.selfDHb4b4coupl, obj);
        })


        .def("selfDHt4t4coupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HQQ, 2}, (const double*) &D.selfDHt4t4coupl, obj);
        })


        .def("selfDHzzcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV, 2}, (const double*) &D.selfDHzzcoupl, obj);
        })


        .def("selfDHwwcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV, 2}, (const double*) &D.selfDHwwcoupl, obj);
        })


        .def("selfDHzzLambda_qsq", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV_LAMBDAQSQ, SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
        })


        .def("selfDHwwLambda_qsq", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV_LAMBDAQSQ, SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
        })


        .def("selfDHzzCLambda_qsq", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzCLambda_qsq, obj);
        })


        .def("selfDHwwCLambda_qsq", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{nSupportedHiggses, SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
        })


        .def("selfDHzzpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_HVV, 2}, (const double*) &D.selfDHzzpcoupl, obj);
        })


        .def("selfDHzpzpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_HVV, 2}, (const double*) &D.selfDHzpzpcoupl, obj);
        })


        .def("selfDZpffcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_Vpff, 2}, (const double*) &D.selfDZpffcoupl, obj);
        })


        .def("selfDHwwpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_HVV, 2}, (const double*) &D.selfDHwwpcoupl, obj);
        })


        .def("selfDHwpwpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_HVV, 2}, (const double*) &D.selfDHwpwpcoupl, obj);
        })


        .def("selfDWpffcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_Vpff, 2}, (const double*) &D.selfDWpffcoupl, obj);
        })


        .def("selfDZqqcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_ZQQ, 2}, (const double*) &D.selfDZqqcoupl, obj);
        })


        .def("selfDZvvcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_ZVV, 2}, (const double*) &D.selfDZvvcoupl, obj);
        })


        .def("selfDGqqcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_GQQ, 2}, (const double*) &D.selfDGqqcoupl, obj);
        })


        .def("selfDGggcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_GGG, 2}, (const double*) &D.selfDGggcoupl, obj);
        })


        .def("selfDGvvcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_GVV, 2}, (const double*) &D.selfDGvvcoupl, obj);
        })


        .def("selfDGvvpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_GVV, 2}, (const double*) &D.selfDGvvpcoupl, obj);
        })


        .def("selfDGvpvpcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_GVV, 2}, (const double*) &D.selfDGvpvpcoupl, obj);
        })


        .def("selfDaTQGCcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_ATQGC, 2}, (const double*) &D.selfDaTQGCcoupl, obj);
        })


        .def("selfDAZffcoupl", [](py::object &obj){
            Mela &D = obj.cast<Mela&>();
            return py::array_t<double>(std::vector<int>{SIZE_AZff, 2}, (const double*) &D.selfDAZffcoupl, obj);
        })


        .def_readwrite("selfDM_Zprime", &Mela::selfDM_Zprime)
        .def_readwrite("selfDGa_Zprime", &Mela::selfDGa_Zprime)
        .def_readwrite("selfDM_Wprime", &Mela::selfDM_Wprime)
        .def_readwrite("selfDGa_Wprime", &Mela::selfDGa_Wprime);

    py::class_<SimpleParticle_t>(m, "SimpleParticle_t")
        .def(py::init(&particle_initializer), py::arg("id"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("e"), py::arg("ptEtaPhi") = false)
        .def_property_readonly("id", [](SimpleParticle_t& P){
            return P.first;
        })
        .def_property_readonly("vector", [](SimpleParticle_t& P){
            P.second.Print();
            return py::make_tuple(P.second.Px(), P.second.Py(), P.second.Pz(), P.second.E());
        })
        .def("__repr__",[](SimpleParticle_t& P){
            return "Particle with id " + std::to_string(P.first);
        });


    py::class_<SimpleParticleCollection_t>(m, "SimpleParticleCollection_t")
        .def(py::init(&collection_initializer_from_column), py::arg("ids"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("e"), py::arg("ptEtaPhi") = false)
        .def(py::init(&collection_initializer))
        .def(py::init())
        .def("add_particle", [](SimpleParticleCollection_t &C, SimpleParticle_t P){
            C.push_back(P);
        })
        .def("__iter__", [](const SimpleParticleCollection_t &C){
            return py::make_iterator(C.begin(), C.end());
        }, py::keep_alive<0, 1>())
        .def("toList", [](SimpleParticleCollection_t &C){
            py::list list_type = py::cast(C);
            return list_type;
        });

    py::enum_<TVar::VerbosityLevel>(m, "VerbosityLevel")
        .value("SILENT",TVar::SILENT)
        .value("ERROR",TVar::ERROR)
        .value("INFO",TVar::INFO)
        .value("DEBUG",TVar::DEBUG)
        .value("DEBUG_VERBOSE",TVar::DEBUG_VERBOSE)
        .value("DEBUG_MECHECK",TVar::DEBUG_MECHECK);

    py::enum_<TVar::MatrixElement>(m, "MatrixElement")
        .value("MCFM",TVar::MCFM)
        .value("JHUGen",TVar::JHUGen)
        .value("ANALYTICAL",TVar::ANALYTICAL);
    
    py::enum_<TVar::Production>(m, "Production")
        .value("ZZGG",TVar::ZZGG)
        .value("ZZQQB",TVar::ZZQQB)
        .value("ZZQQB_STU",TVar::ZZQQB_STU) // Should be the same as ZZQQB, just for crosscheck
        .value("ZZINDEPENDENT",TVar::ZZINDEPENDENT)
        .value("ttH",TVar::ttH) // ttH
        .value("bbH",TVar::bbH) // bbH
        .value("JQCD",TVar::JQCD) // ? + 1 jet
        .value("JJQCD",TVar::JJQCD) // SBF
        .value("JJVBF",TVar::JJVBF) // VBF
        .value("JJEW",TVar::JJEW) // VBF+VH (had.)
        .value("JJEWQCD",TVar::JJEWQCD) // VBF+VH+QCD, all hadronic
        .value("Had_ZH",TVar::Had_ZH) // ZH, Z->uu/dd
        .value("Had_WH",TVar::Had_WH) // W(+/-)H, W->ud
        .value("Lep_ZH",TVar::Lep_ZH) // ZH, Z->ll/nunu
        .value("Lep_WH",TVar::Lep_WH) // W(+/-)H, W->lnu
        .value("ZZQQB_S",TVar::ZZQQB_S)
        .value("JJQCD_S",TVar::JJQCD_S)
        .value("JJVBF_S",TVar::JJVBF_S)
        .value("JJEW_S",TVar::JJEW_S)
        .value("JJEWQCD_S",TVar::JJEWQCD_S)
        .value("Had_ZH_S",TVar::Had_ZH_S)
        .value("Had_WH_S",TVar::Had_WH_S)
        .value("Lep_ZH_S",TVar::Lep_ZH_S)
        .value("Lep_WH_S",TVar::Lep_WH_S)
        .value("ZZQQB_TU",TVar::ZZQQB_TU)
        .value("JJQCD_TU",TVar::JJQCD_TU)
        .value("JJVBF_TU",TVar::JJVBF_TU)
        .value("JJEW_TU",TVar::JJEW_TU)
        .value("JJEWQCD_TU",TVar::JJEWQCD_TU)
        .value("Had_ZH_TU",TVar::Had_ZH_TU)
        .value("Had_WH_TU",TVar::Had_WH_TU)
        .value("Lep_ZH_TU",TVar::Lep_ZH_TU)
        .value("Lep_WH_TU",TVar::Lep_WH_TU)
        .value("GammaH",TVar::GammaH) // gammaH, stable A (could implement S and TU in the future
        .value("nProductions",TVar::nProductions);
    
    py::enum_<pymela::HQQ_indices>(m, "HQQ_indices")
        .value("gHIGGS_KAPPA", pymela::gHIGGS_KAPPA)
        .value("gHIGGS_KAPPA_TILDE", pymela::gHIGGS_KAPPA_TILDE)
        .value("SIZE_HQQ", pymela::SIZE_HQQ);

    py::enum_<pymela::HVV_indices>(m, "HVV_indices")
        .value("gHIGGS_VV_1", pymela::gHIGGS_VV_1)
        .value("gHIGGS_VV_2", pymela::gHIGGS_VV_2)
        .value("gHIGGS_VV_3", pymela::gHIGGS_VV_3)
        .value("gHIGGS_VV_4", pymela::gHIGGS_VV_4)
        .value("gHIGGS_ZA_2", pymela::gHIGGS_ZA_2)
        .value("gHIGGS_ZA_3", pymela::gHIGGS_ZA_3)
        .value("gHIGGS_ZA_4", pymela::gHIGGS_ZA_4)
        .value("gHIGGS_AA_2", pymela::gHIGGS_AA_2)
        .value("gHIGGS_AA_3", pymela::gHIGGS_AA_3)
        .value("gHIGGS_AA_4", pymela::gHIGGS_AA_4)
        .value("gHIGGS_VV_1_PRIME", pymela::gHIGGS_VV_1_PRIME)
        .value("gHIGGS_VV_1_PRIME2", pymela::gHIGGS_VV_1_PRIME2)
        .value("gHIGGS_VV_1_PRIME3", pymela::gHIGGS_VV_1_PRIME3)
        .value("gHIGGS_VV_1_PRIME4", pymela::gHIGGS_VV_1_PRIME4)
        .value("gHIGGS_VV_1_PRIME5", pymela::gHIGGS_VV_1_PRIME5)
        .value("gHIGGS_VV_2_PRIME", pymela::gHIGGS_VV_2_PRIME)
        .value("gHIGGS_VV_2_PRIME2", pymela::gHIGGS_VV_2_PRIME2)
        .value("gHIGGS_VV_2_PRIME3", pymela::gHIGGS_VV_2_PRIME3)
        .value("gHIGGS_VV_2_PRIME4", pymela::gHIGGS_VV_2_PRIME4)
        .value("gHIGGS_VV_2_PRIME5", pymela::gHIGGS_VV_2_PRIME5)
        .value("gHIGGS_VV_3_PRIME", pymela::gHIGGS_VV_3_PRIME)
        .value("gHIGGS_VV_3_PRIME2", pymela::gHIGGS_VV_3_PRIME2)
        .value("gHIGGS_VV_3_PRIME3", pymela::gHIGGS_VV_3_PRIME3)
        .value("gHIGGS_VV_3_PRIME4", pymela::gHIGGS_VV_3_PRIME4)
        .value("gHIGGS_VV_3_PRIME5", pymela::gHIGGS_VV_3_PRIME5)
        .value("gHIGGS_VV_4_PRIME", pymela::gHIGGS_VV_4_PRIME)
        .value("gHIGGS_VV_4_PRIME2", pymela::gHIGGS_VV_4_PRIME2)
        .value("gHIGGS_VV_4_PRIME3", pymela::gHIGGS_VV_4_PRIME3)
        .value("gHIGGS_VV_4_PRIME4", pymela::gHIGGS_VV_4_PRIME4)
        .value("gHIGGS_VV_4_PRIME5", pymela::gHIGGS_VV_4_PRIME5)
        .value("gHIGGS_ZA_1_PRIME2", pymela::gHIGGS_ZA_1_PRIME2)
        .value("gHIGGS_VV_1_PRIME6", pymela::gHIGGS_VV_1_PRIME6)
        .value("gHIGGS_VV_1_PRIME7", pymela::gHIGGS_VV_1_PRIME7)
        .value("gHIGGS_VV_2_PRIME6", pymela::gHIGGS_VV_2_PRIME6)
        .value("gHIGGS_VV_2_PRIME7", pymela::gHIGGS_VV_2_PRIME7)
        .value("gHIGGS_VV_3_PRIME6", pymela::gHIGGS_VV_3_PRIME6)
        .value("gHIGGS_VV_3_PRIME7", pymela::gHIGGS_VV_3_PRIME7)
        .value("gHIGGS_VV_4_PRIME6", pymela::gHIGGS_VV_4_PRIME6)
        .value("gHIGGS_VV_4_PRIME7", pymela::gHIGGS_VV_4_PRIME7)
        .value("SIZE_HVV", pymela::SIZE_HVV);

    py::enum_<pymela::HVV_LAMBDAQSQ_indices>(m, "HVV_LAMBDAQSQ_indices")
        .value("LambdaHIGGS_QSQ_VV_1", pymela::LambdaHIGGS_QSQ_VV_1)
        .value("LambdaHIGGS_QSQ_VV_2", pymela::LambdaHIGGS_QSQ_VV_2)
        .value("LambdaHIGGS_QSQ_VV_3", pymela::LambdaHIGGS_QSQ_VV_3)
        .value("LambdaHIGGS_QSQ_VV_4", pymela::LambdaHIGGS_QSQ_VV_4)
        .value("SIZE_HVV_LAMBDAQSQ", pymela::SIZE_HVV_LAMBDAQSQ);

    py::enum_<pymela::HVV_CQSQ_indices>(m, "HVV_CQSQ_indices")
        .value("cLambdaHIGGS_VV_QSQ1", pymela::cLambdaHIGGS_VV_QSQ1)
        .value("cLambdaHIGGS_VV_QSQ2", pymela::cLambdaHIGGS_VV_QSQ2)
        .value("cLambdaHIGGS_VV_QSQ12", pymela::cLambdaHIGGS_VV_QSQ12)
        .value("SIZE_HVV_CQSQ", pymela::SIZE_HVV_CQSQ);

    py::enum_<pymela::VPff_indices>(m, "VPff_indices")
        .value("gHIGGS_Vp_El_left", pymela::gHIGGS_Vp_El_left)
        .value("gHIGGS_Vp_El_right", pymela::gHIGGS_Vp_El_right)
        .value("gHIGGS_Vp_Mu_left", pymela::gHIGGS_Vp_Mu_left)
        .value("gHIGGS_Vp_Mu_right", pymela::gHIGGS_Vp_Mu_right)
        .value("gHIGGS_Vp_Ta_left", pymela::gHIGGS_Vp_Ta_left)
        .value("gHIGGS_Vp_Ta_right", pymela::gHIGGS_Vp_Ta_right)
        .value("gHIGGS_Vp_NuE_left", pymela::gHIGGS_Vp_NuE_left)
        .value("gHIGGS_Vp_NuE_right", pymela::gHIGGS_Vp_NuE_right)
        .value("gHIGGS_Vp_Dn_left", pymela::gHIGGS_Vp_Dn_left)
        .value("gHIGGS_Vp_Dn_right", pymela::gHIGGS_Vp_Dn_right)
        .value("gHIGGS_Vp_Up_left", pymela::gHIGGS_Vp_Up_left)
        .value("gHIGGS_Vp_Up_right", pymela::gHIGGS_Vp_Up_right)
        .value("gHIGGS_Vp_Str_left", pymela::gHIGGS_Vp_Str_left)
        .value("gHIGGS_Vp_Str_right", pymela::gHIGGS_Vp_Str_right)
        .value("gHIGGS_Vp_Chm_left", pymela::gHIGGS_Vp_Chm_left)
        .value("gHIGGS_Vp_Chm_right", pymela::gHIGGS_Vp_Chm_right)
        .value("gHIGGS_Vp_Bot_left", pymela::gHIGGS_Vp_Bot_left)
        .value("gHIGGS_Vp_Bot_right", pymela::gHIGGS_Vp_Bot_right)
        .value("gHIGGS_Vp_Top_left", pymela::gHIGGS_Vp_Top_left)
        .value("gHIGGS_Vp_Top_right", pymela::gHIGGS_Vp_Top_right)
        .value("SIZE_Vpff", pymela::SIZE_Vpff);

    py::enum_<pymela::ZQQ_indices>(m, "ZQQ_indices")
        .value("gZPRIME_QQ_LEFT", pymela::gZPRIME_QQ_LEFT)
        .value("gZPRIME_QQ_RIGHT", pymela::gZPRIME_QQ_RIGHT)
        .value("SIZE_ZQQ", pymela::SIZE_ZQQ);

    py::enum_<pymela::ZVV_indices>(m, "ZVV_indices")
        .value("gZPRIME_VV_1", pymela::gZPRIME_VV_1)
        .value("gZPRIME_VV_2", pymela::gZPRIME_VV_2)
        .value("SIZE_ZVV", pymela::SIZE_ZVV);

    py::enum_<pymela::GQQ_indices>(m, "GQQ_indices")
        .value("gGRAVITON_QQ_LEFT", pymela::gGRAVITON_QQ_LEFT)
        .value("gGRAVITON_QQ_RIGHT", pymela::gGRAVITON_QQ_RIGHT)
        .value("SIZE_GQQ", pymela::SIZE_GQQ);

    py::enum_<pymela::GGG_indices>(m, "GGG_indices")
        .value("gGRAVITON_GG_1", pymela::gGRAVITON_GG_1)
        .value("gGRAVITON_GG_2", pymela::gGRAVITON_GG_2)
        .value("gGRAVITON_GG_3", pymela::gGRAVITON_GG_3)
        .value("gGRAVITON_GG_4", pymela::gGRAVITON_GG_4)
        .value("gGRAVITON_GG_5", pymela::gGRAVITON_GG_5)
        .value("SIZE_GGG", pymela::SIZE_GGG);

    py::enum_<pymela::GVV_indices>(m, "GVV_indices")
        .value("gGRAVITON_VV_1", pymela::gGRAVITON_VV_1)
        .value("gGRAVITON_VV_2", pymela::gGRAVITON_VV_2)
        .value("gGRAVITON_VV_3", pymela::gGRAVITON_VV_3)
        .value("gGRAVITON_VV_4", pymela::gGRAVITON_VV_4)
        .value("gGRAVITON_VV_5", pymela::gGRAVITON_VV_5)
        .value("gGRAVITON_VV_6", pymela::gGRAVITON_VV_6)
        .value("gGRAVITON_VV_7", pymela::gGRAVITON_VV_7)
        .value("gGRAVITON_VV_8", pymela::gGRAVITON_VV_8)
        .value("gGRAVITON_VV_9", pymela::gGRAVITON_VV_9)
        .value("gGRAVITON_VV_10", pymela::gGRAVITON_VV_10)
        .value("gGRAVITON_ZA_1", pymela::gGRAVITON_ZA_1)
        .value("gGRAVITON_ZA_2", pymela::gGRAVITON_ZA_2)
        .value("gGRAVITON_ZA_3", pymela::gGRAVITON_ZA_3)
        .value("gGRAVITON_ZA_4", pymela::gGRAVITON_ZA_4)
        .value("gGRAVITON_ZA_8", pymela::gGRAVITON_ZA_8)
        .value("gGRAVITON_AA_1", pymela::gGRAVITON_AA_1)
        .value("gGRAVITON_AA_2", pymela::gGRAVITON_AA_2)
        .value("gGRAVITON_AA_3", pymela::gGRAVITON_AA_3)
        .value("gGRAVITON_AA_4", pymela::gGRAVITON_AA_4)
        .value("gGRAVITON_AA_8", pymela::gGRAVITON_AA_8)
        .value("SIZE_GVV", pymela::SIZE_GVV);

    py::enum_<pymela::ATQGC_indices>(m, "ATQGC_indices")
        .value("gATQGC_dVA", pymela::gATQGC_dVA)
        .value("gATQGC_dPA", pymela::gATQGC_dPA)
        .value("gATQGC_dMA", pymela::gATQGC_dMA)
        .value("gATQGC_dFourA", pymela::gATQGC_dFourA)
        .value("gATQGC_dVZ", pymela::gATQGC_dVZ)
        .value("gATQGC_dPZ", pymela::gATQGC_dPZ)
        .value("gATQGC_dMZ", pymela::gATQGC_dMZ)
        .value("gATQGC_dFourZ", pymela::gATQGC_dFourZ)
        .value("gATQGC_dAAWpWm", pymela::gATQGC_dAAWpWm)
        .value("gATQGC_dZAWpWm", pymela::gATQGC_dZAWpWm)
        .value("gATQGC_dZZWpWm", pymela::gATQGC_dZZWpWm)
        .value("SIZE_ATQGC", pymela::SIZE_ATQGC);

    py::enum_<pymela::AZff_indices>(m, "AZff_indices")
        .value("gAZff_ZllRH", pymela::gAZff_ZllRH)
        .value("gAZff_ZllLH", pymela::gAZff_ZllLH)
        .value("gAZff_ZuuRH", pymela::gAZff_ZuuRH)
        .value("gAZff_ZuuLH", pymela::gAZff_ZuuLH)
        .value("gAZff_ZddRH", pymela::gAZff_ZddRH)
        .value("gAZff_ZddLH", pymela::gAZff_ZddLH)
        .value("gAZff_ZnunuRH", pymela::gAZff_ZnunuRH)
        .value("gAZff_ZnunuLH", pymela::gAZff_ZnunuLH)
        .value("gAZff_uZRH", pymela::gAZff_uZRH)
        .value("gAZff_uZLH", pymela::gAZff_uZLH)
        .value("gAZff_dZRH", pymela::gAZff_dZRH)
        .value("gAZff_dZLH", pymela::gAZff_dZLH)
        .value("SIZE_AZff", pymela::SIZE_AZff);


}