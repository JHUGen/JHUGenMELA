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

    enum HGG_indices{
        gHIGGS_GG_2,
        gHIGGS_GG_3,
        gHIGGS_GG_4,
    
        SIZE_HGG
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
        .def(py::init<double, double>())
        .def(py::init<double>())
        .def(py::init<>())
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
            return py::array_t<int>(std::vector<int>{nSupportedHiggses, SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
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
        .def_readwrite("selfDGa_Wprime", &Mela::selfDGa_Wprime)


        //Here be couplings
        .def_property(
            "ghg2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHggcoupl[0][gHIGGS_GG_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHggcoupl[0][gHIGGS_GG_2][0] = coupl[0];
                    D.selfDHggcoupl[0][gHIGGS_GG_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghg3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHggcoupl[0][gHIGGS_GG_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHggcoupl[0][gHIGGS_GG_3][0] = coupl[0];
                    D.selfDHggcoupl[0][gHIGGS_GG_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghg4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHggcoupl[0][gHIGGS_GG_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHggcoupl[0][gHIGGS_GG_4][0] = coupl[0];
                    D.selfDHggcoupl[0][gHIGGS_GG_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_ZA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_ZA_2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_ZA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_ZA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_ZA_3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_ZA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_ZA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_ZA_4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_ZA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghgsgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_AA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_AA_2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_AA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghgsgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_AA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_AA_3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_AA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghgsgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_AA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_AA_4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_AA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME5][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME5][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME5][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME3][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME4][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME5][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_4_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzgs1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME6][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz1_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME7][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME6][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz2_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME7][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_2_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz3_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghz4_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzzcoupl[0][gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cz_q1sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z11", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z12", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z13", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z14", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cz_q2sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z21", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z22", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z23", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z24", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cz_q12sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<int>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z01", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z02", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z03", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_z04", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME2][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME4][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME5][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME2][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME3][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME4][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME5][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME2][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME3][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME4][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME5][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME2][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME3][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME4][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME5][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_4_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME6][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw1_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME7][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME6][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw2_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME7][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_2_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw3_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghw4_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHwwcoupl[0][gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cw_q1sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w11", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w12", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w13", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w14", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ1) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cw_q2sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w21", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w22", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w23", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w24", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ2) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cw_q12sq", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    return array_val.at(0,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwCLambda_qsq, obj);
                    array_val.mutable_at(0,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w01", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_1,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w02", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_2,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w03", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_3,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "Lambda_w04", 
            py::cpp_function(
                [](py::object &obj){
                    Mela& D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    return array_val.at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12);
                }),
            py::cpp_function(
                [](py::object &obj, double coupl){
                    Mela &D = obj.cast<Mela&>();
                    py::array_t array_val = py::array_t<double>(std::vector<int>{0,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj);
                    array_val.mutable_at(0,LambdaHIGGS_QSQ_VV_4,cLambdaHIGGS_VV_QSQ12) = coupl;
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHqqcoupl[0][gHIGGS_KAPPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHqqcoupl[0][gHIGGS_KAPPA][0] = coupl[0];
                    D.selfDHqqcoupl[0][gHIGGS_KAPPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_tilde", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHqqcoupl[0][gHIGGS_KAPPA_TILDE], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHqqcoupl[0][gHIGGS_KAPPA_TILDE][0] = coupl[0];
                    D.selfDHqqcoupl[0][gHIGGS_KAPPA_TILDE][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_top", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHttcoupl[0][gHIGGS_KAPPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHttcoupl[0][gHIGGS_KAPPA][0] = coupl[0];
                    D.selfDHttcoupl[0][gHIGGS_KAPPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_top_tilde", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHttcoupl[0][gHIGGS_KAPPA_TILDE], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHttcoupl[0][gHIGGS_KAPPA_TILDE][0] = coupl[0];
                    D.selfDHttcoupl[0][gHIGGS_KAPPA_TILDE][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_bot", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHbbcoupl[0][gHIGGS_KAPPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHbbcoupl[0][gHIGGS_KAPPA][0] = coupl[0];
                    D.selfDHbbcoupl[0][gHIGGS_KAPPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_bot_tilde", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHbbcoupl[0][gHIGGS_KAPPA_TILDE], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHbbcoupl[0][gHIGGS_KAPPA_TILDE][0] = coupl[0];
                    D.selfDHbbcoupl[0][gHIGGS_KAPPA_TILDE][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_4gen_top", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHt4t4coupl[0][gHIGGS_KAPPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHt4t4coupl[0][gHIGGS_KAPPA][0] = coupl[0];
                    D.selfDHt4t4coupl[0][gHIGGS_KAPPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_4gen_top_tilde", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHt4t4coupl[0][gHIGGS_KAPPA_TILDE], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHt4t4coupl[0][gHIGGS_KAPPA_TILDE][0] = coupl[0];
                    D.selfDHt4t4coupl[0][gHIGGS_KAPPA_TILDE][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_4gen_bot", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHb4b4coupl[0][gHIGGS_KAPPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHb4b4coupl[0][gHIGGS_KAPPA][0] = coupl[0];
                    D.selfDHb4b4coupl[0][gHIGGS_KAPPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "kappa_4gen_bot_tilde", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHb4b4coupl[0][gHIGGS_KAPPA_TILDE], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHb4b4coupl[0][gHIGGS_KAPPA_TILDE][0] = coupl[0];
                    D.selfDHb4b4coupl[0][gHIGGS_KAPPA_TILDE][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_ZA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_ZA_2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_ZA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_ZA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_ZA_3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_ZA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_ZA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_ZA_4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_ZA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME5][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME5][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME5][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME3][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME4][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME5][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_4_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpgs1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_ZA_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_ZA_1_PRIME2][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_ZA_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME6][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp1_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME7][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_1_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME6][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp2_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME7][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_2_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp3_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzzp4_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzzpcoupl[gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME2][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME3][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME4][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME5][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME2][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME3][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME4][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME5][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME2][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME3][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME4][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME5][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME2][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME3][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME4][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME5][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_4_PRIME5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME6][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp1_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME7][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_1_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME6][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp2_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME7][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_2_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp3_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghzpzp4_prime7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7][0] = coupl[0];
                    D.selfDHzpzpcoupl[gHIGGS_VV_3_PRIME7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_El_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_El_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_El_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_El_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_El_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_El_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_El_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_El_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Mu_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Mu_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Mu_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Mu_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Mu_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Mu_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Ta_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Ta_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Ta_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Ta_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Ta_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Ta_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_NuE_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_NuE_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_NuE_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_NuE_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_NuE_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_NuE_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_NuE_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Dn_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Dn_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Dn_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Dn_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Dn_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Dn_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Up_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Up_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Up_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Up_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Up_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Up_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Up_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Up_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Str_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Str_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Str_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Str_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Str_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Str_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Str_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Str_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Chm_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Chm_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Chm_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Chm_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Chm_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Chm_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Bot_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Bot_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Bot_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Bot_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Bot_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Bot_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Top_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Top_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Top_left][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Top_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ezp_Top_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZpffcoupl[gHIGGS_Vp_Top_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZpffcoupl[gHIGGS_Vp_Top_right][0] = coupl[0];
                    D.selfDZpffcoupl[gHIGGS_Vp_Top_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghwwp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwwpcoupl[gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwwpcoupl[gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHwwpcoupl[gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ghwpwp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDHwpwpcoupl[gHIGGS_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDHwpwpcoupl[gHIGGS_VV_1][0] = coupl[0];
                    D.selfDHwpwpcoupl[gHIGGS_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_El_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_El_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_El_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_El_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_El_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_El_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_El_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_El_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Mu_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Mu_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Mu_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Mu_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Mu_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Mu_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Mu_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Ta_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Ta_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Ta_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Ta_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Ta_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Ta_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Ta_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Up_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Up_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Up_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Up_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Up_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Up_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Up_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Up_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Chm_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Chm_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Chm_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Chm_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Chm_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Chm_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Chm_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Top_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Top_left], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Top_left][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Top_left][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "ewp_Top_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDWpffcoupl[gHIGGS_Vp_Top_right], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDWpffcoupl[gHIGGS_Vp_Top_right][0] = coupl[0];
                    D.selfDWpffcoupl[gHIGGS_Vp_Top_right][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "zprime_qq_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZqqcoupl[gZPRIME_QQ_LEFT], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZqqcoupl[gZPRIME_QQ_LEFT][0] = coupl[0];
                    D.selfDZqqcoupl[gZPRIME_QQ_LEFT][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "zprime_qq_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZqqcoupl[gZPRIME_QQ_RIGHT], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZqqcoupl[gZPRIME_QQ_RIGHT][0] = coupl[0];
                    D.selfDZqqcoupl[gZPRIME_QQ_RIGHT][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "zprime_zz_1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZvvcoupl[gZPRIME_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZvvcoupl[gZPRIME_VV_1][0] = coupl[0];
                    D.selfDZvvcoupl[gZPRIME_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "zprime_zz_2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDZvvcoupl[gZPRIME_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDZvvcoupl[gZPRIME_VV_2][0] = coupl[0];
                    D.selfDZvvcoupl[gZPRIME_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "graviton_qq_left", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGqqcoupl[gGRAVITON_QQ_LEFT], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGqqcoupl[gGRAVITON_QQ_LEFT][0] = coupl[0];
                    D.selfDGqqcoupl[gGRAVITON_QQ_LEFT][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "graviton_qq_right", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGqqcoupl[gGRAVITON_QQ_RIGHT], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGqqcoupl[gGRAVITON_QQ_RIGHT][0] = coupl[0];
                    D.selfDGqqcoupl[gGRAVITON_QQ_RIGHT][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "a1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGggcoupl[gGRAVITON_GG_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGggcoupl[gGRAVITON_GG_1][0] = coupl[0];
                    D.selfDGggcoupl[gGRAVITON_GG_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "a2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGggcoupl[gGRAVITON_GG_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGggcoupl[gGRAVITON_GG_2][0] = coupl[0];
                    D.selfDGggcoupl[gGRAVITON_GG_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "a3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGggcoupl[gGRAVITON_GG_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGggcoupl[gGRAVITON_GG_3][0] = coupl[0];
                    D.selfDGggcoupl[gGRAVITON_GG_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "a4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGggcoupl[gGRAVITON_GG_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGggcoupl[gGRAVITON_GG_4][0] = coupl[0];
                    D.selfDGggcoupl[gGRAVITON_GG_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "a5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGggcoupl[gGRAVITON_GG_5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGggcoupl[gGRAVITON_GG_5][0] = coupl[0];
                    D.selfDGggcoupl[gGRAVITON_GG_5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_1][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_2][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_3][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_4][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_5][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_6][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_7][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_8][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b9", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_9], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_9][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_9][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "b10", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_VV_10], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_VV_10][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_VV_10][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzgs1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_ZA_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_ZA_1][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_ZA_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_ZA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_ZA_2][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_ZA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_ZA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_ZA_3][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_ZA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_ZA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_ZA_4][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_ZA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzgs8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_ZA_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_ZA_8][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_ZA_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bgsgs1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_AA_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_AA_1][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_AA_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bgsgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_AA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_AA_2][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_AA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bgsgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_AA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_AA_3][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_AA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bgsgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_AA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_AA_4][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_AA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bgsgs8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvcoupl[gGRAVITON_AA_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvcoupl[gGRAVITON_AA_8][0] = coupl[0];
                    D.selfDGvvcoupl[gGRAVITON_AA_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_1][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_2][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_3][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_4][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_5][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_6][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_7][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_8][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp9", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_9], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_9][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_9][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzzp10", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_VV_10], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_VV_10][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_VV_10][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpgs1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_ZA_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_ZA_1][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_ZA_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpgs2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_ZA_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_ZA_2][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_ZA_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpgs3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_ZA_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_ZA_3][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_ZA_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpgs4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_ZA_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_ZA_4][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_ZA_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpgs8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvvpcoupl[gGRAVITON_ZA_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvvpcoupl[gGRAVITON_ZA_8][0] = coupl[0];
                    D.selfDGvvpcoupl[gGRAVITON_ZA_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp1", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_1], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_1][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_1][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp2", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_2], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_2][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_2][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp3", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_3], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_3][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_3][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp4", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_4], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_4][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_4][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp5", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_5], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_5][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_5][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp6", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_6], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_6][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_6][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp7", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_7], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_7][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_7][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp8", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_8], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_8][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_8][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp9", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_9], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_9][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_9][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "bzpzp10", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDGvpvpcoupl[gGRAVITON_VV_10], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDGvpvpcoupl[gGRAVITON_VV_10][0] = coupl[0];
                    D.selfDGvpvpcoupl[gGRAVITON_VV_10][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dV_A", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dVA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dVA][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dVA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dP_A", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dPA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dPA][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dPA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dM_A", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dMA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dMA][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dMA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dFour_A", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dFourA], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dFourA][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dFourA][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dV_Z", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dVZ], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dVZ][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dVZ][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dP_Z", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dPZ], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dPZ][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dPZ][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dM_Z", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dMZ], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dMZ][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dMZ][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dFour_Z", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dFourZ], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dFourZ][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dFourZ][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dAAWpWm", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dAAWpWm], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dAAWpWm][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dAAWpWm][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dZAWpWm", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dZAWpWm], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dZAWpWm][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dZAWpWm][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "dZZWpWm", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDaTQGCcoupl[gATQGC_dZZWpWm], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDaTQGCcoupl[gATQGC_dZZWpWm][0] = coupl[0];
                    D.selfDaTQGCcoupl[gATQGC_dZZWpWm][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "reZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZllRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZllRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZllRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "leZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZllLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZllLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZllLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "rquZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZuuRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZuuRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZuuRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "lquZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZuuLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZuuLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZuuLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "rqdZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZddRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZddRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZddRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "lqdZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZddLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZddLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZddLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "rnZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZnunuRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZnunuRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZnunuRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "lnZ", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_ZnunuLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_ZnunuLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_ZnunuLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cranou", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_uZRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_uZRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_uZRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "clanou", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_uZLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_uZLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_uZLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "cranod", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_dZRH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_dZRH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_dZRH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        )

        .def_property(
            "clanod", 
            py::cpp_function(
                [](py::object &obj){
                    Mela &D = obj.cast<Mela&>();
                    return py::array_t<double>(std::vector<int>{2}, (const double*) &D.selfDAZffcoupl[gAZff_dZLH], obj);
                }, py::keep_alive<0, 1>()),
            py::cpp_function(
                [](Mela &D, std::array<double, 2> coupl){
                    D.selfDAZffcoupl[gAZff_dZLH][0] = coupl[0];
                    D.selfDAZffcoupl[gAZff_dZLH][1] = coupl[1];
                }, py::keep_alive<0, 1>())
        );



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
    
    py::enum_<TVar::Process>(m, "Process")
        .value("HSMHiggs",TVar::HSMHiggs) 
        .value("H0_g1prime2",TVar::H0_g1prime2)
        .value("H0hplus",TVar::H0hplus)
        .value("H0minus",TVar::H0minus)
        .value("H0_Zgsg1prime2",TVar::H0_Zgsg1prime2)
        .value("H0_Zgs",TVar::H0_Zgs)
        .value("H0_Zgs_PS",TVar::H0_Zgs_PS)
        .value("H0_gsgs",TVar::H0_gsgs)
        .value("H0_gsgs_PS",TVar::H0_gsgs_PS)
        .value("D_g1g1prime2",TVar::D_g1g1prime2)
        .value("D_g1g2",TVar::D_g1g2)
        .value("D_g1g2_pi_2",TVar::D_g1g2_pi_2)
        .value("D_g1g4",TVar::D_g1g4)
        .value("D_g1g4_pi_2",TVar::D_g1g4_pi_2)
        .value("D_zzzg",TVar::D_zzzg)
        .value("D_zzgg",TVar::D_zzgg)
        .value("D_zzzg_PS",TVar::D_zzzg_PS)
        .value("D_zzgg_PS",TVar::D_zzgg_PS)
        .value("D_zzzg_g1prime2",TVar::D_zzzg_g1prime2)
        .value("D_zzzg_g1prime2_pi_2",TVar::D_zzzg_g1prime2_pi_2)
        .value("H1minus",TVar::H1minus) 
        .value("H1plus",TVar::H1plus) 
        .value("H2_g1",TVar::H2_g1) 
        .value("H2_g2",TVar::H2_g2) 
        .value("H2_g3",TVar::H2_g3) 
        .value("H2_g4",TVar::H2_g4) 
        .value("H2_g5",TVar::H2_g5) 
        .value("H2_g1g5",TVar::H2_g1g5) 
        .value("H2_g6",TVar::H2_g6) 
        .value("H2_g7",TVar::H2_g7) 
        .value("H2_g8",TVar::H2_g8) 
        .value("H2_g9",TVar::H2_g9) 
        .value("H2_g10",TVar::H2_g10) 
        .value("bkgGammaGamma",TVar::bkgGammaGamma) 
        .value("bkgZGamma",TVar::bkgZGamma) 
        .value("bkgZJets",TVar::bkgZJets) 
        .value("bkgZZ",TVar::bkgZZ) 
        .value("bkgWW",TVar::bkgWW) 
        .value("bkgWWZZ",TVar::bkgWWZZ) 
        .value("bkgZZ_SMHiggs",TVar::bkgZZ_SMHiggs) 
        .value("bkgWW_SMHiggs",TVar::bkgWW_SMHiggs) 
        .value("bkgWWZZ_SMHiggs",TVar::bkgWWZZ_SMHiggs) 
        .value("HSMHiggs_WWZZ",TVar::HSMHiggs_WWZZ) 
        .value("D_gg10",TVar::D_gg10)
        .value("SelfDefine_spin0",TVar::SelfDefine_spin0)
        .value("SelfDefine_spin1",TVar::SelfDefine_spin1)
        .value("SelfDefine_spin2",TVar::SelfDefine_spin2)
        .value("nProcesses",TVar::nProcesses);
    
    py::enum_<pymela::HGG_indices>(m, "HGG_indices")
        .value("gHIGGS_GG_2", pymela::gHIGGS_GG_2)
        .value("gHIGGS_GG_3", pymela::gHIGGS_GG_3)
        .value("gHIGGS_GG_4", pymela::gHIGGS_GG_4)
        .value("SIZE_HGG", pymela::SIZE_HGG);

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