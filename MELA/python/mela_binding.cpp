#include "Mela.h"
#include "TVar.hh"
#include "TCouplingsBase.hh"
#include "TMCFM.hh"
#include "TUtil.hh"
#include "MELACandidate.h"
#include "TLorentzVector.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/operators.h"
namespace py = pybind11;
using namespace std;

/**
 * @defgroup Pychanges Novel Python Functions
 * @brief These are additions to the MELA code that allow for the Python to operate
 * 
 * The changes made provide for a more usable/readable version of the Python bindings. 
 * They do not change any of the functionality of MELA, 
 * but simply provide either compatibility fixes or syntactic sugar.
 * For all other PyMELA related things, refer to the subpage @PyMela_page "here".
 * @{
*/

/**
 * @defgroup ReferenceGroup Python Functions that are pass-by-reference in C++
 * @brief These functions are pass-by-reference in C++ and turned into 
 * functions that return values for usage in Python
 * @{
*/

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeP
/// @param mela Mela object instance (for python function calling using obj.<function>)
/// @param useConstant Boolean for using a multiplicative constant
/// @return the result of Mela::computeP
float computeP(Mela& mela, bool useConstant=true){
    float result;
    mela.computeP(result, useConstant);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeProdP
/// @param mela Mela object instance (for python function calling using obj.<function>)
/// @param useConstant Boolean for using a multiplicative constant
/// @return the result of Mela::computeProdP
float computeProdP(Mela& mela, bool useConstant=true) {
    float result;
    mela.computeProdP(result, useConstant);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of MelThey cover a slew of uses, from converting a::computeProdDecP
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param useConstant Boolean for using a multiplicative constant 
/// @return the result of Mela::computeProdDecP
float computeProdDecP(Mela& mela, bool useConstant=true) {
    float result;
    mela.computeProdDecP(result, useConstant);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeD_CP
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param myME The Matrix element to be used
/// @param myType The process to be used in the calculation
/// @return the result of Mela::computeD_CP
float computeD_CP(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
    float result;
    mela.computeD_CP(myME, myType, result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeProdP_VH
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param includeHiggsDecay Whether you would like to include Higgs decay in your calculation
/// @param useConstant Boolean for using a multiplicative constant 
/// @return the result of Mela::computePropP_VH
float computeProdP_VH(Mela& mela, bool includeHiggsDecay, bool useConstant) {
    float result;
    mela.computeProdP_VH(result, includeHiggsDecay, useConstant);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeProdP_ttH
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param topProcess 
/// @param topDecay 
/// @param useConstant Boolean for using a multiplicative constant 
/// @return the result of Mela::computeProdP_ttH
float computeProdP_ttH(Mela& mela, int topProcess, int topDecay, bool useConstant) {
    float result;
    mela.computeProdP_ttH(result, topProcess, topDecay, useConstant);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::compute4FermionWeight
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of Mela::compute4FermionWeight
float compute4FermionWeight(Mela& mela) {
    float result;
    mela.compute4FermionWeight(result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::getXPropagator
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param scheme The Propagator scheme you would like to use
/// @return the result of Mela::getXPropagator
float getXPropagator(Mela& mela, TVar::ResonancePropagatorScheme scheme) {
    float result;
    mela.getXPropagator(scheme, result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computePM4l
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param syst 
/// @return the result of Mela::computePM4l
float computePM4l(Mela& mela, TVar::SuperMelaSyst syst) {
    float result;
    mela.computePM4l(syst, result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeD_gg
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param myME The matrix element you would like to use
/// @param myType The process you would like to use
/// @return the result of Mela::computeD_gg
float computeD_gg(Mela& mela, TVar::MatrixElement myME, TVar::Process myType) {
    float result;
    mela.computeD_gg(myME, myType, result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief the analog of Mela::getConstant
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of Mela::getConstant
float getConstant(Mela& mela) {
    float result;
    mela.getConstant(result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief the analog of Mela::computeDijetConvBW
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param useTrueBW 
/// @return the result of Mela::computeDijetConvBW
float computeDijetConvBW(Mela& mela, bool useTrueBW) {
    float result;
    mela.computeDijetConvBW(result, useTrueBW);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief the analog of MelaIO::getWeightedMEArray
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of MelaIO::getWeightedMEArray
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

/// @ingroup ReferenceGroup
/// @brief the analog of MelaIO::getUnweightedMEArray
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of MelaIO::getUnweightedMEArray
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

/// @ingroup ReferenceGroup
/// @brief the analog of MelaIO::getPartonWeights
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of MelaIO::getPartonWeights
py::tuple getPartonWeights(Mela& mela) {
    std::pair<vector<double>, vector<double>> result;
    result.first.resize(nmsq);
    result.second.resize(nmsq);
    MelaIO* melaio = mela.getIORecord();
    melaio->getPartonWeights(result.first.data(), result.second.data());
    py::tuple py_result = py::make_tuple(py::cast(result.first), py::cast(result.second));
    return py_result;
}

/// @ingroup ReferenceGroup
/// @brief the analog of Mela::getPAux
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of Mela::getPAux
float getPAux(Mela& mela) {
    float result;
    mela.getPAux(result);
    return result;
}

/// @ingroup ReferenceGroup
/// @brief the analog of Mela::computeDecayAngles
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @note The order of the values is mH, m1, m2 (usually Z1 and Z2), cos(theta1), cos(theta2), phi, cos(theta-star), phi1
/// @return the result of Mela::computeDecayAngles
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

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeVBFAngles
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of Mela::computeVBFAngles
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

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeVBFAngles_ComplexBoost
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @return the result of Mela::computeVBFAngles_ComplexBoost
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

/// @ingroup ReferenceGroup
/// @brief analog of Mela::computeVHAngles
/// @param mela Mela object instance (for python function calling using obj.<function>) 
/// @param prod the result of Mela::computeVHAngles
/// @return 
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
/** @} */

/**
 * @defgroup Constructors Functions that construct objects differently than the C++
 * @brief These functions are factories for datatypes that either need to be created differently
 * in Python bindings or make it much easier to do so
 * @{
*/

/// @ingroup Constructors
/// @brief This function intializes a single SimpleParticle_t in the Python.
/// @anchor particle_initializer
/// @param id The particle PDG ID
/// @param x Either the pX or the pT of the vector
/// @param y Either the pY or the Eta of the vector
/// @param z Either the pZ or the Phi of the Vector
/// @param e Either the Energy or the Mass of the vector
/// @param ptEtaPhi if true, interpret the vector using <pT, Eta, Phi, m>. Otherwise, use <pX, pY, pZ, E>
/// @return A TVar::Simpleparticle_t with the provided id and Lorentz 4-vector.
SimpleParticle_t particle_initializer(int id, float x, float y, float z, float e, bool ptEtaPhi=false){
    TLorentzVector vec = TLorentzVector();
    if(ptEtaPhi){
        vec.SetPtEtaPhiM(x, y, z, e);
    } else{
        vec.SetPxPyPzE(x, y, z, e);
    }
    return SimpleParticle_t(id, vec);
}

/// @ingroup Constructors
/// @anchor collection_initializer_from_column
/// @brief This function initializes a single SimpleParticleCollection_t (or a list of SimpleParticle_t) in the Python from a series of lists
/// @note This function is useful for data in columnar form! One can imagine dumping columns of vector quantities into this function.
/// @param ids a list of PDG ids
/// @param x A list of either the pX or the pT of the vectors
/// @param y A list of either the pY or the Eta of the vectors
/// @param z A list of either the pZ or the Phi of the Vectors
/// @param e A list of either the Energy or the Mass of the vectors
/// @param ptEtaPhi if true, interpret the vector using <pT, Eta, Phi, m>. Otherwise, use <pX, pY, pZ, E>
/// @return A TVar::SimpleParticleCollection_t object
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

/// @ingroup Constructors
/// @brief This function initializes a single SimpleParticleCollection_t (or a list of SimpleParticle_t) in the Python from a single list
/// @anchor collection_initializer
/// @param listOfParticles A list of SimpleParticle_t objects
/// @return A TVar::SimpleParticleCollection_t object
SimpleParticleCollection_t collection_initializer(py::list listOfParticles){

    SimpleParticleCollection_t collection = SimpleParticleCollection_t();
    for(py::handle P : listOfParticles){
        SimpleParticle_t particle = P.cast<SimpleParticle_t>();
        collection.push_back(particle);
    }
    return collection;
}
/** @} */

/**
 * @defgroup Macros Coupling definition Macros
 * @anchor py_macros
 * @brief These are C++ macros that define named couplings in the Python code
 * @{
*/

/** 
 * @ingroup macros
 * @param arrayName This is the name of the array containing the coupling
 * @param size This is the size of the array
 * @param arrType The data type of the array
 * @brief Generates the array for spin 0 values in JHUGen and JHUGen-MCFM
 * @note These coupling entries have both a real and an imaginary component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<arrayName>()[<couplingIndex>][0] = <real>
 * m.<arrayName>()[<couplingIndex>][1] = <imag>
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_ARR_SPIN_ZERO(arrayName, size, arrType)\
        .def(#arrayName, [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<arrType>(std::vector<arrType>{nSupportedHiggses, size, 2}, (const arrType*) &D.arrayName, obj); \
        })

/** 
 * @ingroup macros
 * @param arrayName This is the name of the array containing the coupling
 * @param size This is the size of the array
 * @brief Generates the array for spin 1 and spin 2 values in JHUGen
 * @note These coupling entries have both a real and an imaginary component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<arrayName>()[<couplingIndex>][0] = <real>
 * m.<arrayName>()[<couplingIndex>][1] = <imag>
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_ARR_SPIN_ONETWO(arrayName, size)\
        .def(#arrayName, [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{size, 2}, (const double*) &D.arrayName, obj); \
        })

/** 
 * @ingroup macros
 * @brief Generates the couplings for spin 0 values in JHUGen and JHUGen-MCFM
 * @param arrayName This is the name of the array containing the coupling
 * @param couplingName This is the name of the coupling to set (should be unique for every couplingIndex/higgsIndex combination)
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @param higgsIndex This is the index of which Higgs to use (MCFM supports 2 resonances)
 * @note These coupling entries have both a real and an imaginary component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = [real, imag]
 * print(m.<couplingName>) #will print a list of [real, imag]
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(arrayName, couplingName, couplingIndex, higgsIndex)\
        .def_property(\
            #couplingName,\
            py::cpp_function(\
                [](py::object &obj){\
                    Mela &D = obj.cast<Mela&>();\
                    return py::array_t<double>(std::vector<double>{2}, (const double*) &D.arrayName[higgsIndex][couplingIndex], obj);\
                }, py::keep_alive<0, 1>()),\
            py::cpp_function(\
                [](Mela &D, std::array<double, 2> coupl){\
                    D.arrayName[higgsIndex][couplingIndex][0] = coupl[0];\
                    D.arrayName[higgsIndex][couplingIndex][1] = coupl[1];\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the couplings for real-valued spin 0 values in JHUGen and JHUGen-MCFM
 * @param arrayName This is the name of the array containing the coupling
 * @param couplingName This is the name of the coupling to set (should be unique for every couplingIndex/higgsIndex combination)
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @param higgsIndex This is the index of which Higgs to use (MCFM supports 2 resonances)
 * @note These coupling entries have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = real
 * print(m.<couplingName>) #will print a real value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_REAL_SPIN_ZERO(arrayName, couplingName, couplingIndex, higgsIndex)\
        .def_property(\
            #couplingName,\
            py::cpp_function(\
                [](py::object &obj){\
                    Mela &D = obj.cast<Mela&>();\
                    return py::array_t<double>(std::vector<double>{2}, (const double*) &D.arrayName[higgsIndex][couplingIndex], obj);\
                }, py::keep_alive<0, 1>()),\
            py::cpp_function(\
                [](Mela &D, double coupl){\
                    D.arrayName[higgsIndex][couplingIndex] = coupl;\
                }, py::keep_alive<0, 1>())\
        )
/** 
 * @ingroup macros
 * @brief Generates the couplings for spin 1 and spin 2 values in JHUGen
 * @param arrayName This is the name of the array containing the coupling
 * @param couplingName This is the name of the coupling to set (should be unique for every couplingIndex/higgsIndex combination)
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @note These coupling entries have both a real and an imaginary component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = [real, imag]
 * print(m.<couplingName>) #will print a list of [real, imag]
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(arrayName, couplingName, couplingIndex)\
        .def_property(\
            #couplingName,\
            py::cpp_function(\
                [](py::object &obj){\
                    Mela &D = obj.cast<Mela&>();\
                    return py::array_t<double>(std::vector<double>{2}, (const double*) &D.arrayName[couplingIndex], obj);\
                }, py::keep_alive<0, 1>()),\
            py::cpp_function(\
                [](Mela &D, std::array<double, 2> coupl){\
                    D.arrayName[couplingIndex][0] = coupl[0];\
                    D.arrayName[couplingIndex][1] = coupl[1];\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the couplings for C Lambda values in JHUGen and JHUGen-MCFM
 * @param arrayName This is the name of the array containing the coupling
 * @param couplingName This is the name of the coupling to set (should be unique for every couplingIndex/higgsIndex combination)
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @param higgsIndex This is the index of which Higgs to use (MCFM supports 2 resonances)
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_C_LAMBDA(arrayName, couplingName, couplingIndex, higgsIndex)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<int>(std::vector<int>{nSupportedHiggses,SIZE_HVV_CQSQ}, (const int*) &D.arrayName, obj);\
                    return array_val.at(higgsIndex,couplingIndex);\
                }),\
            py::cpp_function(\
                [](py::object &obj, double coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<int>(std::vector<int>{nSupportedHiggses,SIZE_HVV_CQSQ}, (const int*) &D.arrayName, obj);\
                    array_val.mutable_at(higgsIndex,couplingIndex) = coupl;\
                }, py::keep_alive<0, 1>())\
        )


/** 
 * @ingroup macros
 * @brief Generates the couplings for Lambda values in JHUGen and JHUGen-MCFM
 * @param arrayName This is the name of the array containing the coupling
 * @param couplingName This is the name of the coupling to set (should be unique for every couplingIndex/higgsIndex combination)
 * @param couplingIndex_1 This is the index at which the name corresponds to the array
 * @param couplingIndex_2 This is the index at which the name corresponds to the array
 * @param higgsIndex This is the index of which Higgs to use (MCFM supports 2 resonances)
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_LAMBDA(arrayName, couplingName, couplingIndex_1, couplingIndex_2, higgsIndex)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{nSupportedHiggses,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.arrayName, obj);\
                    return array_val.at(higgsIndex,couplingIndex_1, couplingIndex_2);\
                }),\
            py::cpp_function(\
                [](py::object &obj, double coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{nSupportedHiggses,SIZE_HVV_LAMBDAQSQ,SIZE_HVV_CQSQ}, (const double*) &D.arrayName, obj);\
                    array_val.mutable_at(higgsIndex,couplingIndex_1, couplingIndex_2) = coupl;\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the couplings for SMEFTSim Wilson Coefficients in MadMELA
 * @param couplingName This is the name of the Wilson Coefficient to set
 * @param couplingIndex_1 This is the index at which the name corresponds to the array
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_MADMELA(couplingName, couplingIndex_1)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{SIZE_SMEFT}, (const double*) &D.selfDSMEFTSimcoupl, obj);\
                    return array_val.at(couplingIndex_1);\
                }),\
            py::cpp_function(\
                [](py::object &obj, double coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{SIZE_SMEFT}, (const double*) &D.selfDSMEFTSimcoupl, obj);\
                    array_val.mutable_at(couplingIndex_1) = coupl;\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the couplings for gluon fusion Higgs self-coupling
 * @param couplingName This is the name of the self-coupling coefficient to set
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_HHH(couplingName, couplingIndex)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{SIZE_HHH}, (const double*) &D.selfDHHHcoupl, obj);\
                    return array_val.at(couplingIndex);\
                }),\
            py::cpp_function(\
                [](py::object &obj, double coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{SIZE_HHH}, (const double*) &D.selfDHHHcoupl, obj);\
                    array_val.mutable_at(couplingIndex) = coupl;\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the form factors for the Higgs
 * @param couplingName This is the name of the form factor to set
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_LAMBDA_FF(couplingName, couplingIndex, higgsIndex)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{nSupportedHiggses,SIZE_HVV_LAMBDAFF}, (const double*) &D.selfDHvvLambda_ff, obj);\
                    return array_val.at(higgsIndex,couplingIndex);\
                }),\
            py::cpp_function(\
                [](py::object &obj, double coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<double>(std::vector<double>{nSupportedHiggses,SIZE_HVV_LAMBDAFF}, (const double*) &D.selfDHvvLambda_ff, obj);\
                    array_val.mutable_at(higgsIndex,couplingIndex) = coupl;\
                }, py::keep_alive<0, 1>())\
        )

/** 
 * @ingroup macros
 * @brief Generates the form factors for the Higgs
 * @param couplingName This is the name of the form factor to set
 * @param couplingIndex This is the index at which the name corresponds to the array
 * @note These coupling entries only have a real component. You can set them via:
 * ~~~~~~~~~~~~~{.py}
 * import Mela
 * m = Mela.Mela()
 * m.<couplingName> = value
 * print(m.<couplingName>) #will print the value
 * ~~~~~~~~~~~~~
 */
#define MAKE_COUPLING_N_FF(couplingName, couplingIndex, higgsIndex)\
        .def_property(\
            #couplingName, \
            py::cpp_function(\
                [](py::object &obj){\
                    Mela& D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<int>(std::vector<int>{nSupportedHiggses,SIZE_HVV_NFF}, (const int*) &D.selfDHvvn_ff, obj);\
                    return array_val.at(higgsIndex,couplingIndex);\
                }),\
            py::cpp_function(\
                [](py::object &obj, int coupl){\
                    Mela &D = obj.cast<Mela&>();\
                    py::array_t array_val = py::array_t<int>(std::vector<int>{nSupportedHiggses,SIZE_HVV_NFF}, (const int*) &D.selfDHvvn_ff, obj);\
                    array_val.mutable_at(higgsIndex,couplingIndex) = coupl;\
                }, py::keep_alive<0, 1>())\
        )

PYBIND11_MAKE_OPAQUE(SimpleParticle_t)
PYBIND11_MAKE_OPAQUE(SimpleParticleCollection_t)
/// @brief The actual binding code for MELA
/// @anchor PyMELA
PYBIND11_MODULE(Mela, m) {
    py::class_<SimpleParticle_t>(m, "SimpleParticle_t")
        .def(py::init(&particle_initializer), py::arg("id"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("e"), py::arg("ptEtaPhi") = false)
        .def_property("id", 
        [](SimpleParticle_t& P){
            return P.first;
        },
        [](SimpleParticle_t& P, int id){
            P.first = id;
        })
        .def_property_readonly("PxPyPzE_vector", [](SimpleParticle_t& P){
            return py::make_tuple(P.second.Px(), P.second.Py(), P.second.Pz(), P.second.E());
        })
        .def_property_readonly("PtEtaPhiM_vector", [](SimpleParticle_t& P){
            return py::make_tuple(P.second.Pt(), P.second.Eta(), P.second.Phi(), P.second.M());
        })
        .def("setVector", [](SimpleParticle_t& P, double Px, double Py, double Pz, double E){
            P.second = TLorentzVector(Px, Py, Pz, E);
        })
        .def("__repr__",[](SimpleParticle_t& P){
            return "SimpleParticle(id=" + std::to_string(P.first) + ",P4=<" + std::to_string(P.second.Px()) + ", " + std::to_string(P.second.Py()) + ", " + std::to_string(P.second.Pz()) + ", " + std::to_string(P.second.E()) + ">)";
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
        })
        .def("__getitem__", [](SimpleParticleCollection_t &C, int idx){
            return &(C.at(idx));
        })
        .def("__setitem__", [](SimpleParticleCollection_t &C, int idx, SimpleParticle_t &P){
            C.at(idx) = P;
        })
        .def("Sum", [](SimpleParticleCollection_t &C){
            TLorentzVector sum = TLorentzVector(0,0,0,0);
            for(SimpleParticle_t P : C){
                sum += P.second;
            }
            return py::make_tuple(sum.Px(), sum.Py(), sum.Pz(), sum.E());
        })
        .def("MTotal", [](SimpleParticleCollection_t &C){
            TLorentzVector sum = TLorentzVector(0,0,0,0);
            for(SimpleParticle_t P : C){
                sum += P.second;
            }
            return sum.M();
        })
        .def(py::pickle(
            [](const SimpleParticleCollection_t& C){
                py::list pickleable;
                for(int i = 0; i < (int)C.size(); i++){
                    pickleable.append(py::make_tuple(C[i].first, C[i].second.Px(), C[i].second.Py(), C[i].second.Pz(), C[i].second.E()));
                }
                return py::cast<py::tuple>(pickleable);
            },
            [](py::tuple t){
                std::vector<int> ids;
                std::vector<double> x;
                std::vector<double> y;
                std::vector<double> z;
                std::vector<double> e;
                // for(auto it = t.begin(); it != t.end(); it++){
                for(int i = 0; i < (int)t.size(); i++){
                    py::tuple o = t[i];
                    ids.push_back(o[0].cast<int>());
                    x.push_back(o[1].cast<double>());
                    y.push_back(o[2].cast<double>());
                    z.push_back(o[3].cast<double>());
                    e.push_back(o[4].cast<double>());
                }
                return collection_initializer_from_column(ids, x, y, z, e);
            }
        ))
        .def("__repr__", [](SimpleParticleCollection_t& C){
            return "SimpleParticleCollection_t(length=" + to_string(C.size()) + ")";
        });
    
    py::class_<MELAParticle>(m, "MELAParticle")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init([](int id_, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return MELAParticle(id_, vec);
        }))
        .def(py::init<MELAParticle const &>())
        .def("assign", &MELAParticle::operator=)
        .def("__iadd__", [](MELAParticle& mp, MELAParticle* part){
            mp += part;
        })
        .def("addVec", [](MELAParticle& mp, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            mp += vec;
        })
        .def_readwrite("id", &MELAParticle::id)
        .def_property_readonly("p4", 
        [](MELAParticle& mp){
            TLorentzVector p4 = mp.p4;
            return py::make_tuple(p4.Px(), p4.Py(), p4.Pz(), p4.E());
        })
        .def("setP4", [](MELAParticle& mp, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            mp.p4 = vec;
        })
        .def_readwrite("passSelection", &MELAParticle::passSelection)
        .def_readwrite("genStatus", &MELAParticle::genStatus)
        .def_readwrite("lifetime", &MELAParticle::lifetime)

        .def("swap", &MELAParticle::swap)
        .def("setSelected", &MELAParticle::setSelected)
        .def("setGenStatus", &MELAParticle::setGenStatus)
        .def("setLifetime", &MELAParticle::setLifetime)
        .def("addMother", &MELAParticle::addMother)
        .def("addDaughter", &MELAParticle::addDaughter)
        .def("getNMothers", &MELAParticle::getNMothers)
        .def("getNDaughters", &MELAParticle::getNDaughters)
        .def("getDaughterIds", &MELAParticle::getDaughterIds)
        .def("getMother", &MELAParticle::getMother)
        .def("getDaughter", &MELAParticle::getDaughter)
        .def("getRelatedParticles", &MELAParticle::getRelatedParticles)
        .def("getDaughterParticles", &MELAParticle::getDaughterParticles)

        .def("charge", &MELAParticle::charge)
        .def("m", &MELAParticle::m)
        .def("x", &MELAParticle::x)
        .def("y", &MELAParticle::y)
        .def("z", &MELAParticle::z)
        .def("t", &MELAParticle::t)
        .def("p", &MELAParticle::p)
        .def("pt", &MELAParticle::pt)
        .def("eta", &MELAParticle::eta)
        .def("phi", &MELAParticle::phi)
        .def("rapidity", &MELAParticle::rapidity)
        .def("dot", [](MELAParticle& m, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return m.dot(vec);
        })
        .def("dot", [](MELAParticle& m, MELAParticle& other){
            return m.dot(other);
        })
        .def("euclidean_dot", [](MELAParticle& m, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return m.euclidean_dot(vec);
        })
        .def("euclidean_dot", [](MELAParticle& m, MELAParticle& other){
            return m.dot(other);
        })
        .def("deltaR", [](MELAParticle& m, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return m.deltaR(vec);
        })
        .def("deltaR", [](MELAParticle& m, MELAParticle& other){
            return m.deltaR(other);
        })
        .def("boost", [](MELAParticle& m, double Px, double Py, double Pz, bool boostAll){
            TVector3 v = TVector3(Px, Py, Pz);
            m.boost(v, boostAll);
            return m;
        })
        .def("vect", [](MELAParticle& m){
            TVector3 v = m.vect();
            return py::make_tuple(v.X(), v.Y(), v.Z());
        })
        .def("calculateTotalDisplacement", [](MELAParticle& m){
            TVector3 v = m.calculateTotalDisplacement();
            return py::make_tuple(v.X(), v.Y(), v.Z());
        })

        .def("__repr__", [](MELAParticle &mp){
            return "MELAParticle(id=" + to_string(mp.id) + ",passSelection=" + to_string(mp.passSelection) + ",genStatus=" + to_string(mp.genStatus) + ",lifetime=" + to_string(mp.lifetime) + ")";
        });

    py::class_<MELAThreeBodyDecayCandidate, MELAParticle>(m, "MELAThreeBodyDecayCandidate")
        .def(py::init<>())
        .def(py::init([](int id_, double Px, double Py, double Pz, double E){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return MELAThreeBodyDecayCandidate(id_, vec);
        }))
        .def(py::init<MELAParticle*, MELAParticle*, MELAParticle*>())
        .def(py::init<MELAThreeBodyDecayCandidate const &>())
        .def("assign", &MELAThreeBodyDecayCandidate::operator=)
        .def("swap", &MELAThreeBodyDecayCandidate::swap)
        .def("setPartnerParticle", &MELAThreeBodyDecayCandidate::setPartnerParticle)
        .def("setWFermion", &MELAThreeBodyDecayCandidate::setWFermion)
        .def("setWAntifermion", &MELAThreeBodyDecayCandidate::setWAntifermion)
        .def("getPartnerParticle", [](const MELAThreeBodyDecayCandidate& mt){
            return mt.getPartnerParticle();
            })
        .def("getWFermion", [](const MELAThreeBodyDecayCandidate& mt){
            return mt.getWFermion();
            })
        .def("getWAntifermion", [](const MELAThreeBodyDecayCandidate& mt){
            return mt.getWAntifermion();
            })
        .def("testPreSelectedDaughters", &MELAThreeBodyDecayCandidate::testPreSelectedDaughters)
        .def("getWmass", &MELAThreeBodyDecayCandidate::getWmass)
        .def("checkCandidateExists", &MELAThreeBodyDecayCandidate::checkCandidateExists);

    py::class_<MELACandidate, MELAParticle>(m, "MELACandidate")
        .def(py::init<>())
        .def(py::init<int, bool>(), py::arg("id_"), py::arg("associatedByHighestPt_")=false)
        .def(py::init([](int id_, double Px, double Py, double Pz, double E, bool associatedByHighestPt_){
            TLorentzVector vec = TLorentzVector();
            vec.SetPxPyPzE(Px, Py, Pz, E);
            return MELACandidate(id_, vec, associatedByHighestPt_);
        }), py::arg("id_"), py::arg("Px"), py::arg("Py"), py::arg("Pz"), py::arg("E"), py::arg("associatedByHighestPt_")=false)
        .def(py::init<MELACandidate const &>())
        .def("assign", &MELACandidate::operator=)
        .def("shallowCopy", &MELACandidate::shallowCopy)
        .def("swap", &MELACandidate::swap)
        .def("getSortedDaughter", &MELACandidate::getSortedDaughter)
        .def("getSortedV", &MELACandidate::getSortedV)
        .def("getAssociatedLepton", &MELACandidate::getAssociatedLepton)
        .def("getAssociatedNeutrino", &MELACandidate::getAssociatedNeutrino)
        .def("getAssociatedPhoton", &MELACandidate::getAssociatedPhoton)
        .def("getAssociatedJet", &MELACandidate::getAssociatedJet)
        .def("getAssociatedTop", &MELACandidate::getAssociatedTop)
        .def("getSortedDaughters", [](const MELACandidate& mc){
            return mc.getSortedDaughters();
            })
        .def("getSortedVs", [](const MELACandidate& mc){
            return mc.getSortedVs();
            })
        .def("getAssociatedLeptons", [](const MELACandidate& mc){
            return mc.getAssociatedLeptons();
            })
        .def("getAssociatedNeutrinos", [](const MELACandidate& mc){
            return mc.getAssociatedNeutrinos();
            })
        .def("getAssociatedPhotons", [](const MELACandidate& mc){
            return mc.getAssociatedPhotons();
            })
        .def("getAssociatedJets", [](const MELACandidate& mc){
            return mc.getAssociatedJets();
            })
        .def("getAssociatedTops", [](const MELACandidate& mc){
            return mc.getAssociatedTops();
            })
        .def("getAssociatedSortedVs", [](const MELACandidate& mc){
            return mc.getAssociatedSortedVs();
            })
        .def("getRelatedParticles", &MELACandidate::getRelatedParticles)
        .def("getDaughterParticles", &MELACandidate::getDaughterParticles)
        .def("getNAssociatedLeptons", &MELACandidate::getNAssociatedLeptons)
        .def("getNAssociatedNeutrinos", &MELACandidate::getNAssociatedNeutrinos)
        .def("getNAssociatedPhotons", &MELACandidate::getNAssociatedPhotons)
        .def("getNAssociatedJets", &MELACandidate::getNAssociatedJets)
        .def("getNAssociatedTops", &MELACandidate::getNAssociatedTops)
        .def("getNSortedVs", &MELACandidate::getNSortedVs)
        .def("addAssociatedLepton", &MELACandidate::addAssociatedLepton)
        .def("addAssociatedNeutrino", &MELACandidate::addAssociatedNeutrino)
        .def("addAssociatedPhoton", &MELACandidate::addAssociatedPhoton)
        .def("addAssociatedJet", &MELACandidate::addAssociatedJet)
        .def("addAssociatedTop", &MELACandidate::addAssociatedTop)
        .def("addSortedV", &MELACandidate::addSortedV)
        .def("addAssociatedVs", &MELACandidate::addAssociatedVs)
        .def("resetVs", &MELACandidate::resetVs)
        .def("recreateVs", &MELACandidate::recreateVs)
        .def("sortDaughters", &MELACandidate::sortDaughters)
        .def("testPreSelectedDaughters", &MELACandidate::testPreSelectedDaughters)
        .def("testShallowCopy", &MELACandidate::testShallowCopy)
        .def("daughtersInterfere", &MELACandidate::daughtersInterfere)
        .def("setDecayMode", &MELACandidate::setDecayMode)
        .def("setAddAssociatedByHighestPt", &MELACandidate::setAddAssociatedByHighestPt)
        .def("setShallowCopy", &MELACandidate::setShallowCopy)
        .def("addUnordered", [](MELACandidate& mc, MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
            mc.addUnordered(myParticle, particleArray);
            })
        .def("addUnordered", [](MELACandidate& mc, MELAThreeBodyDecayCandidate* myParticle, std::vector<MELAThreeBodyDecayCandidate*>& particleArray){
            mc.addUnordered(myParticle, particleArray);
            })
        .def("addByHighestPt", [](MELACandidate& mc, MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
            mc.addByHighestPt(myParticle, particleArray);
            })
        .def("addByHighestPt", [](MELACandidate& mc, MELAThreeBodyDecayCandidate* myParticle, std::vector<MELAThreeBodyDecayCandidate*>& particleArray){
            mc.addByHighestPt(myParticle, particleArray);
            });

    py::class_<TVar::event_scales_type>(m, "event_scales_type")
        .def(py::init<TVar::EventScaleScheme, TVar::EventScaleScheme, double, double>())
        .def_readwrite("renomalizationScheme", &TVar::event_scales_type::renomalizationScheme)
        .def_readwrite("factorizationScheme", &TVar::event_scales_type::factorizationScheme)
        .def_readwrite("ren_scale_factor", &TVar::event_scales_type::ren_scale_factor)
        .def_readwrite("fac_scale_factor", &TVar::event_scales_type::fac_scale_factor);

    py::class_<TVar::simple_event_record>(m, "simple_event_record")
        .def(py::init<>())
        .def_readwrite("AssociationCode", &TVar::simple_event_record::AssociationCode)
        .def_readwrite("AssociationVCompatibility", &TVar::simple_event_record::AssociationVCompatibility)
        .def_readwrite("nRequested_AssociatedJets", &TVar::simple_event_record::nRequested_AssociatedJets)
        .def_readwrite("nRequested_AssociatedLeptons", &TVar::simple_event_record::nRequested_AssociatedLeptons)
        .def_readwrite("nRequested_AssociatedPhotons", &TVar::simple_event_record::nRequested_AssociatedPhotons)
        .def_readwrite("nRequested_Tops", &TVar::simple_event_record::nRequested_Tops)
        .def_readwrite("nRequested_Antitops", &TVar::simple_event_record::nRequested_Antitops)

        .def_readwrite("intermediateVid", &TVar::simple_event_record::intermediateVid)
        .def_readwrite("pDaughters", &TVar::simple_event_record::pDaughters)
        .def_readwrite("pAssociated", &TVar::simple_event_record::pAssociated)
        .def_readwrite("pMothers", &TVar::simple_event_record::pMothers)
        .def_readwrite("pTopDaughters", &TVar::simple_event_record::pTopDaughters)
        .def_readwrite("pAntitopDaughters", &TVar::simple_event_record::pAntitopDaughters)
        .def_readwrite("pStableTops", &TVar::simple_event_record::pStableTops)
        .def_readwrite("pStableAntitops", &TVar::simple_event_record::pStableAntitops);

    m.def("PrintCandidateSummary", [](TVar::simple_event_record curCand){
        TVar::simple_event_record* curRecord = &curCand;
        TUtil::PrintCandidateSummary(curRecord);
    });
    m.def("PrintCandidateSummary", [](MELACandidate curCand){
        MELACandidate* curRecord = &curCand;
        TUtil::PrintCandidateSummary(curRecord);
    });
    py::class_<Mela>(m, "Mela")
        .def(py::init<double, double, TVar::VerbosityLevel>())
        .def(py::init<double, double>())
        .def(py::init<double>())
        .def(py::init<>())
        .def(py::init<Mela const &>())
        .def("__repr__", [](Mela& me){
            return "Mela";
        })
        .def("setProcess", &Mela::setProcess)
        .def("setVerbosity", &Mela::setVerbosity)
        .def("setInputEvent", &Mela::setInputEvent, py::arg("pDaughters"), py::arg("pAssociated")=nullptr, py::arg("pMothers")=nullptr, py::arg("isGen")=false)
        .def("setCandidateDecayMode", &Mela::setCandidateDecayMode)
        .def("setMelaHiggsMass", &Mela::setMelaHiggsMass, py::arg("myHiggsMass"), py::arg("index")=0)
        .def("setMelaHiggsWidth", &Mela::setMelaHiggsWidth, py::arg("myHiggsWidth")=-1, py::arg("index")=0)
        .def("setMelaHiggsMassWidth", &Mela::setMelaHiggsMassWidth, py::arg("myHiggsMass"), py::arg("myHiggsWidth"), py::arg("index"))
        .def("setMelaLeptonInterference", &Mela::setMelaLeptonInterference)
        .def("setRenFacScaleMode", &Mela::setRenFacScaleMode)
        .def("SetMadgraphCKMElements", &Mela::SetMadgraphCKMElements, py::arg("ckmlambda")=0.2265, py::arg("ckma")=0.79, py::arg("ckmrho")=0.141, py::arg("ckmeta")=0.357, py::arg("force_refresh")=false)

        .def("resetInputEvent", &Mela::resetInputEvent)
        .def("resetMass", &Mela::resetMass)
        .def("resetYukawaMass", &Mela::resetYukawaMass)
        .def("resetWidth", &Mela::resetWidth)
        .def("resetQuarkMasses", &Mela::resetQuarkMasses)
        .def("resetMCFM_EWKParameters", &Mela::resetMCFM_EWKParameters)

        .def("getPrimaryMass", &Mela::getPrimaryMass)
        .def("getPrimaryWidth", &Mela::getPrimaryWidth)
        .def("getHiggsWidthAtPoleMass", &Mela::getHiggsWidthAtPoleMass)
        .def("GetMadgraphCKMElement", &Mela::GetMadgraphCKMElement, py::arg("iquark"), py::arg("jquark"))

        .def("getIORecord", &Mela::getIORecord)
        .def("getWeightedMEArray", &getWeightedMEArray)
        .def("getUnweightedMEArray", &getUnweightedMEArray)
        .def("getPartonWeights", &getPartonWeights)
        .def("getPAux", &getPAux)
        .def("getRenFacScaleMode", &Mela::getRenFacScaleMode)


        .def("computeP", &computeP, py::arg("useConstant"))
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
        .def("PrintCurrentCandidateSummary", [](Mela& mela){
            MELACandidate* curCand = mela.getCurrentCandidate();
            TUtil::PrintCandidateSummary(curCand);
        })
        .def("cleanLinkedFiles", &Mela::cleanLinkedFiles)
        .def("calculate4Momentum", &Mela::calculate4Momentum)

        .def("computeDecayAngles", &computeDecayAngles)
        .def("computeVBFAngles", &computeVBFAngles)
        .def("computeVBFAngles_ComplexBoost", &computeVBFAngles_ComplexBoost)
        .def("computeVHAngles", &computeVHAngles)

        .def_readwrite("differentiate_HWW_HZZ", &Mela::differentiate_HWW_HZZ)

        //Raw coupling arrays
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHggcoupl,SIZE_HGG,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHg4g4coupl,SIZE_HGG,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHqqcoupl,SIZE_HQQ,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHbbcoupl,SIZE_HQQ,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHttcoupl,SIZE_HQQ,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHb4b4coupl,SIZE_HQQ,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHt4t4coupl,SIZE_HQQ,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHzzcoupl,SIZE_HVV,double)
        MAKE_COUPLING_ARR_SPIN_ZERO(selfDHwwcoupl,SIZE_HVV,double)
        .def("selfDHHHcoupl", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{SIZE_HHH}, (const double*) &D.selfDHHHcoupl, obj); \
        })
        .def("selfDHzzLambda_qsq", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{nSupportedHiggses, SIZE_HVV_LAMBDAQSQ, SIZE_HVV_CQSQ}, (const double*) &D.selfDHzzLambda_qsq, obj); \
        })
        .def("selfDHwwLambda_qsq", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{nSupportedHiggses, SIZE_HVV_LAMBDAQSQ, SIZE_HVV_CQSQ}, (const double*) &D.selfDHwwLambda_qsq, obj); \
        })
        .def("selfDHzzCLambda_qsq", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<int>(std::vector<int>{nSupportedHiggses, SIZE_HVV_CQSQ}, (const int*) &D.selfDHzzCLambda_qsq, obj); \
        })
        .def("selfDHwwCLambda_qsq", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<int>(std::vector<int>{nSupportedHiggses, SIZE_HVV_CQSQ}, (const int*) &D.selfDHwwCLambda_qsq, obj); \
        })
        .def("selfDHvvLambda_ff", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{nSupportedHiggses, SIZE_HVV_LAMBDAFF}, (const double*) &D.selfDHvvLambda_ff, obj); \
        })
        .def("selfDHvvn_ff", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<int>(std::vector<int>{nSupportedHiggses, SIZE_HVV_NFF}, (const int*) &D.selfDHvvn_ff, obj); \
        })
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDHzzpcoupl,SIZE_HVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDHzpzpcoupl,SIZE_HVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDZpffcoupl,SIZE_Vpff)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDHwwpcoupl,SIZE_HVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDHwpwpcoupl,SIZE_HVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDWpffcoupl,SIZE_Vpff)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDZqqcoupl,SIZE_ZQQ)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDZvvcoupl,SIZE_ZVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDGqqcoupl,SIZE_GQQ)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDGggcoupl,SIZE_GGG)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDGvvcoupl,SIZE_GVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDGvvpcoupl,SIZE_GVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDGvpvpcoupl,SIZE_GVV)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDaTQGCcoupl,SIZE_ATQGC)
        MAKE_COUPLING_ARR_SPIN_ONETWO(selfDAZffcoupl,SIZE_AZff)
        .def("selfDSMEFTSimcoupl", [](py::object &obj){ \
            Mela &D = obj.cast<Mela&>(); \
            return py::array_t<double>(std::vector<double>{SIZE_SMEFT}, (const double*) &D.selfDSMEFTSimcoupl, obj); \
        })
        .def_readwrite("M_Zprime", &Mela::selfDM_Zprime)
        .def_readwrite("Ga_Zprime", &Mela::selfDGa_Zprime)
        .def_readwrite("M_Wprime", &Mela::selfDM_Wprime)
        .def_readwrite("Ga_Wprime", &Mela::selfDGa_Wprime)


        //Here be couplings
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, ghg2, gHIGGS_GG_2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, gh2g2, gHIGGS_GG_2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, ghg3, gHIGGS_GG_3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, gh2g3, gHIGGS_GG_3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, ghg4, gHIGGS_GG_4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHggcoupl, gh2g4, gHIGGS_GG_4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1, gHIGGS_VV_1, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1, gHIGGS_VV_1, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2, gHIGGS_VV_2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2, gHIGGS_VV_2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3, gHIGGS_VV_3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3, gHIGGS_VV_3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4, gHIGGS_VV_4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4, gHIGGS_VV_4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghzgs2, gHIGGS_ZA_2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2zgs2, gHIGGS_ZA_2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghzgs3, gHIGGS_ZA_3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2zgs3, gHIGGS_ZA_3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghzgs4, gHIGGS_ZA_4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2zgs4, gHIGGS_ZA_4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghgsgs2, gHIGGS_AA_2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2gsgs2, gHIGGS_AA_2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghgsgs3, gHIGGS_AA_3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2gsgs3, gHIGGS_AA_3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghgsgs4, gHIGGS_AA_4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2gsgs4, gHIGGS_AA_4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime, gHIGGS_VV_1_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime, gHIGGS_VV_1_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime2, gHIGGS_VV_1_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime2, gHIGGS_VV_1_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime3, gHIGGS_VV_1_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime3, gHIGGS_VV_1_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime4, gHIGGS_VV_1_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime4, gHIGGS_VV_1_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime5, gHIGGS_VV_1_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime5, gHIGGS_VV_1_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime, gHIGGS_VV_2_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime, gHIGGS_VV_2_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime2, gHIGGS_VV_2_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime2, gHIGGS_VV_2_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime3, gHIGGS_VV_2_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime3, gHIGGS_VV_2_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime4, gHIGGS_VV_2_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime4, gHIGGS_VV_2_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime5, gHIGGS_VV_2_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime5, gHIGGS_VV_2_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime, gHIGGS_VV_3_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime, gHIGGS_VV_3_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime2, gHIGGS_VV_3_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime2, gHIGGS_VV_3_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime3, gHIGGS_VV_3_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime3, gHIGGS_VV_3_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime4, gHIGGS_VV_3_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime4, gHIGGS_VV_3_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime5, gHIGGS_VV_3_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime5, gHIGGS_VV_3_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime, gHIGGS_VV_4_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime, gHIGGS_VV_4_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime2, gHIGGS_VV_4_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime2, gHIGGS_VV_4_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime3, gHIGGS_VV_4_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime3, gHIGGS_VV_4_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime4, gHIGGS_VV_4_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime4, gHIGGS_VV_4_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime5, gHIGGS_VV_4_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime5, gHIGGS_VV_4_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghzgs1_prime2, gHIGGS_ZA_1_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2zgs1_prime2, gHIGGS_ZA_1_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime6, gHIGGS_VV_1_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime6, gHIGGS_VV_1_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz1_prime7, gHIGGS_VV_1_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z1_prime7, gHIGGS_VV_1_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime6, gHIGGS_VV_2_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime6, gHIGGS_VV_2_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz2_prime7, gHIGGS_VV_2_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z2_prime7, gHIGGS_VV_2_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime6, gHIGGS_VV_3_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime6, gHIGGS_VV_3_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz3_prime7, gHIGGS_VV_3_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z3_prime7, gHIGGS_VV_3_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime6, gHIGGS_VV_3_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime6, gHIGGS_VV_3_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, ghz4_prime7, gHIGGS_VV_3_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHzzcoupl, gh2z4_prime7, gHIGGS_VV_3_PRIME7, 1)

        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q1sq,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q1sq_h2,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z11,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z11_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z12,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z12_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z13,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z13_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z14,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z14_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q2sq,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q2sq_h2,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z21,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z21_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z22,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z22_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z23,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z23_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z24,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z24_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q12sq,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHzzCLambda_qsq, cz_q12sq_h2,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z01,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z01_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z02,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z02_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z03,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z03_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z04,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHzzLambda_qsq, Lambda_z04_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1, gHIGGS_VV_1, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1, gHIGGS_VV_1, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2, gHIGGS_VV_2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2, gHIGGS_VV_2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3, gHIGGS_VV_3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3, gHIGGS_VV_3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4, gHIGGS_VV_4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4, gHIGGS_VV_4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime, gHIGGS_VV_1_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime, gHIGGS_VV_1_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime2, gHIGGS_VV_1_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime2, gHIGGS_VV_1_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime3, gHIGGS_VV_1_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime3, gHIGGS_VV_1_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime4, gHIGGS_VV_1_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime4, gHIGGS_VV_1_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime5, gHIGGS_VV_1_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime5, gHIGGS_VV_1_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime, gHIGGS_VV_2_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime, gHIGGS_VV_2_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime2, gHIGGS_VV_2_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime2, gHIGGS_VV_2_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime3, gHIGGS_VV_2_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime3, gHIGGS_VV_2_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime4, gHIGGS_VV_2_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime4, gHIGGS_VV_2_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime5, gHIGGS_VV_2_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime5, gHIGGS_VV_2_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime, gHIGGS_VV_3_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime, gHIGGS_VV_3_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime2, gHIGGS_VV_3_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime2, gHIGGS_VV_3_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime3, gHIGGS_VV_3_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime3, gHIGGS_VV_3_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime4, gHIGGS_VV_3_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime4, gHIGGS_VV_3_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime5, gHIGGS_VV_3_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime5, gHIGGS_VV_3_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime, gHIGGS_VV_4_PRIME, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime, gHIGGS_VV_4_PRIME, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime2, gHIGGS_VV_4_PRIME2, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime2, gHIGGS_VV_4_PRIME2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime3, gHIGGS_VV_4_PRIME3, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime3, gHIGGS_VV_4_PRIME3, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime4, gHIGGS_VV_4_PRIME4, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime4, gHIGGS_VV_4_PRIME4, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime5, gHIGGS_VV_4_PRIME5, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime5, gHIGGS_VV_4_PRIME5, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime6, gHIGGS_VV_1_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime6, gHIGGS_VV_1_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw1_prime7, gHIGGS_VV_1_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w1_prime7, gHIGGS_VV_1_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime6, gHIGGS_VV_2_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime6, gHIGGS_VV_2_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw2_prime7, gHIGGS_VV_2_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w2_prime7, gHIGGS_VV_2_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime6, gHIGGS_VV_3_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime6, gHIGGS_VV_3_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw3_prime7, gHIGGS_VV_3_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w3_prime7, gHIGGS_VV_3_PRIME7, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime6, gHIGGS_VV_3_PRIME6, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime6, gHIGGS_VV_3_PRIME6, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, ghw4_prime7, gHIGGS_VV_3_PRIME7, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHwwcoupl, gh2w4_prime7, gHIGGS_VV_3_PRIME7, 1)

        MAKE_COUPLING_HHH(c6, gHIGGS_HH_c6)
        MAKE_COUPLING_HHH(t1, gHIGGS_HH_t1)
        MAKE_COUPLING_HHH(t2, gHIGGS_HH_t2)
        MAKE_COUPLING_HHH(t3, gHIGGS_HH_t3)
        MAKE_COUPLING_HHH(t4, gHIGGS_HH_t4)
        MAKE_COUPLING_HHH(t5, gHIGGS_HH_t5)
        MAKE_COUPLING_HHH(t6, gHIGGS_HH_t6)
        MAKE_COUPLING_HHH(t7, gHIGGS_HH_t7)
        MAKE_COUPLING_HHH(w1, gHIGGS_HH_w1)
        MAKE_COUPLING_HHH(w2, gHIGGS_HH_w2)
        MAKE_COUPLING_HHH(w3, gHIGGS_HH_w3)
        MAKE_COUPLING_HHH(w4, gHIGGS_HH_w4)
        MAKE_COUPLING_HHH(w5, gHIGGS_HH_w5)

        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q1sq,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q1sq_h2,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w11,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w11_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w12,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w12_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w13,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w13_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w14,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ1, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w14_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ1, 1)

        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q2sq,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q2sq_h2,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w21,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w21_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w22,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w22_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w23,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w23_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w24,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ2, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w24_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ2, 1)

        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q12sq,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_C_LAMBDA(selfDHwwCLambda_qsq, cw_q12sq_h2,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w01,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w01_h2,  LambdaHIGGS_QSQ_VV_1,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w02,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w02_h2,  LambdaHIGGS_QSQ_VV_2,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w03,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w03_h2,  LambdaHIGGS_QSQ_VV_3,  cLambdaHIGGS_VV_QSQ12, 1)

        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w04,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ12, 0)
        MAKE_COUPLING_LAMBDA(selfDHwwLambda_qsq, Lambda_w04_h2,  LambdaHIGGS_QSQ_VV_4,  cLambdaHIGGS_VV_QSQ12, 1)
        
        MAKE_COUPLING_N_FF(n_ff1,  nHIGGS_VV_FF1, 0)
        MAKE_COUPLING_N_FF(n2_ff1,  nHIGGS_VV_FF1, 1)
        
        MAKE_COUPLING_N_FF(n_ff2,  nHIGGS_VV_FF2, 0)
        MAKE_COUPLING_N_FF(n2_ff2,  nHIGGS_VV_FF2, 1)

        MAKE_COUPLING_LAMBDA_FF(Lambda_ff1, LambdaHIGGS_VV_FF1, 0)
        MAKE_COUPLING_LAMBDA_FF(Lambda2_ff1, LambdaHIGGS_VV_FF1, 1)

        MAKE_COUPLING_LAMBDA_FF(Lambda_ff2, LambdaHIGGS_VV_FF2, 0)
        MAKE_COUPLING_LAMBDA_FF(Lambda2_ff2, LambdaHIGGS_VV_FF2, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHqqcoupl, kappa, gHIGGS_KAPPA, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHqqcoupl, kappa_h2, gHIGGS_KAPPA, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHqqcoupl, kappa_tilde, gHIGGS_KAPPA_TILDE, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHqqcoupl, kappa_tilde_h2, gHIGGS_KAPPA_TILDE, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHttcoupl, kappa_top, gHIGGS_KAPPA, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHttcoupl, kappa_top_h2, gHIGGS_KAPPA, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHttcoupl, kappa_top_tilde, gHIGGS_KAPPA_TILDE, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHttcoupl, kappa_top_tilde_h2, gHIGGS_KAPPA_TILDE, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHbbcoupl, kappa_bot, gHIGGS_KAPPA, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHbbcoupl, kappa_bot_h2, gHIGGS_KAPPA, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHbbcoupl, kappa_bot_tilde, gHIGGS_KAPPA_TILDE, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHbbcoupl, kappa_bot_tilde_h2, gHIGGS_KAPPA_TILDE, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHt4t4coupl, kappa_4gen_top, gHIGGS_KAPPA, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHt4t4coupl, kappa_4gen_top_h2, gHIGGS_KAPPA, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHt4t4coupl, kappa_4gen_top_tilde, gHIGGS_KAPPA_TILDE, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHt4t4coupl, kappa_4gen_top_tilde_h2, gHIGGS_KAPPA_TILDE, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHb4b4coupl, kappa_4gen_bot, gHIGGS_KAPPA, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHb4b4coupl, kappa_4gen_bot_h2, gHIGGS_KAPPA, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHb4b4coupl, kappa_4gen_bot_tilde, gHIGGS_KAPPA_TILDE, 0)
        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO(selfDHb4b4coupl, kappa_4gen_bot_tilde_h2, gHIGGS_KAPPA_TILDE, 1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1, gHIGGS_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2, gHIGGS_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3, gHIGGS_VV_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4, gHIGGS_VV_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzpgs2, gHIGGS_ZA_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzpgs3, gHIGGS_ZA_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzpgs4, gHIGGS_ZA_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime, gHIGGS_VV_1_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime2, gHIGGS_VV_1_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime3, gHIGGS_VV_1_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime4, gHIGGS_VV_1_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime5, gHIGGS_VV_1_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime, gHIGGS_VV_2_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime2, gHIGGS_VV_2_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime3, gHIGGS_VV_2_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime4, gHIGGS_VV_2_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime5, gHIGGS_VV_2_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime, gHIGGS_VV_3_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime2, gHIGGS_VV_3_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime3, gHIGGS_VV_3_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime4, gHIGGS_VV_3_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime5, gHIGGS_VV_3_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime, gHIGGS_VV_4_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime2, gHIGGS_VV_4_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime3, gHIGGS_VV_4_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime4, gHIGGS_VV_4_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime5, gHIGGS_VV_4_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzpgs1_prime2, gHIGGS_ZA_1_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime6, gHIGGS_VV_1_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp1_prime7, gHIGGS_VV_1_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime6, gHIGGS_VV_2_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp2_prime7, gHIGGS_VV_2_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime6, gHIGGS_VV_3_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp3_prime7, gHIGGS_VV_3_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime6, gHIGGS_VV_3_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzzpcoupl, ghzzp4_prime7, gHIGGS_VV_3_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1, gHIGGS_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2, gHIGGS_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3, gHIGGS_VV_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4, gHIGGS_VV_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime, gHIGGS_VV_1_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime2, gHIGGS_VV_1_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime3, gHIGGS_VV_1_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime4, gHIGGS_VV_1_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime5, gHIGGS_VV_1_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime, gHIGGS_VV_2_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime2, gHIGGS_VV_2_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime3, gHIGGS_VV_2_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime4, gHIGGS_VV_2_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime5, gHIGGS_VV_2_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime, gHIGGS_VV_3_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime2, gHIGGS_VV_3_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime3, gHIGGS_VV_3_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime4, gHIGGS_VV_3_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime5, gHIGGS_VV_3_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime, gHIGGS_VV_4_PRIME)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime2, gHIGGS_VV_4_PRIME2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime3, gHIGGS_VV_4_PRIME3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime4, gHIGGS_VV_4_PRIME4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime5, gHIGGS_VV_4_PRIME5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime6, gHIGGS_VV_1_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp1_prime7, gHIGGS_VV_1_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime6, gHIGGS_VV_2_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp2_prime7, gHIGGS_VV_2_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime6, gHIGGS_VV_3_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp3_prime7, gHIGGS_VV_3_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime6, gHIGGS_VV_3_PRIME6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHzpzpcoupl, ghzpzp4_prime7, gHIGGS_VV_3_PRIME7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_El_left, gHIGGS_Vp_El_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_El_right, gHIGGS_Vp_El_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Mu_left, gHIGGS_Vp_Mu_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Mu_right, gHIGGS_Vp_Mu_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Ta_left, gHIGGS_Vp_Ta_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Ta_right, gHIGGS_Vp_Ta_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_NuE_left, gHIGGS_Vp_NuE_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_NuE_right, gHIGGS_Vp_NuE_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Dn_left, gHIGGS_Vp_Dn_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Dn_right, gHIGGS_Vp_Dn_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Up_left, gHIGGS_Vp_Up_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Up_right, gHIGGS_Vp_Up_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Str_left, gHIGGS_Vp_Str_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Str_right, gHIGGS_Vp_Str_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Chm_left, gHIGGS_Vp_Chm_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Chm_right, gHIGGS_Vp_Chm_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Bot_left, gHIGGS_Vp_Bot_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Bot_right, gHIGGS_Vp_Bot_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Top_left, gHIGGS_Vp_Top_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZpffcoupl, ezp_Top_right, gHIGGS_Vp_Top_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHwwpcoupl, ghwwp1, gHIGGS_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDHwpwpcoupl, ghwpwp1, gHIGGS_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_El_left, gHIGGS_Vp_El_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_El_right, gHIGGS_Vp_El_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Mu_left, gHIGGS_Vp_Mu_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Mu_right, gHIGGS_Vp_Mu_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Ta_left, gHIGGS_Vp_Ta_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Ta_right, gHIGGS_Vp_Ta_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Up_left, gHIGGS_Vp_Up_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Up_right, gHIGGS_Vp_Up_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Chm_left, gHIGGS_Vp_Chm_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Chm_right, gHIGGS_Vp_Chm_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Top_left, gHIGGS_Vp_Top_left)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDWpffcoupl, ewp_Top_right, gHIGGS_Vp_Top_right)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZqqcoupl, zprime_qq_left, gZPRIME_QQ_LEFT)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZqqcoupl, zprime_qq_right, gZPRIME_QQ_RIGHT)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZvvcoupl, zprime_zz_1, gZPRIME_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDZvvcoupl, zprime_zz_2, gZPRIME_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGqqcoupl, graviton_qq_left, gGRAVITON_QQ_LEFT)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGqqcoupl, graviton_qq_right, gGRAVITON_QQ_RIGHT)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGggcoupl, a1, gGRAVITON_GG_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGggcoupl, a2, gGRAVITON_GG_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGggcoupl, a3, gGRAVITON_GG_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGggcoupl, a4, gGRAVITON_GG_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGggcoupl, a5, gGRAVITON_GG_5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b1, gGRAVITON_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b2, gGRAVITON_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b3, gGRAVITON_VV_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b4, gGRAVITON_VV_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b5, gGRAVITON_VV_5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b6, gGRAVITON_VV_6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b7, gGRAVITON_VV_7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b8, gGRAVITON_VV_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b9, gGRAVITON_VV_9)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, b10, gGRAVITON_VV_10)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bzgs1, gGRAVITON_ZA_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bzgs2, gGRAVITON_ZA_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bzgs3, gGRAVITON_ZA_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bzgs4, gGRAVITON_ZA_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bzgs8, gGRAVITON_ZA_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bgsgs1, gGRAVITON_AA_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bgsgs2, gGRAVITON_AA_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bgsgs3, gGRAVITON_AA_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bgsgs4, gGRAVITON_AA_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvcoupl, bgsgs8, gGRAVITON_AA_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp1, gGRAVITON_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp2, gGRAVITON_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp3, gGRAVITON_VV_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp4, gGRAVITON_VV_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp5, gGRAVITON_VV_5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp6, gGRAVITON_VV_6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp7, gGRAVITON_VV_7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp8, gGRAVITON_VV_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp9, gGRAVITON_VV_9)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzzp10, gGRAVITON_VV_10)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzpgs1, gGRAVITON_ZA_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzpgs2, gGRAVITON_ZA_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzpgs3, gGRAVITON_ZA_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzpgs4, gGRAVITON_ZA_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvvpcoupl, bzpgs8, gGRAVITON_ZA_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp1, gGRAVITON_VV_1)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp2, gGRAVITON_VV_2)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp3, gGRAVITON_VV_3)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp4, gGRAVITON_VV_4)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp5, gGRAVITON_VV_5)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp6, gGRAVITON_VV_6)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp7, gGRAVITON_VV_7)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp8, gGRAVITON_VV_8)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp9, gGRAVITON_VV_9)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDGvpvpcoupl, bzpzp10, gGRAVITON_VV_10)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dV_A, gATQGC_dVA)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dP_A, gATQGC_dPA)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dM_A, gATQGC_dMA)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dFour_A, gATQGC_dFourA)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dV_Z, gATQGC_dVZ)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dP_Z, gATQGC_dPZ)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dM_Z, gATQGC_dMZ)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dFour_Z, gATQGC_dFourZ)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dAAWpWm, gATQGC_dAAWpWm)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dZAWpWm, gATQGC_dZAWpWm)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDaTQGCcoupl, dZZWpWm, gATQGC_dZZWpWm)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, reZ, gAZff_ZllRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, leZ, gAZff_ZllLH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, rquZ, gAZff_ZuuRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, lquZ, gAZff_ZuuLH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, rqdZ, gAZff_ZddRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, lqdZ, gAZff_ZddLH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, rnZ, gAZff_ZnunuRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, lnZ, gAZff_ZnunuLH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, cranou, gAZff_uZRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, clanou, gAZff_uZLH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, cranod, gAZff_dZRH)

        MAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO(selfDAZffcoupl, clanod, gAZff_dZLH)

        MAKE_COUPLING_MADMELA(mdl_ch, gMDL_ch)
        MAKE_COUPLING_MADMELA(mdl_chbox, gMDL_chbox)
        MAKE_COUPLING_MADMELA(mdl_chdd, gMDL_chdd)
        MAKE_COUPLING_MADMELA(mdl_chg, gMDL_chg)
        MAKE_COUPLING_MADMELA(mdl_chw, gMDL_chw)
        MAKE_COUPLING_MADMELA(mdl_chb, gMDL_chb)
        MAKE_COUPLING_MADMELA(mdl_chwb, gMDL_chwb)
        MAKE_COUPLING_MADMELA(mdl_cehre, gMDL_cehre)
        MAKE_COUPLING_MADMELA(mdl_cuhre, gMDL_cuhre)
        MAKE_COUPLING_MADMELA(mdl_cdhre, gMDL_cdhre)
        MAKE_COUPLING_MADMELA(mdl_cewre, gMDL_cewre)
        MAKE_COUPLING_MADMELA(mdl_cebre, gMDL_cebre)
        MAKE_COUPLING_MADMELA(mdl_cugre, gMDL_cugre)
        MAKE_COUPLING_MADMELA(mdl_cuwre, gMDL_cuwre)
        MAKE_COUPLING_MADMELA(mdl_cubre, gMDL_cubre)
        MAKE_COUPLING_MADMELA(mdl_cdgre, gMDL_cdgre)
        MAKE_COUPLING_MADMELA(mdl_cdwre, gMDL_cdwre)
        MAKE_COUPLING_MADMELA(mdl_cdbre, gMDL_cdbre)
        MAKE_COUPLING_MADMELA(mdl_chl1, gMDL_chl1)
        MAKE_COUPLING_MADMELA(mdl_chl3, gMDL_chl3)
        MAKE_COUPLING_MADMELA(mdl_che, gMDL_che)
        MAKE_COUPLING_MADMELA(mdl_chq1, gMDL_chq1)
        MAKE_COUPLING_MADMELA(mdl_chq3, gMDL_chq3)
        MAKE_COUPLING_MADMELA(mdl_chu, gMDL_chu)
        MAKE_COUPLING_MADMELA(mdl_chd, gMDL_chd)
        MAKE_COUPLING_MADMELA(mdl_chudre, gMDL_chudre)
        MAKE_COUPLING_MADMELA(mdl_cll, gMDL_cll)
        MAKE_COUPLING_MADMELA(mdl_cll1, gMDL_cll1)
        MAKE_COUPLING_MADMELA(mdl_cqq1, gMDL_cqq1)
        MAKE_COUPLING_MADMELA(mdl_cqq11, gMDL_cqq11)
        MAKE_COUPLING_MADMELA(mdl_cqq3, gMDL_cqq3)
        MAKE_COUPLING_MADMELA(mdl_cqq31, gMDL_cqq31)
        MAKE_COUPLING_MADMELA(mdl_clq1, gMDL_clq1)
        MAKE_COUPLING_MADMELA(mdl_clq3, gMDL_clq3)
        MAKE_COUPLING_MADMELA(mdl_cee, gMDL_cee)
        MAKE_COUPLING_MADMELA(mdl_cuu, gMDL_cuu)
        MAKE_COUPLING_MADMELA(mdl_cuu1, gMDL_cuu1)
        MAKE_COUPLING_MADMELA(mdl_cdd, gMDL_cdd)
        MAKE_COUPLING_MADMELA(mdl_cdd1, gMDL_cdd1)
        MAKE_COUPLING_MADMELA(mdl_ceu, gMDL_ceu)
        MAKE_COUPLING_MADMELA(mdl_ced, gMDL_ced)
        MAKE_COUPLING_MADMELA(mdl_cud1, gMDL_cud1)
        MAKE_COUPLING_MADMELA(mdl_cud8, gMDL_cud8)
        MAKE_COUPLING_MADMELA(mdl_cle, gMDL_cle)
        MAKE_COUPLING_MADMELA(mdl_clu, gMDL_clu)
        MAKE_COUPLING_MADMELA(mdl_cld, gMDL_cld)
        MAKE_COUPLING_MADMELA(mdl_cqe, gMDL_cqe)
        MAKE_COUPLING_MADMELA(mdl_cqu1, gMDL_cqu1)
        MAKE_COUPLING_MADMELA(mdl_cqu8, gMDL_cqu8)
        MAKE_COUPLING_MADMELA(mdl_cqd1, gMDL_cqd1)
        MAKE_COUPLING_MADMELA(mdl_cqd8, gMDL_cqd8)
        MAKE_COUPLING_MADMELA(mdl_cledqre, gMDL_cledqre)
        MAKE_COUPLING_MADMELA(mdl_cquqd1re, gMDL_cquqd1re)
        MAKE_COUPLING_MADMELA(mdl_cquqd11re, gMDL_cquqd11re)
        MAKE_COUPLING_MADMELA(mdl_cquqd8re, gMDL_cquqd8re)
        MAKE_COUPLING_MADMELA(mdl_cquqd81re, gMDL_cquqd81re)
        MAKE_COUPLING_MADMELA(mdl_clequ1re, gMDL_clequ1re)
        MAKE_COUPLING_MADMELA(mdl_clequ3re, gMDL_clequ3re)
        MAKE_COUPLING_MADMELA(mdl_cgtil, gMDL_cgtil)
        MAKE_COUPLING_MADMELA(mdl_cwtil, gMDL_cwtil)
        MAKE_COUPLING_MADMELA(mdl_chgtil, gMDL_chgtil)
        MAKE_COUPLING_MADMELA(mdl_chwtil, gMDL_chwtil)
        MAKE_COUPLING_MADMELA(mdl_chbtil, gMDL_chbtil)
        MAKE_COUPLING_MADMELA(mdl_chwbtil, gMDL_chwbtil)
        MAKE_COUPLING_MADMELA(mdl_cewim, gMDL_cewim)
        MAKE_COUPLING_MADMELA(mdl_cebim, gMDL_cebim)
        MAKE_COUPLING_MADMELA(mdl_cugim, gMDL_cugim)
        MAKE_COUPLING_MADMELA(mdl_cuwim, gMDL_cuwim)
        MAKE_COUPLING_MADMELA(mdl_cubim, gMDL_cubim)
        MAKE_COUPLING_MADMELA(mdl_cdgim, gMDL_cdgim)
        MAKE_COUPLING_MADMELA(mdl_cdwim, gMDL_cdwim)
        MAKE_COUPLING_MADMELA(mdl_cdbim, gMDL_cdbim)
        MAKE_COUPLING_MADMELA(mdl_chudim, gMDL_chudim)
        MAKE_COUPLING_MADMELA(mdl_cehim, gMDL_cehim)
        MAKE_COUPLING_MADMELA(mdl_cuhim, gMDL_cuhim)
        MAKE_COUPLING_MADMELA(mdl_cdhim, gMDL_cdhim)
        MAKE_COUPLING_MADMELA(mdl_cledqim, gMDL_cledqim)
        MAKE_COUPLING_MADMELA(mdl_cquqd1im, gMDL_cquqd1im)
        MAKE_COUPLING_MADMELA(mdl_cquqd8im, gMDL_cquqd8im)
        MAKE_COUPLING_MADMELA(mdl_cquqd11im, gMDL_cquqd11im)
        MAKE_COUPLING_MADMELA(mdl_cquqd81im, gMDL_cquqd81im)
        MAKE_COUPLING_MADMELA(mdl_clequ1im, gMDL_clequ1im)
        MAKE_COUPLING_MADMELA(mdl_clequ3im, gMDL_clequ3im);

    py::enum_<TVar::VerbosityLevel>(m, "VerbosityLevel", py::arithmetic())
        .value("SILENT",TVar::SILENT)
        .value("ERROR",TVar::ERROR)
        .value("INFO",TVar::INFO)
        .value("DEBUG",TVar::DEBUG)
        .value("DEBUG_VERBOSE",TVar::DEBUG_VERBOSE)
        .value("DEBUG_MECHECK",TVar::DEBUG_MECHECK);

    py::enum_<TVar::MatrixElement>(m, "MatrixElement")
        .value("MCFM",TVar::MCFM)
        .value("JHUGen",TVar::JHUGen)
        .value("ANALYTICAL",TVar::ANALYTICAL)
        .value("MADGRAPH",TVar::MADGRAPH);

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

    py::enum_<TVar::ResonancePropagatorScheme>(m, "ResonancePropagatorScheme")
        .value("NoPropagator", TVar::NoPropagator)
        .value("RunningWidth", TVar::RunningWidth)
        .value("FixedWidth", TVar::FixedWidth)
        .value("CPS", TVar::CPS)
        .value("AltRunningWidth", TVar::AltRunningWidth);

    py::enum_<TVar::EventScaleScheme>(m, "EventScaleScheme")
        .value("DefaultScaleScheme", TVar::DefaultScaleScheme)
        .value("Fixed_mH", TVar::Fixed_mH)
        .value("Fixed_mW", TVar::Fixed_mW)
        .value("Fixed_mZ", TVar::Fixed_mZ)
        .value("Fixed_mWPlusmH", TVar::Fixed_mWPlusmH)
        .value("Fixed_mZPlusmH", TVar::Fixed_mZPlusmH)
        .value("Fixed_TwomtPlusmH", TVar::Fixed_TwomtPlusmH)
        .value("Fixed_mtPlusmH", TVar::Fixed_mtPlusmH)
        .value("Dynamic_qH", TVar::Dynamic_qH)
        .value("Dynamic_qJJH", TVar::Dynamic_qJJH)
        .value("Dynamic_qJJ_qH", TVar::Dynamic_qJJ_qH)
        .value("Dynamic_qJ_qJ_qH", TVar::Dynamic_qJ_qJ_qH)
        .value("Dynamic_HT", TVar::Dynamic_HT)
        .value("Dynamic_Leading_pTJ", TVar::Dynamic_Leading_pTJ)
        .value("Dynamic_Softest_pTJ", TVar::Dynamic_Softest_pTJ)
        .value("Dynamic_RandomUniform_Constrained", TVar::Dynamic_RandomUniform_Constrained)
        .value("nEventScaleSchemes", TVar::nEventScaleSchemes);

    py::enum_<TVar::CandidateDecayMode>(m, "CandidateDecayMode")
        .value("CandidateDecay_Stable", TVar::CandidateDecay_Stable)
        .value("CandidateDecay_ff", TVar::CandidateDecay_ff)
        .value("CandidateDecay_WW", TVar::CandidateDecay_WW)
        .value("CandidateDecay_ZZ", TVar::CandidateDecay_ZZ)
        .value("CandidateDecay_ZW", TVar::CandidateDecay_ZW)
        .value("CandidateDecay_ZG", TVar::CandidateDecay_ZG)
        .value("CandidateDecay_WG", TVar::CandidateDecay_WG)
        .value("CandidateDecay_GG", TVar::CandidateDecay_GG);
    
    py::enum_<TVar::SuperMelaSyst>(m, "SuperMelaSyst")
        .value("SMSyst_None", TVar::SMSyst_None)
        .value("SMSyst_ScaleUp", TVar::SMSyst_ScaleUp)
        .value("SMSyst_ScaleDown", TVar::SMSyst_ScaleDown)
        .value("SMSyst_ResUp", TVar::SMSyst_ResUp)
        .value("SMSyst_ResDown", TVar::SMSyst_ResDown);

    py::enum_<TVar::LeptonInterference>(m, "LeptonInterference")
        .value("DefaultLeptonInterf", TVar::DefaultLeptonInterf)
        .value("InterfOn", TVar::InterfOn)
        .value("InterfOff", TVar::InterfOff);

    py::enum_<TVar::FermionMassRemoval>(m, "FermionMassRemoval")
        .value("NoRemoval", TVar::NoRemoval)
        .value("ConserveDifermionMass", TVar::ConserveDifermionMass)
        .value("MomentumToEnergy", TVar::MomentumToEnergy)
        .value("nFermionMassRemovalSchemes", TVar::nFermionMassRemovalSchemes);

    py::enum_<CouplingIndex_HQQ>(m, "CouplingIndex_HQQ")
        .value("gHIGGS_KAPPA", gHIGGS_KAPPA)
        .value("gHIGGS_KAPPA_TILDE", gHIGGS_KAPPA_TILDE)
        .value("SIZE_HQQ", SIZE_HQQ);

    py::enum_<CouplingIndex_HGG>(m, "CouplingIndex_HGG")
        .value("gHIGGS_GG_2", gHIGGS_GG_2)
        .value("gHIGGS_GG_3", gHIGGS_GG_3)
        .value("gHIGGS_GG_4", gHIGGS_GG_4)
        .value("SIZE_HGG", SIZE_HGG);

    py::enum_<CouplingIndex_HVV>(m, "CouplingIndex_HVV")
        .value("gHIGGS_VV_1", gHIGGS_VV_1)
        .value("gHIGGS_VV_2", gHIGGS_VV_2)
        .value("gHIGGS_VV_3", gHIGGS_VV_3)
        .value("gHIGGS_VV_4", gHIGGS_VV_4)
        .value("gHIGGS_ZA_2", gHIGGS_ZA_2)
        .value("gHIGGS_ZA_3", gHIGGS_ZA_3)
        .value("gHIGGS_ZA_4", gHIGGS_ZA_4)
        .value("gHIGGS_AA_2", gHIGGS_AA_2)
        .value("gHIGGS_AA_3", gHIGGS_AA_3)
        .value("gHIGGS_AA_4", gHIGGS_AA_4)
        .value("gHIGGS_VV_1_PRIME", gHIGGS_VV_1_PRIME)
        .value("gHIGGS_VV_1_PRIME2", gHIGGS_VV_1_PRIME2)
        .value("gHIGGS_VV_1_PRIME3", gHIGGS_VV_1_PRIME3)
        .value("gHIGGS_VV_1_PRIME4", gHIGGS_VV_1_PRIME4)
        .value("gHIGGS_VV_1_PRIME5", gHIGGS_VV_1_PRIME5)
        .value("gHIGGS_VV_2_PRIME", gHIGGS_VV_2_PRIME)
        .value("gHIGGS_VV_2_PRIME2", gHIGGS_VV_2_PRIME2)
        .value("gHIGGS_VV_2_PRIME3", gHIGGS_VV_2_PRIME3)
        .value("gHIGGS_VV_2_PRIME4", gHIGGS_VV_2_PRIME4)
        .value("gHIGGS_VV_2_PRIME5", gHIGGS_VV_2_PRIME5)
        .value("gHIGGS_VV_3_PRIME", gHIGGS_VV_3_PRIME)
        .value("gHIGGS_VV_3_PRIME2", gHIGGS_VV_3_PRIME2)
        .value("gHIGGS_VV_3_PRIME3", gHIGGS_VV_3_PRIME3)
        .value("gHIGGS_VV_3_PRIME4", gHIGGS_VV_3_PRIME4)
        .value("gHIGGS_VV_3_PRIME5", gHIGGS_VV_3_PRIME5)
        .value("gHIGGS_VV_4_PRIME", gHIGGS_VV_4_PRIME)
        .value("gHIGGS_VV_4_PRIME2", gHIGGS_VV_4_PRIME2)
        .value("gHIGGS_VV_4_PRIME3", gHIGGS_VV_4_PRIME3)
        .value("gHIGGS_VV_4_PRIME4", gHIGGS_VV_4_PRIME4)
        .value("gHIGGS_VV_4_PRIME5", gHIGGS_VV_4_PRIME5)
        .value("gHIGGS_ZA_1_PRIME2", gHIGGS_ZA_1_PRIME2)
        .value("gHIGGS_VV_1_PRIME6", gHIGGS_VV_1_PRIME6)
        .value("gHIGGS_VV_1_PRIME7", gHIGGS_VV_1_PRIME7)
        .value("gHIGGS_VV_2_PRIME6", gHIGGS_VV_2_PRIME6)
        .value("gHIGGS_VV_2_PRIME7", gHIGGS_VV_2_PRIME7)
        .value("gHIGGS_VV_3_PRIME6", gHIGGS_VV_3_PRIME6)
        .value("gHIGGS_VV_3_PRIME7", gHIGGS_VV_3_PRIME7)
        .value("gHIGGS_VV_4_PRIME6", gHIGGS_VV_4_PRIME6)
        .value("gHIGGS_VV_4_PRIME7", gHIGGS_VV_4_PRIME7)
        .value("SIZE_HVV", SIZE_HVV);

    py::enum_<CouplingIndex_LAMBDAQSQ>(m, "CouplingIndex_LAMBDAQSQ")
        .value("LambdaHIGGS_QSQ_VV_1", LambdaHIGGS_QSQ_VV_1)
        .value("LambdaHIGGS_QSQ_VV_2", LambdaHIGGS_QSQ_VV_2)
        .value("LambdaHIGGS_QSQ_VV_3", LambdaHIGGS_QSQ_VV_3)
        .value("LambdaHIGGS_QSQ_VV_4", LambdaHIGGS_QSQ_VV_4)
        .value("SIZE_HVV_LAMBDAQSQ", SIZE_HVV_LAMBDAQSQ);

    py::enum_<CouplingIndex_HVV_CQSQ>(m, "CouplingIndex_HVV_CQSQ")
        .value("cLambdaHIGGS_VV_QSQ1", cLambdaHIGGS_VV_QSQ1)
        .value("cLambdaHIGGS_VV_QSQ2", cLambdaHIGGS_VV_QSQ2)
        .value("cLambdaHIGGS_VV_QSQ12", cLambdaHIGGS_VV_QSQ12)
        .value("SIZE_HVV_CQSQ", SIZE_HVV_CQSQ);

    py::enum_<CouplingIndex_Vpff>(m, "CouplingIndex_Vpff")
        .value("gHIGGS_Vp_El_left", gHIGGS_Vp_El_left)
        .value("gHIGGS_Vp_El_right", gHIGGS_Vp_El_right)
        .value("gHIGGS_Vp_Mu_left", gHIGGS_Vp_Mu_left)
        .value("gHIGGS_Vp_Mu_right", gHIGGS_Vp_Mu_right)
        .value("gHIGGS_Vp_Ta_left", gHIGGS_Vp_Ta_left)
        .value("gHIGGS_Vp_Ta_right", gHIGGS_Vp_Ta_right)
        .value("gHIGGS_Vp_NuE_left", gHIGGS_Vp_NuE_left)
        .value("gHIGGS_Vp_NuE_right", gHIGGS_Vp_NuE_right)
        .value("gHIGGS_Vp_Dn_left", gHIGGS_Vp_Dn_left)
        .value("gHIGGS_Vp_Dn_right", gHIGGS_Vp_Dn_right)
        .value("gHIGGS_Vp_Up_left", gHIGGS_Vp_Up_left)
        .value("gHIGGS_Vp_Up_right", gHIGGS_Vp_Up_right)
        .value("gHIGGS_Vp_Str_left", gHIGGS_Vp_Str_left)
        .value("gHIGGS_Vp_Str_right", gHIGGS_Vp_Str_right)
        .value("gHIGGS_Vp_Chm_left", gHIGGS_Vp_Chm_left)
        .value("gHIGGS_Vp_Chm_right", gHIGGS_Vp_Chm_right)
        .value("gHIGGS_Vp_Bot_left", gHIGGS_Vp_Bot_left)
        .value("gHIGGS_Vp_Bot_right", gHIGGS_Vp_Bot_right)
        .value("gHIGGS_Vp_Top_left", gHIGGS_Vp_Top_left)
        .value("gHIGGS_Vp_Top_right", gHIGGS_Vp_Top_right)
        .value("SIZE_Vpff", SIZE_Vpff);

    py::enum_<CouplingIndex_ZQQ>(m, "CouplingIndex_ZQQ")
        .value("gZPRIME_QQ_LEFT", gZPRIME_QQ_LEFT)
        .value("gZPRIME_QQ_RIGHT", gZPRIME_QQ_RIGHT)
        .value("SIZE_ZQQ", SIZE_ZQQ);

    py::enum_<CouplingIndex_ZVV>(m, "CouplingIndex_ZVV")
        .value("gZPRIME_VV_1", gZPRIME_VV_1)
        .value("gZPRIME_VV_2", gZPRIME_VV_2)
        .value("SIZE_ZVV", SIZE_ZVV);

    py::enum_<CouplingIndex_GQQ>(m, "CouplingIndex_GQQ")
        .value("gGRAVITON_QQ_LEFT", gGRAVITON_QQ_LEFT)
        .value("gGRAVITON_QQ_RIGHT", gGRAVITON_QQ_RIGHT)
        .value("SIZE_GQQ", SIZE_GQQ);

    py::enum_<CouplingIndex_GGG>(m, "CouplingIndex_GGG")
        .value("gGRAVITON_GG_1", gGRAVITON_GG_1)
        .value("gGRAVITON_GG_2", gGRAVITON_GG_2)
        .value("gGRAVITON_GG_3", gGRAVITON_GG_3)
        .value("gGRAVITON_GG_4", gGRAVITON_GG_4)
        .value("gGRAVITON_GG_5", gGRAVITON_GG_5)
        .value("SIZE_GGG", SIZE_GGG);

    py::enum_<CouplingIndex_GVV>(m, "CouplingIndex_GVV")
        .value("gGRAVITON_VV_1", gGRAVITON_VV_1)
        .value("gGRAVITON_VV_2", gGRAVITON_VV_2)
        .value("gGRAVITON_VV_3", gGRAVITON_VV_3)
        .value("gGRAVITON_VV_4", gGRAVITON_VV_4)
        .value("gGRAVITON_VV_5", gGRAVITON_VV_5)
        .value("gGRAVITON_VV_6", gGRAVITON_VV_6)
        .value("gGRAVITON_VV_7", gGRAVITON_VV_7)
        .value("gGRAVITON_VV_8", gGRAVITON_VV_8)
        .value("gGRAVITON_VV_9", gGRAVITON_VV_9)
        .value("gGRAVITON_VV_10", gGRAVITON_VV_10)
        .value("gGRAVITON_ZA_1", gGRAVITON_ZA_1)
        .value("gGRAVITON_ZA_2", gGRAVITON_ZA_2)
        .value("gGRAVITON_ZA_3", gGRAVITON_ZA_3)
        .value("gGRAVITON_ZA_4", gGRAVITON_ZA_4)
        .value("gGRAVITON_ZA_8", gGRAVITON_ZA_8)
        .value("gGRAVITON_AA_1", gGRAVITON_AA_1)
        .value("gGRAVITON_AA_2", gGRAVITON_AA_2)
        .value("gGRAVITON_AA_3", gGRAVITON_AA_3)
        .value("gGRAVITON_AA_4", gGRAVITON_AA_4)
        .value("gGRAVITON_AA_8", gGRAVITON_AA_8)
        .value("SIZE_GVV", SIZE_GVV);

    py::enum_<CouplingIndex_ATQGC>(m, "CouplingIndex_ATQGC")
        .value("gATQGC_dVA", gATQGC_dVA)
        .value("gATQGC_dPA", gATQGC_dPA)
        .value("gATQGC_dMA", gATQGC_dMA)
        .value("gATQGC_dFourA", gATQGC_dFourA)
        .value("gATQGC_dVZ", gATQGC_dVZ)
        .value("gATQGC_dPZ", gATQGC_dPZ)
        .value("gATQGC_dMZ", gATQGC_dMZ)
        .value("gATQGC_dFourZ", gATQGC_dFourZ)
        .value("gATQGC_dAAWpWm", gATQGC_dAAWpWm)
        .value("gATQGC_dZAWpWm", gATQGC_dZAWpWm)
        .value("gATQGC_dZZWpWm", gATQGC_dZZWpWm)
        .value("SIZE_ATQGC", SIZE_ATQGC);

    py::enum_<CouplingIndex_AZff>(m, "CouplingIndex_AZff")
        .value("gAZff_ZllRH", gAZff_ZllRH)
        .value("gAZff_ZllLH", gAZff_ZllLH)
        .value("gAZff_ZuuRH", gAZff_ZuuRH)
        .value("gAZff_ZuuLH", gAZff_ZuuLH)
        .value("gAZff_ZddRH", gAZff_ZddRH)
        .value("gAZff_ZddLH", gAZff_ZddLH)
        .value("gAZff_ZnunuRH", gAZff_ZnunuRH)
        .value("gAZff_ZnunuLH", gAZff_ZnunuLH)
        .value("gAZff_uZRH", gAZff_uZRH)
        .value("gAZff_uZLH", gAZff_uZLH)
        .value("gAZff_dZRH", gAZff_dZRH)
        .value("gAZff_dZLH", gAZff_dZLH)
        .value("SIZE_AZff", SIZE_AZff);

    py::enum_<CouplingIndex_SMEFT>(m, "CouplingIndex_SMEFT")
        .value("gMDL_ch", gMDL_ch)
        .value("gMDL_chbox", gMDL_chbox)
        .value("gMDL_chdd", gMDL_chdd)
        .value("gMDL_chg", gMDL_chg)
        .value("gMDL_chw", gMDL_chw)
        .value("gMDL_chb", gMDL_chb)
        .value("gMDL_chwb", gMDL_chwb)
        .value("gMDL_cehre", gMDL_cehre)
        .value("gMDL_cuhre", gMDL_cuhre)
        .value("gMDL_cdhre", gMDL_cdhre)
        .value("gMDL_cewre", gMDL_cewre)
        .value("gMDL_cebre", gMDL_cebre)
        .value("gMDL_cugre", gMDL_cugre)
        .value("gMDL_cuwre", gMDL_cuwre)
        .value("gMDL_cubre", gMDL_cubre)
        .value("gMDL_cdgre", gMDL_cdgre)
        .value("gMDL_cdwre", gMDL_cdwre)
        .value("gMDL_cdbre", gMDL_cdbre)
        .value("gMDL_chl1", gMDL_chl1)
        .value("gMDL_chl3", gMDL_chl3)
        .value("gMDL_che", gMDL_che)
        .value("gMDL_chq1", gMDL_chq1)
        .value("gMDL_chq3", gMDL_chq3)
        .value("gMDL_chu", gMDL_chu)
        .value("gMDL_chd", gMDL_chd)
        .value("gMDL_chudre", gMDL_chudre)
        .value("gMDL_cll", gMDL_cll)
        .value("gMDL_cll1", gMDL_cll1)
        .value("gMDL_cqq1", gMDL_cqq1)
        .value("gMDL_cqq11", gMDL_cqq11)
        .value("gMDL_cqq3", gMDL_cqq3)
        .value("gMDL_cqq31", gMDL_cqq31)
        .value("gMDL_clq1", gMDL_clq1)
        .value("gMDL_clq3", gMDL_clq3)
        .value("gMDL_cee", gMDL_cee)
        .value("gMDL_cuu", gMDL_cuu)
        .value("gMDL_cuu1", gMDL_cuu1)
        .value("gMDL_cdd", gMDL_cdd)
        .value("gMDL_cdd1", gMDL_cdd1)
        .value("gMDL_ceu", gMDL_ceu)
        .value("gMDL_ced", gMDL_ced)
        .value("gMDL_cud1", gMDL_cud1)
        .value("gMDL_cud8", gMDL_cud8)
        .value("gMDL_cle", gMDL_cle)
        .value("gMDL_clu", gMDL_clu)
        .value("gMDL_cld", gMDL_cld)
        .value("gMDL_cqe", gMDL_cqe)
        .value("gMDL_cqu1", gMDL_cqu1)
        .value("gMDL_cqu8", gMDL_cqu8)
        .value("gMDL_cqd1", gMDL_cqd1)
        .value("gMDL_cqd8", gMDL_cqd8)
        .value("gMDL_cledqre", gMDL_cledqre)
        .value("gMDL_cquqd1re", gMDL_cquqd1re)
        .value("gMDL_cquqd11re", gMDL_cquqd11re)
        .value("gMDL_cquqd8re", gMDL_cquqd8re)
        .value("gMDL_cquqd81re", gMDL_cquqd81re)
        .value("gMDL_clequ1re", gMDL_clequ1re)
        .value("gMDL_clequ3re", gMDL_clequ3re)
        .value("gMDL_cgtil", gMDL_cgtil)
        .value("gMDL_cwtil", gMDL_cwtil)
        .value("gMDL_chgtil", gMDL_chgtil)
        .value("gMDL_chwtil", gMDL_chwtil)
        .value("gMDL_chbtil", gMDL_chbtil)
        .value("gMDL_chwbtil", gMDL_chwbtil)
        .value("gMDL_cewim", gMDL_cewim)
        .value("gMDL_cebim", gMDL_cebim)
        .value("gMDL_cugim", gMDL_cugim)
        .value("gMDL_cuwim", gMDL_cuwim)
        .value("gMDL_cubim", gMDL_cubim)
        .value("gMDL_cdgim", gMDL_cdgim)
        .value("gMDL_cdwim", gMDL_cdwim)
        .value("gMDL_cdbim", gMDL_cdbim)
        .value("gMDL_chudim", gMDL_chudim)
        .value("gMDL_cehim", gMDL_cehim)
        .value("gMDL_cuhim", gMDL_cuhim)
        .value("gMDL_cdhim", gMDL_cdhim)
        .value("gMDL_cledqim", gMDL_cledqim)
        .value("gMDL_cquqd1im", gMDL_cquqd1im)
        .value("gMDL_cquqd8im", gMDL_cquqd8im)
        .value("gMDL_cquqd11im", gMDL_cquqd11im)
        .value("gMDL_cquqd81im", gMDL_cquqd81im)
        .value("gMDL_clequ1im", gMDL_clequ1im)
        .value("gMDL_clequ3im", gMDL_clequ3im)
        .value("SIZE_SMEFT", SIZE_SMEFT);
}
