#include <iostream>
#include <fstream>
#include "Mela.h"
#include "TVar.hh"
#include "TCouplingsBase.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
using namespace std;

TLorentzVector ptEtaPhiVector(double pt, double eta, double phi, double m){ //constructs pt, eta, phi TLorentzVector
    TLorentzVector vec = TLorentzVector();
    vec.SetPtEtaPhiM(pt, eta, phi, m);
    return vec;
}


int main(int argc, char const *argv[])
{
    // TLorentzVector g1 = TLorentzVector(0.0000000E+00,  0.0000000E+00,  0.5000000E+03, 0.5000000E+03);
    // TLorentzVector g2 = TLorentzVector(0.0000000E+00,  0.0000000E+00,  -0.5000000E+03, 0.5000000E+03);

    // SimpleParticleCollection_t* mothers = new SimpleParticleCollection_t();
    // mothers->push_back(SimpleParticle_t(21, g1));
    // mothers->push_back(SimpleParticle_t(21, g2));

    // TLorentzVector l1 = TLorentzVector(-0.2210069E+02,  0.4008035E+02, -0.7580543E+02, 0.8855133E+02);
    // TLorentzVector l2 = TLorentzVector(-0.1038496E+03, -0.3019338E+03,  0.7649492E+02, 0.3283294E+03);
    // TLorentzVector l3 = TLorentzVector(-0.1058810E+03, -0.9770964E+02,  0.4954839E+02, 0.1523581E+03);
    // TLorentzVector l4 = TLorentzVector( 0.2318313E+03,  0.3595630E+03, -0.5023788E+02, 0.4307611E+03);

    // SimpleParticleCollection_t* daughters = new SimpleParticleCollection_t();
    // daughters->push_back(SimpleParticle_t(11, l2));
    // daughters->push_back(SimpleParticle_t(-11, l1));
    // daughters->push_back(SimpleParticle_t(11, l4));
    // daughters->push_back(SimpleParticle_t(-11, l3));


    TFile* dataFile = TFile::Open("/eos/home-m/msrivast/CMSSW_14_0_0/src/JHUGenMELA/MELA/test/SM_HIGGS_JHUGEN.root");

    TTreeReader myReader("tree", dataFile);

    TTreeReaderValue<vector<short>> LHEDaughterId(myReader,   "LHEDaughterId");
    TTreeReaderValue<vector<float>> LHEDaughterPt(myReader,   "LHEDaughterPt");
    TTreeReaderValue<vector<float>> LHEDaughterEta(myReader,  "LHEDaughterEta");
    TTreeReaderValue<vector<float>> LHEDaughterPhi(myReader,  "LHEDaughterPhi");
    TTreeReaderValue<vector<float>> LHEDaughterMass(myReader, "LHEDaughterMass");
    
    TTreeReaderValue<vector<short>> LHEAssociatedParticleId(myReader,   "LHEAssociatedParticleId");
    TTreeReaderValue<vector<float>> LHEAssociatedParticlePt(myReader,   "LHEAssociatedParticlePt");
    TTreeReaderValue<vector<float>> LHEAssociatedParticleEta(myReader,  "LHEAssociatedParticleEta");
    TTreeReaderValue<vector<float>> LHEAssociatedParticlePhi(myReader,  "LHEAssociatedParticlePhi");
    TTreeReaderValue<vector<float>> LHEAssociatedParticleMass(myReader, "LHEAssociatedParticleMass");
    
    TTreeReaderValue<vector<short>> LHEMotherId(myReader, "LHEMotherId");
    TTreeReaderValue<vector<float>> LHEMotherPz(myReader, "LHEMotherPz");
    TTreeReaderValue<vector<float>> LHEMotherE(myReader,  "LHEMotherE");

    TTreeReaderValue<float> M4L(myReader, "M4L");

    Mela m = Mela(13, 125, TVar::SILENT);


    vector<float> madprobs;
    vector<float> jhugenprobs;
    vector<float> m4l;

    while(myReader.Next()){
        for(int mat_el : {1, 3}){
            SimpleParticleCollection_t* mother_collection = new SimpleParticleCollection_t();
            SimpleParticleCollection_t* daughter_collection = new SimpleParticleCollection_t();
            SimpleParticleCollection_t* associated_collection = new SimpleParticleCollection_t();

            vector<TLorentzVector> motherVecs;
            vector<TLorentzVector> daughterVecs;
            vector<TLorentzVector> associatedVecs;

            int i = 0;
            for(i = 0; i < (*LHEMotherId).size(); i++){
                mother_collection->push_back(
                    SimpleParticle_t((int)(*LHEMotherId)[i], ptEtaPhiVector(0, 0, (*LHEMotherPz)[i], (*LHEMotherE)[i]))
                );
            }

            for(i = 0; i < (*LHEDaughterId).size(); i++){
                daughter_collection->push_back(
                    SimpleParticle_t((int)(*LHEDaughterId)[i],ptEtaPhiVector((*LHEDaughterPt)[i], (*LHEDaughterEta)[i], (*LHEDaughterPhi)[i], (*LHEDaughterMass)[i]))
                );
            }

            for(i = 0; i < (*LHEAssociatedParticleId).size(); i++){
                associated_collection->push_back(
                    SimpleParticle_t((int)(*LHEAssociatedParticleId)[i],ptEtaPhiVector((*LHEAssociatedParticlePt)[i], (*LHEAssociatedParticleEta)[i], (*LHEAssociatedParticlePhi)[i], (*LHEAssociatedParticleMass)[i]))
                );
            }



            m.setProcess(TVar::SelfDefine_spin0, static_cast<TVar::MatrixElement>(mat_el), TVar::ZZINDEPENDENT);
            m.setInputEvent(daughter_collection, 0, mother_collection, true);

            if(mat_el == 1){
                m.selfDHggcoupl[0][gHIGGS_GG_2][0] = 1;
                m.selfDHzzcoupl[0][gHIGGS_VV_1][0] = 1;
            }
            float ans = 0;
            m.computeP(ans);
            m.resetInputEvent();

            if(mat_el == 1){
                jhugenprobs.push_back(ans);
            } else if(mat_el == 3){
                madprobs.push_back(ans);
            }
        }
        m4l.push_back(*M4L);
    }

    dataFile->Close();
    ofstream outputFile;
    outputFile.open("probs_output.csv");
    outputFile << "madprob, jhugenprob, m4l\n";
    for(int j = 0; j < m4l.size(); j++){
        outputFile << madprobs[j] << ", " << jhugenprobs[j] << ", " << m4l[j] << "\n";
    }
    outputFile.close();

    return 0;
}