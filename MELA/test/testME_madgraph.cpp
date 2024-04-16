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
#include "TTreeReaderArray.h"
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

    TFile* dataFile = TFile::Open("/eos/home-m/msrivast/CMSSW_14_0_0/src/JHUGenMELA/MELA/test/PS_HIGGS_MADGRAPH.root");
    // TFile* dataFile = TFile::Open("/eos/home-m/msrivast/CMSSW_14_0_0/src/HexUtils/SimulationTools/lhe_tools/HIGGS_MAD.root");


    TTreeReader myReader("tree", dataFile);

    TTreeReaderArray<float> LHEDaughterId(myReader,   "LHEDaughterId");
    TTreeReaderArray<float> LHEDaughterPt(myReader,   "LHEDaughterPt");
    TTreeReaderArray<float> LHEDaughterEta(myReader,  "LHEDaughterEta");
    TTreeReaderArray<float> LHEDaughterPhi(myReader,  "LHEDaughterPhi");
    TTreeReaderArray<float> LHEDaughterMass(myReader, "LHEDaughterMass");
    
    TTreeReaderArray<float> LHEAssociatedParticleId(myReader,   "LHEAssociatedParticleId");
    TTreeReaderArray<float> LHEAssociatedParticlePt(myReader,   "LHEAssociatedParticlePt");
    TTreeReaderArray<float> LHEAssociatedParticleEta(myReader,  "LHEAssociatedParticleEta");
    TTreeReaderArray<float> LHEAssociatedParticlePhi(myReader,  "LHEAssociatedParticlePhi");
    TTreeReaderArray<float> LHEAssociatedParticleMass(myReader, "LHEAssociatedParticleMass");
    
    TTreeReaderArray<float> LHEMotherId(myReader, "LHEMotherId");
    TTreeReaderArray<float> LHEMotherPx(myReader, "LHEMotherPx");
    TTreeReaderArray<float> LHEMotherPy(myReader,  "LHEMotherPy");
    TTreeReaderArray<float> LHEMotherPz(myReader,  "LHEMotherPz");
    TTreeReaderArray<float> LHEMotherE(myReader,  "LHEMotherE");


    TTreeReaderValue<float> M4L(myReader, "M4L");

    Mela m = Mela(13, 125, TVar::SILENT);


    vector<vector<float>> madprobs;
    madprobs.push_back(vector<float>());
    madprobs.push_back(vector<float>());

    vector<vector<float>> jhugenprobs;
    jhugenprobs.push_back(vector<float>());
    jhugenprobs.push_back(vector<float>());
    vector<float> m4l;

    int counter = 0;
    int NEVENT = 10000;
    while(myReader.Next() && counter < NEVENT){
        if((counter % 100) == 0){
            cerr << counter << endl;
        }
        SimpleParticleCollection_t* mother_collection = new SimpleParticleCollection_t();
        SimpleParticleCollection_t* daughter_collection = new SimpleParticleCollection_t();
        SimpleParticleCollection_t* associated_collection = new SimpleParticleCollection_t();

        int i = 0;
        for(i = 0; i < 2; i++){
            mother_collection->push_back(
                SimpleParticle_t((int)(LHEMotherId[i]), TLorentzVector(LHEMotherPx[i], LHEMotherPy[i], LHEMotherPz[i], LHEMotherE[i]))
            );
        }

        for(i = 0; i < 4; i++){
            daughter_collection->push_back(
                SimpleParticle_t((int)(LHEDaughterId[i]),ptEtaPhiVector(LHEDaughterPt[i], LHEDaughterEta[i], LHEDaughterPhi[i], LHEDaughterMass[i]))
            );
        }

        for(int mat_el : {1, 3}){
            for(int setup : {0,1}){
                m.setProcess(TVar::SelfDefine_spin0, static_cast<TVar::MatrixElement>(mat_el), TVar::ZZGG);
                m.setInputEvent(daughter_collection, 0, mother_collection, true, (bool)(mat_el == 3));

                if(mat_el == 1){
                    if(setup == 0){
                        m.selfDHggcoupl[0][gHIGGS_GG_2][0] = 1;
                        m.selfDHzzcoupl[0][gHIGGS_VV_1][0] = 2;
                        m.selfDHzzcoupl[0][gHIGGS_VV_4][0] = 0;
                    } else if(setup == 1){
                        m.selfDHggcoupl[0][gHIGGS_GG_2][0] = 1;
                        m.selfDHzzcoupl[0][gHIGGS_VV_4][0] = 1;
                        m.selfDHzzcoupl[0][gHIGGS_VV_1][0] = 0;
                    }
                } else if(mat_el == 3){
                    if(setup == 1){
                        // From Lexicon
                        m.mdl_chwtil = -6.34078;
                        m.mdl_chbtil = -1.90674;
                        m.mdl_chg = -8.24752;
                        m.mdl_chwbtil = -6.9542;
                        m.mdl_chbox = -16.495;
                    }
                }
                float ans = 0;
                m.computeP(ans);
                if(mat_el == 1){
                    jhugenprobs[setup].push_back(ans);
                } else if(mat_el == 3){
                    madprobs[setup].push_back(ans);
                }
                m.resetInputEvent();
            }
        }
        m4l.push_back(*M4L);
        counter++;
        delete mother_collection;
        delete daughter_collection;
        delete associated_collection;
    }

    dataFile->Close();
    ofstream outputFile;
    cout << "WRITTEN SUCCESSFULLY" << endl;
    cout << m4l.size() << " EVENTS" << endl;
    outputFile.open("probs_output.csv");
    outputFile << "jhugenprob_g1, jhugenprob_g4, madprob_g1, madprob_g4, m4l\n";
    for(int j = 0; j < m4l.size(); j++){
        outputFile << jhugenprobs[0][j] << ", " << jhugenprobs[1][j] << ", " << madprobs[0][j] << ", " << madprobs[1][j] << ", " << m4l[j] << "\n";
    }
    outputFile.close();

    return 0;
}