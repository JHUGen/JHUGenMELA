import Mela
import uproot
import numpy
import tqdm
# particles = [
#         [0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03 ],
#         [0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03 ],
#         [0.8855133E+02, -0.2210069E+02,  0.4008035E+02, -0.7580543E+02 ],
#         [0.3283294E+03, -0.1038496E+03, -0.3019338E+03,  0.7649492E+02 ],
#         [0.1523581E+03, -0.1058810E+03, -0.9770964E+02,  0.4954839E+02 ],
#         [0.4307611E+03,  0.2318313E+03,  0.3595630E+03, -0.5023788E+02 ]
#     ]
# particles = [
#         [5.37086667033E+02, 0.00000000000E+00,  0.00000000000E+00,  5.37086667033E+02  ],
#         [7.27293062837E+00, 0.00000000000E+00,  0.00000000000E+00, -7.27293062837E+00  ],
#         [2.11966451221E+01,  1.78742924427E+01, -1.12102312034E+01, -2.03424438331E+00 ],
#         [7.02133017289E+01,  3.86248305415E+01, -3.12040320381E+01, -4.96421050079E+01 ],
#         [9.29989707133E+00,  5.37669110106E+00, -4.15860490296E-01, -7.57669704638E+00 ],
#         [1.23934598003E+02, -4.84989558626E+01, -4.60240644895E+01, -1.04352390245E+02 ]
#     ]

# ids = [21, 21, -11, 11, -11, 11]

m = Mela.Mela(13, 125, Mela.VerbosityLevel.SILENT)
m.setProcess(Mela.Process.SelfDefine_spin0, Mela.MatrixElement.MADGRAPH, Mela.Production.ZZGG)

dataFile = uproot.open("SM_HIGGS_JHUGEN.root")
tree = dataFile["tree"]

NEVENT = tree.num_entries

data = tree.arrays(
    [
        "LHEDaughterId", "LHEDaughterPt", "LHEDaughterEta", "LHEDaughterPhi", "LHEDaughterMass", 
        "LHEAssociatedParticleId", "LHEAssociatedParticlePt", "LHEAssociatedParticleEta", "LHEAssociatedParticlePhi", "LHEAssociatedParticleMass",
        "LHEMotherId", "LHEMotherPz", "LHEMotherE", "M4L", "Phid", "costheta1d", "costheta2d", "Phi1d"
    ],
    library='np', entry_stop=NEVENT + 1
)

f = open("probs_output.csv", "w+")
probs = [
    [[], []], #jhugen
    [[], []]  #madgraph
]

jhugen_couplings = [
    {
        "ghg2":[1,0],
        "ghz1":[1,0]
    },
    {
        "ghg2":[1,0],
        "ghz4":[1,0]
    }
]

madgraph_couplings = [
    {
        
    },
    {
        "mdl_chbtil":1,
        "mdl_chwtil":1,
        "mdl_chwbtil":1
    }
]
for n, (jhucoupl, madcoupl) in enumerate(zip(jhugen_couplings, madgraph_couplings)):
    for mat_el in [Mela.MatrixElement.JHUGen, Mela.MatrixElement.MADGRAPH]:
        if mat_el.value == 1:
            coupl = jhucoupl
        elif mat_el.value == 3:
            coupl = madcoupl

        for i in tqdm.tqdm(range(NEVENT), total=NEVENT):   
            mothers = Mela.SimpleParticleCollection_t(
                data["LHEMotherId"][i], 
                [0]*len(data["LHEMotherId"][i]), 
                [0]*len(data["LHEMotherId"][i]), 
                data["LHEMotherPz"][i], 
                data["LHEMotherE"][i], 
                False
            )
            
            daughters = Mela.SimpleParticleCollection_t(
                data["LHEDaughterId"][i], 
                data["LHEDaughterPt"][i], 
                data["LHEDaughterEta"][i], 
                data["LHEDaughterPhi"][i], 
                data["LHEDaughterMass"][i], 
                True
            )
            
            # associated = Mela.SimpleParticleCollection_t(
            #     data["LHEAssociatedParticleId"][i], 
            #     data["LHEAssociatedParticlePt"][i], 
            #     data["LHEAssociatedParticleEta"][i], 
            #     data["LHEAssociatedParticlePhi"][i], 
            #     data["LHEAssociatedParticleMass"][i], 
            #     True
            # )
            
            associated = Mela.SimpleParticleCollection_t()
            
            m.setProcess(Mela.Process.SelfDefine_spin0, mat_el, Mela.Production.ZZGG)
            m.setInputEvent(daughters, associated, mothers, True, (mat_el.value == 3))
            
            for coupling, value in coupl.items():
                setattr(m, coupling, value)
            
            if mat_el.value == 1:
                probs[0][n].append(m.computeP(False))
            elif mat_el.value == 3:
                probs[1][n].append(m.computeP(False))
            else:
                raise ValueError("ERROR!")
            m.resetInputEvent()

f.write("jhugenprob_g1, jhugenprob_g4, madprob_g1, madprob_g4, m4l\n")

# for n, (j,m) in enumerate(zip(probs[0], probs[1])):
#     f.write(f"{j[0]}, {j[1]}, {m[0]}, {j[1]}, {data['M4L'][n]} \n")

for i in range(NEVENT):
    f.write(f"{probs[0][0][i]}, {probs[0][1][i]}, {probs[1][0][i]}, {probs[1][1][i]}, {data['M4L'][i]} \n")

f.close()