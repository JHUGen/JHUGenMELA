import Mela
import numpy as np
import cProfile
import pstats

associated_test = Mela.SimpleParticleCollection_t(
    [2, 1],
    [491.08166, 274.43945],
    [-0.714883,  1.2373503],
    [2.9122967, -2.589596],
    [0.0006989, 0.0005849],
    True #ptEtaPhiM vector
)

daughters_test = Mela.SimpleParticleCollection_t(
    [13, -13, 11, -11],
    [63.553325, 921.65033, 674.59143, 125.21847],
    [-1.177799, -1.180188, -0.235953, -0.099288],
    [0.6205077, 0.9990979, -1.406654, -1.148683],
    # [0, 0, 0, 0],
    [0.1056602, 0.1057852, 1.7767999, 1.7768000],
    True
)

mothers_test = Mela.SimpleParticleCollection_t(
    [1,2],
    [0,0],
    [0,0],
    [1067.5159, -2641.061],
    [1067.5159, 2641.0612],
    False
)

m = Mela.Mela()

m.setInputEvent(
    daughters_test,
    associated_test,
    mothers_test
)

m.setProcess(
    Mela.Process.bkgZZ,
    Mela.MatrixElement.MCFM,
    Mela.Production.JJEW,
)
MCFM_native = m.computeProdDecP(False)

print("SM value:", MCFM_native)

m.alpha_SMEW=[1,0]
m.dAAWpWm=[1,0]
m.dZAWpWm=[1,0]
m.dZZWpWm=[1,0]
m.dM_A=[1,0]
m.dM_Z=[1,0]
m.dP_A=[1,0] 
m.dV_A=[1,0]
m.dV_Z=[1,0]
m.dP_Z=[1,0]
m.dFour_Z=[1,0]
m.dFour_A=[0,0]
m.cW=[0,0]
MCFM_cW_1 = m.computeProdDecP(False)

m.alpha_SMEW=[1,0]
m.dAAWpWm=[1,0]
m.dZAWpWm=[1,0]
m.dZZWpWm=[1,0]
m.dM_A=[1,0]
m.dM_Z=[1,0]
m.dP_A=[1,0] 
m.dV_A=[1,0]
m.dV_Z=[1,0]
m.dP_Z=[1,0]
m.dFour_Z=[1,0]
m.dFour_A=[0,0]
m.cW=[0,0]
MCFM_cW_m1 = m.computeProdDecP(False)

W_cW_pure = (MCFM_cW_1 + MCFM_cW_m1 - 2*MCFM_native)/2/MCFM_native
print(f"pure MCFM cW=: {W_cW_pure:.3e}")

m.setProcess(
    Mela.Process.bkgZZ,
    Mela.MatrixElement.MADGRAPH,
    Mela.Production.JJEW,
)
m.mdl_cW = 1
madgraph_cW_1 = m.computeP(False)

m.mdl_cHWBtil = -1
madgraph_cW_m1 = m.computeP(False)

madgraph_NATIVE = m.computeP(False)

W_chwbtil_pure = (madgraph_cW_1 + madgraph_cW_m1 - 2*madgraph_NATIVE)/2/madgraph_NATIVE
print(f"pure MADGRAPH cW=: {W_chwbtil_pure:.3e}")