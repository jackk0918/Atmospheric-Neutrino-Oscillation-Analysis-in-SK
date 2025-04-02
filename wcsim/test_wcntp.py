import ROOT
import os
from time import time
import numpy as np

dir = "/node06_storage/sk_user/shjung/wcsim/sk/atmnu/ntp/"
dir2 = "/node06_storage/sk_user/shjung/genie_vector/test/GHEPNTP.ATMNU.1000/"

t = time()

# Create new output ROOT file
fout = ROOT.TFile("wcntp_hists.root", "RECREATE")

tree = ROOT.TChain("wcntp")
tree2 = ROOT.TChain("ghepparT")

# Combine all 1000 trees
for i in range(1, 1001):
    file_name = "wcntp.%05d.root" % i
    file_path = os.path.join(dir, file_name)
    file_name2 = "gheppart.%05d.root" % i
    file_path2 = os.path.join(dir2, file_name2)

    tree.Add(file_path)
    tree2.Add(file_path2)

n = tree.GetEntries()

# Initialise histogram objects
'''
hsumq1 = ROOT.TH1D("sumq_ccqe", "SumQ (CCQE);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig1 = ROOT.TH1D("ntrig_ccqe", "nTrig (CCQE);nTrig;Entries", 6, 0, 6)
henergy1 = ROOT.TH1D("energy_ccqe", "Energy (CCQE); Energy [MeV]; Entries", 100, 100, 100000)

hsumq2 = ROOT.TH1D("sumq_dis", "SumQ (DIS);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig2 = ROOT.TH1D("ntrig_dis", "nTrig (DIS);nTrig;Entries", 6, 0, 6)
henergy2 = ROOT.TH1D("energy_dis", "Energy (DIS); Energy [MeV]; Entries", 100, 100, 100000)

hsumq3 = ROOT.TH1D("sumq_res", "SumQ (RES);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig3 = ROOT.TH1D("ntrig_res", "nTrig (RES);nTrig;Entries", 6, 0, 6)
henergy3 = ROOT.TH1D("energy_res", "Energy (RES); Energy [MeV]; Entries", 100, 100, 100000)

hsumq4 = ROOT.TH1D("sumq_coh", "SumQ (COH);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig4 = ROOT.TH1D("ntrig_coh", "nTrig (COH);nTrig;Entries", 6, 0, 6)
henergy4 = ROOT.TH1D("energy_coh", "Energy (COH); Energy [MeV]; Entries", 100, 100, 100000)

hsumq5 = ROOT.TH1D("sumq_dfr", "SumQ (DFR);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig5 = ROOT.TH1D("ntrig_dfr", "nTrig (DFR);nTrig;Entries", 6, 0, 6)
henergy5 = ROOT.TH1D("energy_dfr", "Energy (DFR); Energy [MeV]; Entries", 100, 100, 100000)

hsumq6 = ROOT.TH1D("sumq_mec", "SumQ (MEC);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hntrig6 = ROOT.TH1D("ntrig_mec", "nTrig (MEC);nTrig;Entries", 6, 0, 6)
henergy6 = ROOT.TH1D("energy_mec", "Energy (MEC); Energy [MeV]; Entries", 100, 100, 100000)

hsumq1elike = ROOT.TH1D("sumq_ccqe_elike", "e-like SumQ (CCQE);SumQ [photoelectrons];Entries", 100, 0, 1000000)
hsumq1mulike = ROOT.TH1D("sumq_ccqe_mulike", "#mu-like SumQ (CCQE);SumQ [photoelectrons];Entries", 100, 0, 1000000)
'''
hsumqe30mev = ROOT.TH1D("sumq_ccqe_e_30MeV", "SumQ of 30 MeV e-like CCQE events;SumQ [photoelectrons];Entries", 100, 0, 0)
hsumq = ROOT.TH1D("sumq_ccqe_e", "SumQ of all e-like CCQE events;SumQ [photoelectrons];Entries", 100, 0, 0)
he = ROOT.TH1D("mom_ccqe_e", "Energy of all e-like CCQE events;Energy [MeV];Entries", 100, 0, 0)

for i in range(n):

    tree.GetEntry(i)
    tree2.GetEntry(i)

    sumq = 0
    energy = 0

    # Getting the branches
    nTrig = getattr(tree, "nTrig")
    SumQ = getattr(tree, "SumQ")
    Npar = getattr(tree, "Npar")
    Energy = getattr(tree, "Energy")

    npar = getattr(tree2, "npar")
    pdg = getattr(tree2, "pdg")
    inttype = getattr(tree2, "inttype")
    scattype = getattr(tree2, "scattype")
    momentum = getattr(tree2, "momentum")
    momentum = np.array(momentum).reshape(npar, 4)
    mom = getattr(tree2, "mom")

    ntrig = nTrig



    for j in range(nTrig):
        sumq += SumQ[j]


    def isLepton(index, lepton):
        if (lepton == 0):
            return (pdg[index] == 11 or pdg[index] == -11)
        elif (lepton == 1):
            return (pdg[index] == 13 or pdg[index] == -13)

    if (inttype != 2 or scattype != 1):
        continue

    # find the primary lepton
    l_index = next((j for j in range(npar) if isLepton(j, 0)), -1)
    if l_index == -1:
        continue

    # Get the lepton energy
    energyE = mom[l_index]



    # Fill the histograms according to the selection criteria

    if (energyE >=0.029 and energyE < 0.031):
        hsumqe30mev.Fill(sumq)
    hsumq.Fill(sumq)
    he.Fill(energyE)


    '''
    if (inttype == 2 and scattype == 1):
        # scattype in the order of from 0: unknown to GENIE / QES / 1Kaon / DIS / RES / COH
        # / DFR / NuEEL / IMD / AMNuGamma / MEC / CEvNS / IBD / GLR / IMDAnh / PhotonCOH /
        # PhotonRES / DMEL / DMDIS / DME / unknown

        hsumq1.Fill(sumq)
        hntrig1.Fill(ntrig)
        henergy1.Fill(energy)

        if (pdg[4] == 11):
            hsumq1elike.Fill(sumq)
        elif (pdg[4] == 13):
            hsumq1mulike.Fill(sumq)

    elif (scattype == 3):
        hsumq2.Fill(sumq)
        hntrig2.Fill(ntrig)
        henergy2.Fill(energy)

    elif (scattype == 4):
        hsumq3.Fill(sumq)
        hntrig3.Fill(ntrig)
        henergy3.Fill(energy)

    elif (scattype == 5):
        hsumq4.Fill(sumq)
        hntrig4.Fill(ntrig)
        henergy4.Fill(energy)

    elif (scattype == 6):
        hsumq5.Fill(sumq)
        hntrig5.Fill(ntrig)
        henergy5.Fill(energy)

    elif (scattype == 10):
        hsumq6.Fill(sumq)
        hntrig6.Fill(ntrig)
        henergy6.Fill(energy)
    '''


fout.cd()
'''
hsumq1.Write()
hntrig1.Write()
henergy1.Write()
hsumq2.Write()
hntrig2.Write()
henergy2.Write()
hsumq3.Write()
hntrig3.Write()
henergy3.Write()
hsumq4.Write()
hntrig4.Write()
henergy4.Write()
hsumq5.Write()
hntrig5.Write()
henergy5.Write()
hsumq6.Write()
hntrig6.Write()
henergy6.Write()

hsumq1elike.Write()
hsumq1mulike.Write()
'''
hsumqe30mev.Write()
hsumq.Write()
he.Write()

fout.Close()


print("operation took %f seconds." % (time()-t))