import ROOT
import os
from time import time
import numpy as np


t = time()

# Create new output ROOT file
fout = ROOT.TFile("e30MeVq.root", "RECREATE")

tree = ROOT.TChain("wcntp")
tree.Add("/node06_storage/reno_usr/jhkang/e30ntp.root")


n = tree.GetEntries()

hsumqe30mev = ROOT.TH1D("sumq_e_30MeV", "30MeV electron SumQ;SumQ [photoelectrons];Entries", 100, 0, 0)
hphoton = ROOT.TH1D("nphoton_e_30MeV", "30MeV electron nPhoton;nPhoton;Entries", 100, 0 ,0)

for i in range(n):

    tree.GetEntry(i)

    sumq = 0

    # Getting the branches
    nTrig = tree.nTrig
    SumQ = tree.SumQ
    nPhoton = tree.nPhoton

    for j in range(nTrig):
        sumq += SumQ[j]

    hsumqe30mev.Fill(sumq)
    hphoton.Fill(nPhoton)


fout.cd()

hsumqe30mev.Write()
hphoton.Write()



fout.Close()


print("operation took %f seconds." % (time()-t))