# Generates e-like & mu-like zenith angle histograms

import ROOT
import os
from time import time
import numpy as np
from osc_prob import *
from tqdm import tqdm  # Import tqdm

def generate1DHistReduced(nbins, xmin, xmax, nfiles):

    dir_ge = "/node06_storage/sk_user/shjung/genie_vector/test/GHEPNTP.ATMNU.1000/"
    dir_wc = "/node06_storage/sk_user/shjung/wcsim/sk/atmnu/ntp/"
    dir_fq = "/node06_storage/sk_user/shjung/fitqun/sk/atmnu/"

    thist = time()

    # Create new output ROOT file
    fout = ROOT.TFile("/home/jack0918/thesis/sensitivity/histograms/histograms_dmsq/36.root", "RECREATE")

    tree_ge = ROOT.TChain("ghepparT")
    tree_wc = ROOT.TChain("wcntp")
    tree_fq = ROOT.TChain("fiTQun")

    # Combine all 1000 trees
    for i in range(1, nfiles + 1):
        file_name_ge = "gheppart.%05d.root" % i
        file_path_ge = os.path.join(dir_ge, file_name_ge)
        file_name_wc = "wcntp.%05d.root" % i
        file_path_wc = os.path.join(dir_wc, file_name_wc)
        file_name_fq = "%05d/%05d.merged.root" % (i, i)
        file_path_fq = os.path.join(dir_fq, file_name_fq)

        tree_ge.Add(file_path_ge)
        tree_wc.Add(file_path_wc)
        tree_fq.Add(file_path_fq)

    # Event Criteria 1: Fiducial Volume (2m away from the detector wall)
    d_height = 1810   # detector half-height [cm]
    d_radius = 1690   # detector radius [cm]

    def fv_criteria(x,y,z):
        return abs(z) < (d_height - 200) and x**2 + y**2 < (d_radius - 200)**2

    # Initialise histograms
    n0e = ROOT.TH1D("e", "Sub-GeV e-like;cos;Events", nbins, xmin, xmax)
    n0mu = ROOT.TH1D("mu", "Sub-GeV mu-like;cos;Events", nbins, xmin, xmax)

    n = tree_ge.GetEntries()
    # -----------------------------------------------
    # Main Loop

    for i in tqdm(range(n), desc="Processing events", unit="event"):

        tree_ge.GetEntry(i)
        tree_wc.GetEntry(i)
        tree_fq.GetEntry(i)

        # gheppart branches for calcNuFact()
        npar = tree_ge.npar
        momentum = np.array(tree_ge.momentum).reshape(npar, 4)
        mom = tree_ge.mom
        trueDir = momentum[0][3]/mom[0]
        trueEnergy = mom[0]

        # wcntp branches for selection criteria (SumQ)
        SumQ = tree_wc.SumQ
        nTrig = tree_wc.nTrig
        sumq = 0
        for j in range(nTrig):
            sumq += SumQ[j]
        Vtx = tree_wc.Vtx

        # fiTQun branches for selection criteria (FV)
        fq1rpos = tree_fq.fq1rpos
        fq1rdir = tree_fq.fq1rdir
        fqnse = tree_fq.fqnse
        fq1rpos = np.array(fq1rpos).reshape(fqnse, 7, 3)
        fq1rdir = np.array(fq1rdir).reshape(fqnse, 7, 3)
        fqmrnring = tree_fq.fqmrnring
        fq1rnll = tree_fq.fq1rnll
        fq1rnll = np.array(fq1rnll).reshape(fqnse, 7)
        fq1rmom = tree_fq.fq1rmom
        fq1rmom = np.array(fq1rmom).reshape(fqnse, 7)

        Emu = np.sqrt(fq1rmom[0][2]**2 + 105.66**2)

        # Reconstructed vertices
        fqXvtxE = fq1rpos[0][1][0]
        fqYvtxE = fq1rpos[0][1][1]
        fqZvtxE = fq1rpos[0][1][2]

        fqXvtxMu = fq1rpos[0][2][0]
        fqYvtxMu = fq1rpos[0][2][1]
        fqZvtxMu = fq1rpos[0][2][2]
        # True vertices
        wcXvtx = Vtx[0]
        wcYvtx = Vtx[1]
        wcZvtx = Vtx[2]

        # Reconstructed directions
        zenithE = fq1rdir[0][1][2]
        zenithMu = fq1rdir[0][2][2]

        # Calculate probabilities one event at a time
        weight_e = calcNuFast(trueDir, trueEnergy, particle=12)
        weight_mu = calcNuFast(trueDir, trueEnergy, particle=14)

        # ---------------------------------------------------------------------------------
        # Fill the histograms
        if (fqmrnring[0] == 1):          # 1-ring cut

            # e-like
            if (fq1rnll[0][1] < fq1rnll[0][2] and fq1rnll[0][1] < fq1rnll[0][3] and
                fv_criteria(fqXvtxE,fqYvtxE,fqZvtxE) and
                sumq > 360 and
                fq1rmom[0][1] > 100): # and fq1rmom[0][1] < 1330): # <---- removed Sub-GeV cut
                n0e.Fill(zenithE, weight_e)

            # mu-like
            elif (fq1rnll[0][2] < fq1rnll[0][1] and
                fv_criteria(fqXvtxMu,fqYvtxMu,fqZvtxMu) and
                sumq > 360 and
                fq1rmom[0][2] > 200): # and Emu < 1330): # <------ removed Sub-GeV cut
                n0mu.Fill(zenithMu, weight_mu)


    factorMu = 10854 / 148213
    factorE = 10294 / 93999

    n0e.Scale(factorE)
    n0mu.Scale(factorMu)

    n0e.SetMinimum(0)
    n0mu.SetMinimum(0)
    n0e.SetMaximum(n0e.GetMaximum() / 0.8)
    n0mu.SetMaximum(n0mu.GetMaximum() / 0.8)

    fout.cd()


    n0e.Write()
    n0mu.Write()
    
    fout.Close()

    print("Flux generation took %f seconds." % (time()-thist))

    return n0e, n0mu

# test (after test, delete fout)
generate1DHistReduced(10, -1.0, 1.0, nfiles=1000)