# Test code to draw any distribution I need to check
# 1. plot ccqe single e, mu momentum & sumq & momentum vs sumq distributions after FV cut
# 2. plot energy (momentum) distributions after each cut

import ROOT
import os
from time import time
import numpy as np

dir_ge = "/node06_storage/sk_user/shjung/genie_vector/test/GHEPNTP.ATMNU.1000/"
dir_wc = "/node06_storage/sk_user/shjung/wcsim/sk/atmnu/ntp/"
dir_fq = "/node06_storage/sk_user/shjung/fitqun/sk/atmnu/"

height= 1810    # detector half-height [cm]
radius = 1690   # detector radius [cm]

t = time()

# Create new output ROOT file
fout = ROOT.TFile("distributions.root", "RECREATE")

tree_ge = ROOT.TChain("ghepparT")
tree_wc = ROOT.TChain("wcntp")
tree_fq = ROOT.TChain("fiTQun")


# Combine all 1000 trees
for i in range(1, 51):
    file_name_ge = "gheppart.%05d.root" % i
    file_path_ge = os.path.join(dir_ge, file_name_ge)
    file_name_wc = "wcntp.%05d.root" % i
    file_path_wc = os.path.join(dir_wc, file_name_wc)
    file_name_fq = "%05d/%05d.merged.root" % (i, i)
    file_path_fq = os.path.join(dir_fq, file_name_fq)

    tree_ge.Add(file_path_ge)
    tree_wc.Add(file_path_wc)
    tree_fq.Add(file_path_fq)

print("Trees combined in %f seconds." % (time()-t))

nbins = 1000
xmin = 0
xmax = 0


h_mom_e = ROOT.TH1D("h_mom_e", "mom_e;momentum (GeV/c);Events", nbins, xmin, xmax)
h_mom_mu  = ROOT.TH1D("h_mom_mu", "mom_mu;momentum (GeV/c);Events", nbins, xmin, xmax)
h_sumq_e = ROOT.TH1D("h_sumq_e", "sumq_e;sumq (p.e.);Events", nbins, xmin, xmax)
h_sumq_mu = ROOT.TH1D("h_sumq_mu", "sumq_mu;sumq (p.e.);Events", nbins, xmin, xmax)
h_ms_e = ROOT.TH2D("h_ms_e", "e momentum vs sumq;sumq (p.e.);momentum (GeV/c)", 100, 0, 0, 100, 0, 0)
h_ms_mu = ROOT.TH2D("h_ms_mu", "mu momentum vs sumq;sumq (p.e.);momentum (GeV/c)", 100, 0, 0, 100, 0, 0)

h_mom_uncut_e  = ROOT.TH1D("h_mom_uncut_e", "electron momentum;momentum (GeV/c);Events", nbins, xmin, xmax)
h_mom_uncut_mu  = ROOT.TH1D("h_mom_uncut_mu", "muon momentum;momentum (GeV/c);Events", nbins, xmin, xmax)
h_mom_cut_e  = ROOT.TH1D("h_mom_cut_e", "mom_cut_e;momentum (GeV/c);Events", nbins, xmin, xmax)
h_mom_cut_mu  = ROOT.TH1D("h_mom_cut_mu", "mom_cut_mu;momentum (GeV/c);Events", nbins, xmin, xmax)

h_Mom_e = ROOT.TH1D("h_Mom_e", "Mom_e;Momentum (MeV/c);Events", nbins, xmin, xmax)
h_Mom_mu = ROOT.TH1D("h_Mom_mu", "Mom_mu;Momentum (MeV/c);Events", nbins, xmin, xmax)
h_Energy_e = ROOT.TH1D("h_Energy_e", "Energy_e;Energy (MeV);Events", nbins, xmin, xmax)
h_Energy_mu = ROOT.TH1D("h_Energy_mu", "Energy_mu;Energy (MeV);Events", nbins, xmin, xmax)

n = tree_ge.GetEntries()


# Event Criteria 1: Fiducial Volume (2m away from the detector wall)
def fv_criteria(x,y,z):
   return abs(z) < (height - 200) and x**2 + y**2 < (radius - 200)**2


# Is electron or muon
def isLepton(index, lepton):
    if (lepton == 0):
        return (pdg[index] == 11 or pdg[index] == -11)
    elif (lepton == 1):
        return (pdg[index] == 13 or pdg[index] == -13)



for i in range(n):

    tree_ge.GetEntry(i)
    tree_wc.GetEntry(i)
    tree_fq.GetEntry(i)

    # Getting the ge branches
    npar = getattr(tree_ge, "npar")
    momentum = getattr(tree_ge, "momentum")
    momentum = np.array(momentum)
    momentum = momentum.reshape(npar, 4)
    mom = getattr(tree_ge, "mom")
    inttype = getattr(tree_ge, "inttype")
    scattype = getattr(tree_ge, "scattype")
    pdg = getattr(tree_ge, "pdg")

    # Getting the wc branches
    Vtx = getattr(tree_wc, "Vtx")
    Npar = getattr(tree_wc, "Npar")
    nTrig = getattr(tree_wc, "nTrig")
    SumQ = getattr(tree_wc, "SumQ")
    sumq = 0
    for j in range(nTrig):
        sumq += SumQ[j]

    # Getting the fq branches
    fq1rpos = getattr(tree_fq, "fq1rpos")
    fq1rdir = getattr(tree_fq, "fq1rdir")
    fqnse = getattr(tree_fq, "fqnse")
    fq1rpos = np.array(fq1rpos).reshape((fqnse, 7, 3))
    fq1rdir = np.array(fq1rdir).reshape((fqnse, 7, 3))
    fqmrnring = getattr(tree_fq, "fqmrnring")
    fq1rnll = getattr(tree_fq, "fq1rnll")
    fq1rnll = np.array(fq1rnll).reshape(fqnse, 7)
    fq1rmom = getattr(tree_fq, "fq1rmom")
    fq1rmom = np.array(fq1rmom).reshape((fqnse, 7))


    # Reconstructed vertices (for the FV criteria)
    fqXvtxE = fq1rpos[0][1][0]
    fqYvtxE = fq1rpos[0][1][1]
    fqZvtxE = fq1rpos[0][1][2]

    fqXvtxMu = fq1rpos[0][2][0]
    fqYvtxMu = fq1rpos[0][2][1]
    fqZvtxMu = fq1rpos[0][2][2]

    # true vertices
    wcXvtx = Vtx[0]
    wcYvtx = Vtx[1]
    wcZvtx = Vtx[2]

    # identify primary lepton index
    l_index = -1
    for j in range(npar):
        if isLepton(j, 0) or isLepton(j, 1):
            l_index = j
            break
    if (l_index == -1):
        continue


    # identify primary lepton index (wcsim)
    L_index = -1
    for j in range(Npar):
        if isLepton(j, 0) or isLepton(j, 1):
            L_index = j
            break
    if (L_index == -1):
        continue


    # ---------------------------------------------------------------------------------
    # Fill the histograms
    if (fqmrnring[0] == 1 and inttype == 2 and scattype == 1):          # 1-ring cut + CCQE cut

        # e-like
        if (fq1rnll[0][1] < fq1rnll[0][2] and fq1rnll[0][1] < fq1rnll[0][3] and isLepton(l_index, 0)):
            h_mom_uncut_e.Fill(mom[l_index])

            # Selection Criteria (FV & SumQ Cut)
            if (fv_criteria(fqXvtxE,fqYvtxE,fqZvtxE) and fv_criteria(wcXvtx,wcYvtx,wcZvtx)):
                h_sumq_e.Fill(sumq)
                h_mom_e.Fill(mom[l_index])
                if (sumq < 4000):
                    h_ms_e.Fill(sumq, mom[l_index])

                if (sumq > 360):        # Criteria 2: SumQ > 360 p.e.
                    h_mom_cut_e.Fill(mom[l_index])
                    h_Mom_e.Fill(fq1rmom[0][1])


        # mu-like
        elif (fq1rnll[0][2] < fq1rnll[0][1] and isLepton(l_index, 1)):
            h_mom_uncut_mu.Fill(mom[l_index])

            # Selection Criteria (FV & SumQ Cut)
            if (fv_criteria(fqXvtxMu,fqYvtxMu,fqZvtxMu) and fv_criteria(wcXvtx,wcYvtx,wcZvtx)):
                h_sumq_mu.Fill(sumq)
                h_mom_mu.Fill(mom[l_index])
                if (sumq < 4000):
                    h_ms_mu.Fill(sumq, mom[l_index])

                if (sumq > 360):
                    h_mom_cut_mu.Fill(mom[l_index])
                    h_Mom_e.Fill(fq1rmom[0][1])


fout.cd()

h_mom_e.Write()
h_mom_mu.Write()
h_sumq_e.Write()
h_sumq_mu.Write()
h_ms_e.Write()
h_ms_mu.Write()

h_Mom_e.Write()
h_Mom_mu.Write()
h_Energy_e.Write()
h_Energy_mu.Write()

'''
# ---------------------------------------------------------------------------------
# draw momentum distributions

# Function to calculate and return mean - 1 sigma
def get_mean_plus_sigma(hist):
    return hist.GetMean() + hist.GetStdDev()



canvasE = ROOT.TCanvas("momentum", "momentum", 1200, 600)

canvasE.Divide(2, 1)  # 2 columns, 1 row

# Set different colors
h_mom_uncut_e.SetLineColor(ROOT.kBlue)
h_mom_e.SetLineColor(ROOT.kRed)
h_mom_cut_e.SetLineColor(ROOT.kGreen)
h_mom_uncut_mu.SetLineColor(ROOT.kBlue)
h_mom_mu.SetLineColor(ROOT.kRed)
h_mom_cut_mu.SetLineColor(ROOT.kGreen)


h_mom_uncut_e.SetStats(0)
h_mom_e.SetStats(0)
h_mom_cut_e.SetStats(0)
h_mom_uncut_mu.SetStats(0)
h_mom_mu.SetStats(0)
h_mom_cut_mu.SetStats(0)

# e energy pad
canvasE.cd(1)
h_mom_uncut_e.Draw()
h_mom_e.Draw("SAME")
h_mom_cut_e.Draw("SAME")
value_ee = get_mean_plus_sigma(h_mom_uncut_e)
value_ee_cut1 = get_mean_plus_sigma(h_mom_e)
value_ee_cut2 = get_mean_plus_sigma(h_mom_cut_e)
label_ee = ROOT.TLatex(0.45, 0.85, f"#mu + 1#sigma = {value_ee:.2f}")
label_ee_cut1 = ROOT.TLatex(0.45, 0.80, f"#mu + 1#sigma (cut1) = {value_ee_cut1:.2f}")
label_ee_cut2 = ROOT.TLatex(0.45, 0.75, f"#mu + 1#sigma (cut2) = {value_ee_cut2:.2f}")
label_ee.SetTextSize(0.025)
label_ee_cut1.SetTextSize(0.025)
label_ee_cut2.SetTextSize(0.025)
label_ee.SetNDC()
label_ee_cut1.SetNDC()
label_ee_cut2.SetNDC()
label_ee.Draw()
label_ee_cut1.Draw()
label_ee_cut2.Draw()

legende1 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legende1.AddEntry(h_mom_uncut_e, "uncut", "l")
legende1.AddEntry(h_mom_e, "FV cut", "l")
legende1.AddEntry(h_mom_cut_e, "FV+SumQ cut", "l")
legende1.Draw()


# mu energy pad
canvasE.cd(2)
h_mom_uncut_mu.Draw()
h_mom_mu.Draw("SAME")
h_mom_cut_mu.Draw("SAME")
value_emu = get_mean_plus_sigma(h_mom_uncut_mu)
value_emu_cut1 = get_mean_plus_sigma(h_mom_mu)
value_emu_cut2 = get_mean_plus_sigma(h_mom_cut_mu)
label_emu = ROOT.TLatex(0.45, 0.85, f"#mu + 1#sigma = {value_emu:.2f}")
label_emu_cut1 = ROOT.TLatex(0.45, 0.80, f"#mu + 1#sigma (cut1) = {value_emu_cut1:.2f}")
label_emu_cut2 = ROOT.TLatex(0.45, 0.75, f"#mu + 1#sigma (cut2) = {value_emu_cut2:.2f}")
label_emu.SetTextSize(0.025)
label_emu_cut1.SetTextSize(0.025)
label_emu_cut2.SetTextSize(0.025)
label_emu.SetNDC()
label_emu_cut1.SetNDC()
label_emu_cut2.SetNDC()
label_emu.Draw()
label_emu_cut1.Draw()
label_emu_cut2.Draw()

legende2 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legende2.AddEntry(h_mom_uncut_mu, "uncut", "l")
legende2.AddEntry(h_mom_mu, "FV cut", "l")
legende2.AddEntry(h_mom_cut_mu, "FV+SumQ cut", "l")
legende2.Draw()


canvasE.Update()
canvasE.Write()

'''


print("operation completed in %f seconds." % (time()-t))