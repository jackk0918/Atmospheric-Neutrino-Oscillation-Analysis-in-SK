import ROOT
import os
from time import time
import numpy as np

dir_ge = "/node06_storage/sk_user/shjung/genie_vector/test/GHEPNTP.ATMNU.1000/"
dir_fq = "/node06_storage/sk_user/shjung/fitqun/sk/atmnu/"

h = 1810        # detector half-height [cm]
r = 1690        # detector radius [cm]

t = time()

# Create new output ROOT file
fout = ROOT.TFile("dir_true_vs_recon.root", "RECREATE")

tree_ge = ROOT.TChain("ghepparT")
tree_fq = ROOT.TChain("fiTQun")


# Combine all 1000 trees
for i in range(1, 501):                                         #1, 501
    file_name_ge = "gheppart.%05d.root" % i
    file_path_ge = os.path.join(dir_ge, file_name_ge)
    file_name_fq = "%05d/%05d.merged.root" % (i, i)
    file_path_fq = os.path.join(dir_fq, file_name_fq)

    tree_ge.Add(file_path_ge)
    tree_fq.Add(file_path_fq)

nbins = 100
xmin = 0
xmax = 0

# Initialise histogram objects
h_fqx_e = ROOT.TH1D("h_fqx_e", "e fq, x;x-dir;Events", nbins, xmin, xmax)
h_fqy_e = ROOT.TH1D("h_fqy_e", "e fq, y;y-dir;Events", nbins, xmin, xmax)
h_fqz_e = ROOT.TH1D("h_fqz_e", "e fq, z;z-dir;Events", nbins, xmin, xmax)

h_fqx_mu = ROOT.TH1D("h_fqx_mu", "mu fq, x;x-dir;Events", nbins, xmin, xmax)
h_fqy_mu = ROOT.TH1D("h_fqy_mu", "mu fq, y;y-dir;Events", nbins, xmin, xmax)
h_fqz_mu = ROOT.TH1D("h_fqz_mu", "mu fq, z;z-dir;Events", nbins, xmin, xmax)

h_gex_e = ROOT.TH1D("h_gex_e", "e ge, x;x-dir;Events", nbins, xmin, xmax)
h_gey_e = ROOT.TH1D("h_gey_e", "e ge, y;y-dir;Events", nbins, xmin, xmax)
h_gez_e = ROOT.TH1D("h_gez_e", "e ge, z;z-dir;Events", nbins, xmin, xmax)

h_gex_mu = ROOT.TH1D("h_gex_mu", "mu ge, x;x-dir;Events", nbins, xmin, xmax)
h_gey_mu = ROOT.TH1D("h_gey_mu", "mu ge, y;y-dir;Events", nbins, xmin, xmax)
h_gez_mu = ROOT.TH1D("h_gez_mu", "mu ge, z;z-dir;Events", nbins, xmin, xmax)


n = tree_ge.GetEntries()
m = tree_fq.GetEntries()
if m != n:
    print("The numbers of entries of WC and FQ do not match.")

for i in range(n):

    tree_ge.GetEntry(i)
    tree_fq.GetEntry(i)

    # Getting the ge branches
    momentum = getattr(tree_ge, "momentum")
    momentum = np.array(momentum).reshape(4, 100)
    mom = getattr(tree_ge, "mom")

    # Getting the fq branches
    fq1rdir = getattr(tree_fq, "fq1rdir")
    fqnse = getattr(tree_fq, "fqnse")
    fq1rdir = np.array(fq1rdir).reshape(fqnse, 7, 3)
    fqmrnring = getattr(tree_fq, "fqmrnring")
    fq1rnll = getattr(tree_fq, "fq1rnll")
    fq1rnll = np.array(fq1rnll).reshape(fqnse, 7)

    # Record the vertices
    fqXdirE = fq1rdir[0][1][0]
    fqYdirE = fq1rdir[0][1][1]
    fqZdirE = fq1rdir[0][1][2]

    fqXdirMu = fq1rdir[0][2][0]
    fqYdirMu = fq1rdir[0][2][1]
    fqZdirMu = fq1rdir[0][2][2]

    geXdir = momentum[1][0]/mom[0]
    geYdir = momentum[2][0]/mom[0]
    geZdir = momentum[3][0]/mom[0]

    # Fill the histograms
    if (fqmrnring[0] == 1):
        # electron-like
        if (fq1rnll[0][1] < fq1rnll[0][2] and fq1rnll[0][1] < fq1rnll[0][3]):
            h_fqx_e.Fill(fqXdirE)
            h_fqy_e.Fill(fqYdirE)
            h_fqz_e.Fill(fqZdirE)

            h_gex_e.Fill(geXdir)
            h_gey_e.Fill(geYdir)
            h_gez_e.Fill(geZdir)

        # muon-like
        elif (fq1rnll[0][2] < fq1rnll[0][1]):           #  and fq1rnll[0][2] < fq1rnll[0][3]
            h_fqx_mu.Fill(fqXdirMu)
            h_fqy_mu.Fill(fqYdirMu)
            h_fqz_mu.Fill(fqZdirMu)

            h_gex_mu.Fill(geXdir)
            h_gey_mu.Fill(geYdir)
            h_gez_mu.Fill(geZdir)


fout.cd()


def setup_canvas(canvas_name, h_fq, h_ge, title):
    c = ROOT.TCanvas(canvas_name, "Canvas", 800, 600)

    h_fq.SetLineColor(ROOT.kRed)
    h_fq.Draw()

    h_ge.SetLineColor(ROOT.kBlue)
    h_ge.Draw("SAME")

    # Get the statistics box for fiTQun
    stats_fq = h_fq.FindObject("stats")
    if not stats_fq or not isinstance(stats_fq, ROOT.TPaveStats):
        # If stats box doesn't exist or is not TPaveStats, create a new one
        stats_fq = ROOT.TPaveStats(0.78, 0.86, 0.93, 0.98, "brNDC")
        stats_fq.SetName("stats_fq")
        stats_fq.SetTextColor(ROOT.kRed)
        stats_fq.SetTextSize(0.03)
        h_fq.GetListOfFunctions().Add(stats_fq)
        h_fq.SetStats(1)  # Ensure stats are on

    # Adjust the position of fiTQun stats box
    stats_fq.SetX1NDC(0.78)
    stats_fq.SetX2NDC(0.98)
    stats_fq.SetY1NDC(0.77)
    stats_fq.SetY2NDC(0.97)

    # Create a new statistics box for genie
    stats_ge = ROOT.TPaveStats(0.78, 0.74, 0.93, 0.86, "brNDC")
    stats_ge.SetName("stats_ge")
    stats_ge.SetBorderSize(1)
    stats_ge.SetFillColor(0)
    stats_ge.SetTextAlign(12)
    stats_ge.SetTextColor(ROOT.kBlue)
    stats_ge.SetTextSize(0.03)
    stats_ge.AddText("h_ge")
    stats_ge.AddText(f"Entries = {h_ge.GetEntries()}")
    stats_ge.AddText(f"Mean = {h_ge.GetMean():.0f}")
    stats_ge.AddText(f"Std Dev = {h_ge.GetStdDev():.0f}")
    h_ge.GetListOfFunctions().Add(stats_ge)
    h_ge.SetStats(1)  # Ensure stats are on

    # Draw both stat boxes
    stats_fq.Draw()
    stats_ge.Draw()

    # Create a legend
    legend = ROOT.TLegend(0.7, 0.5, 0.9, 0.55)
    legend.AddEntry(h_fq, "fiTQun", "l")
    legend.AddEntry(h_ge, "WCSim", "l")
    legend.Draw()

    h_fq.SetTitle(title)

    h_fq.Write()
    h_ge.Write()
    c.Write()


# Create canvases
setup_canvas("cx_e", h_fqx_e, h_gex_e, "Electron X-dir: fiTQun vs WCSim")
setup_canvas("cy_e", h_fqy_e, h_gey_e, "Electron Y-dir: fiTQun vs WCSim")
setup_canvas("cz_e", h_fqz_e, h_gez_e, "Electron Z-dir: fiTQun vs WCSim")
setup_canvas("cx_mu", h_fqx_mu, h_gex_mu, "Muon X-dir: fiTQun vs WCSim")
setup_canvas("cy_mu", h_fqy_mu, h_gey_mu, "Muon Y-dir: fiTQun vs WCSim")
setup_canvas("cz_mu", h_fqz_mu, h_gez_mu, "Muon Z-dir: fiTQun vs WCSim")


fout.Close()


print("operation took %f seconds." % (time()-t))
