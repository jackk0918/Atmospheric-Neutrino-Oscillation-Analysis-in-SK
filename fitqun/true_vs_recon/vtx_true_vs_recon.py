import ROOT
import os
from time import time
import numpy as np

dir_wc = "/node06_storage/sk_user/shjung/wcsim/sk/atmnu/ntp/"
dir_fq = "/node06_storage/sk_user/shjung/fitqun/sk/atmnu/"

h = 1810        # detector half-height [cm]
r = 1690        # detector radius [cm]

t = time()

# Create new output ROOT file
fout = ROOT.TFile("vtx_true_vs_recon.root", "RECREATE")

tree_wc = ROOT.TChain("wcntp")
tree_fq = ROOT.TChain("fiTQun")


# Combine all 1000 trees
for i in range(1, 501):                                         #1, 501
    file_name_wc = "wcntp.%05d.root" % i
    file_path_wc = os.path.join(dir_wc, file_name_wc)
    file_name_fq = "%05d/%05d.merged.root" % (i, i)
    file_path_fq = os.path.join(dir_fq, file_name_fq)

    tree_wc.Add(file_path_wc)
    tree_fq.Add(file_path_fq)

nbins = 100
xmin = -2200
xmax = 2200

# Initialise histogram objects
h_fqx_e = ROOT.TH1D("h_fqx_e", "e fq, x;x-vertex;Events", nbins, xmin, xmax)
h_fqy_e = ROOT.TH1D("h_fqy_e", "e fq, y;y-vertex;Events", nbins, xmin, xmax)
h_fqz_e = ROOT.TH1D("h_fqz_e", "e fq, z;z-vertex;Events", nbins, xmin, xmax)

h_fqx_mu = ROOT.TH1D("h_fqx_mu", "mu fq, x;x-vertex;Events", nbins, xmin, xmax)
h_fqy_mu = ROOT.TH1D("h_fqy_mu", "mu fq, y;y-vertex;Events", nbins, xmin, xmax)
h_fqz_mu = ROOT.TH1D("h_fqz_mu", "mu fq, z;z-vertex;Events", nbins, xmin, xmax)

h_wcx_e = ROOT.TH1D("h_wcx_e", "e wc, x;x-vertex;Events", nbins, xmin, xmax)
h_wcy_e = ROOT.TH1D("h_wcy_e", "e wc, y;y-vertex;Events", nbins, xmin, xmax)
h_wcz_e = ROOT.TH1D("h_wcz_e", "e wc, z;z-vertex;Events", nbins, xmin, xmax)

h_wcx_mu = ROOT.TH1D("h_wcx_mu", "mu wc, x;x-vertex;Events", nbins, xmin, xmax)
h_wcy_mu = ROOT.TH1D("h_wcy_mu", "mu wc, y;y-vertex;Events", nbins, xmin, xmax)
h_wcz_mu = ROOT.TH1D("h_wcz_mu", "mu wc, z;z-vertex;Events", nbins, xmin, xmax)


n = tree_wc.GetEntries()
m = tree_fq.GetEntries()
if m != n:
    print("The numbers of entries of WC and FQ do not match.")

for i in range(n):

    tree_wc.GetEntry(i)
    tree_fq.GetEntry(i)

    # Getting the wc branches
    Vtx = getattr(tree_wc, "Vtx")

    # Getting the fq branches
    fq1rpos = getattr(tree_fq, "fq1rpos")
    fqnse = getattr(tree_fq, "fqnse")
    fq1rpos = np.array(fq1rpos).reshape((fqnse, 7, 3))
    fqmrnring = getattr(tree_fq, "fqmrnring")
    fq1rnll = getattr(tree_fq, "fq1rnll")
    fq1rnll = np.array(fq1rnll).reshape(fqnse, 7)

    fqXvtxE = 0
    fqYvtxE = 0
    fqZvtxE = 0
    fqXvtxMu = 0
    fqYvtxMu = 0
    fqZvtxMu = 0

    # Record the vertices
    fqXvtxE = fq1rpos[0][1][0]
    fqYvtxE = fq1rpos[0][1][1]
    fqZvtxE = fq1rpos[0][1][2]

    fqXvtxMu = fq1rpos[0][2][0]
    fqYvtxMu = fq1rpos[0][2][1]
    fqZvtxMu = fq1rpos[0][2][2]

    wcXvtx = Vtx[0]
    wcYvtx = Vtx[1]
    wcZvtx = Vtx[2]

    # Fill the histograms
    if (fqmrnring[0] == 1):
        if (fq1rnll[0][1] < fq1rnll[0][2] and fq1rnll[0][1] < fq1rnll[0][3]):   # and fq1rnll[0][1] < fq1rnll[0][3]
            h_fqx_e.Fill(fqXvtxE)
            h_fqy_e.Fill(fqYvtxE)
            h_fqz_e.Fill(fqZvtxE)

            h_wcx_e.Fill(wcXvtx)
            h_wcy_e.Fill(wcYvtx)
            h_wcz_e.Fill(wcZvtx)

        elif (fq1rnll[0][2] < fq1rnll[0][1]):           #  and fq1rnll[0][2] < fq1rnll[0][3]
            h_fqx_mu.Fill(fqXvtxMu)
            h_fqy_mu.Fill(fqYvtxMu)
            h_fqz_mu.Fill(fqZvtxMu)

            h_wcx_mu.Fill(wcXvtx)
            h_wcy_mu.Fill(wcYvtx)
            h_wcz_mu.Fill(wcZvtx)


fout.cd()


def setup_canvas(canvas_name, h_fq, h_wc, legend1, legend2, title):
    c = ROOT.TCanvas(canvas_name, "Canvas", 800, 600)

    h_fq.SetLineColor(ROOT.kRed)
    h_fq.Draw()

    h_wc.SetLineColor(ROOT.kBlue)
    h_wc.Draw("SAME")

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


    # Create a new statistics box for WCSim
    stats_wc = ROOT.TPaveStats(0.78, 0.74, 0.93, 0.86, "brNDC")
    stats_wc.SetName("stats_wc")
    stats_wc.SetBorderSize(1)
    stats_wc.SetFillColor(0)
    stats_wc.SetTextAlign(12)
    stats_wc.SetTextColor(ROOT.kBlue)
    stats_wc.SetTextSize(0.03)
    stats_wc.AddText(legend2)
    stats_wc.AddText(f"Entries = {h_wc.GetEntries()}")
    stats_wc.AddText(f"Mean = {h_wc.GetMean():.2f}")
    stats_wc.AddText(f"Std Dev = {h_wc.GetStdDev():.2f}")
    h_wc.GetListOfFunctions().Add(stats_wc)
    h_wc.SetStats(1)  # Ensure stats are on

    # Draw both stat boxes
    stats_fq.Draw()
    stats_wc.Draw()

    # Create a legend
    legend = ROOT.TLegend(0.7, 0.5, 0.9, 0.55)
    legend.AddEntry(h_fq, legend1, "l")
    legend.AddEntry(h_wc, legend2, "l")
    legend.Draw()

    h_fq.SetTitle(title)

    h_fq.Write()
    h_wc.Write()
    c.Write()


# Create canvases
setup_canvas("cx_e", h_fqx_e, h_wcx_e, "fq", "wc", "Electron X-vertex: fiTQun vs WCSim")
setup_canvas("cy_e", h_fqy_e, h_wcy_e, "fq", "wc", "Electron Y-vertex: fiTQun vs WCSim")
setup_canvas("cz_e", h_fqz_e, h_wcz_e, "fq", "wc", "Electron Z-vertex: fiTQun vs WCSim")
setup_canvas("cx_mu", h_fqx_mu, h_wcx_mu, "fq", "wc", "Muon X-vertex: fiTQun vs WCSim")
setup_canvas("cy_mu", h_fqy_mu, h_wcy_mu, "fq", "wc", "Muon Y-vertex: fiTQun vs WCSim")
setup_canvas("cz_mu", h_fqz_mu, h_wcz_mu, "fq", "wc", "Muon Z-vertex: fiTQun vs WCSim")


fout.Close()


print("operation took %f seconds." % (time()-t))
