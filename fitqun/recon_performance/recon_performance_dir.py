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
fout = ROOT.TFile("recon_performance_dir.root", "RECREATE")

tree_ge = ROOT.TChain("ghepparT")
tree_wc = ROOT.TChain("wcntp")
tree_fq = ROOT.TChain("fiTQun")


# Combine all 1000 trees
for i in range(1, 1001):
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

nbins = 100
xmin = 0
xmax = 0

# Initialise histograms
h_derr_e = ROOT.TH1D("h_dir_e", "electron direction error;deg;Events", nbins, xmin, xmax)
h_derr_mu = ROOT.TH1D("h_dir_mu", "muon direction error;deg;Events", nbins, xmin, xmax)
h_derr_e_cut = ROOT.TH1D("h_dir_e_cut", "electron direction error (cut);deg;Events", nbins, xmin, xmax)
h_derr_mu_cut = ROOT.TH1D("h_dir_mu_cut", "muon direction error(cut);deg;Events", nbins, xmin, xmax)

n = tree_ge.GetEntries()
m = tree_fq.GetEntries()
if m != n:
    print("The numbers of entries of GE and FQ do not match.")


# ---------------------------------------------------------------------------------
# define necessary functions

# Event Criteria 1: Fiducial Volume (2m away from the detector wall)
def fv_criteria(x,y,z):
   return abs(z) < (height - 200) and x**2 + y**2 < (radius - 200)**2


# Cos theta formula
def angleDiff(x1, y1, z1, x2, y2, z2):
    # Calculate dot product
    dot_product = x1*x2 + y1*y2 + z1*z2

    # Calculate magnitudes
    magnitude1 = np.sqrt(x1**2 + y1**2 + z1**2)
    magnitude2 = np.sqrt(x2**2 + y2**2 + z2**2)

    # Calculate cosine of the angle
    cos_theta = dot_product / (magnitude1 * magnitude2)

    # Calculate angle in radians
    angle_rad = np.arccos(cos_theta)

    # Convert to degrees
    angle_deg = np.degrees(angle_rad)

    return angle_deg


# Is electron or muon
def isLepton(index, lepton):
    if (lepton == 0):
        return (pdg[index] == 11 or pdg[index] == -11)
    elif (lepton == 1):
        return (pdg[index] == 13 or pdg[index] == -13)


# ---------------------------------------------------------------------------------
# main loop

tloop = time()

for i in range(n):

    tree_ge.GetEntry(i)
    tree_wc.GetEntry(i)
    tree_fq.GetEntry(i)

    # Getting the ge branches
    npar = getattr(tree_ge, "npar")
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
    fq1rpos = np.array(fq1rpos).reshape(fqnse, 7, 3)
    fq1rdir = np.array(fq1rdir).reshape(fqnse, 7, 3)
    fqmrnring = getattr(tree_fq, "fqmrnring")
    fq1rnll = getattr(tree_fq, "fq1rnll")
    fq1rnll = np.array(fq1rnll).reshape(fqnse, 7)
    fq1rmom = getattr(tree_fq, "fq1rmom")
    fq1rmom = np.array(fq1rmom).reshape(fqnse, 7)

    Emu = np.sqrt(fq1rmom[0][2]**2 + 105.66**2)

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

    # Reconstructed directions
    fqXdirE = fq1rdir[0][1][0]
    fqYdirE = fq1rdir[0][1][1]
    fqZdirE = fq1rdir[0][1][2]

    fqXdirMu = fq1rdir[0][2][0]
    fqYdirMu = fq1rdir[0][2][1]
    fqZdirMu = fq1rdir[0][2][2]


    # identify primary lepton index
    l_index = -1
    for j in range(npar):
        if isLepton(j, 0) or isLepton(j, 1):
            l_index = j
            break
    if (l_index == -1):
        continue

    # primary lepton true direction
    geXdir = momentum[l_index][1]/mom[l_index]
    geYdir = momentum[l_index][2]/mom[l_index]
    geZdir = momentum[l_index][3]/mom[l_index]


    # ---------------------------------------------------------------------------------
    # Fill the histograms
    if (fqmrnring[0] == 1 and inttype == 2 and scattype == 1):          # 1-ring cut + CCQE cut

        # e-like
        if (fq1rnll[0][1] < fq1rnll[0][2] and fq1rnll[0][1] < fq1rnll[0][3] and isLepton(l_index, 0)):
            derr_e = angleDiff(fqXdirE,fqYdirE,fqZdirE,geXdir,geYdir,geZdir)
            h_derr_e.Fill(derr_e)

            # Selection Criteria (FV & SumQ Cut)
            if (fv_criteria(fqXvtxE,fqYvtxE,fqZvtxE) and fv_criteria(wcXvtx,wcYvtx,wcZvtx) and
                sumq > 360 and fq1rmom[0][1] > 100 and fq1rmom[0][1] < 1330):   # Criteria 2: SumQ, 3: Sub-GeV
                    derr_e_cut = angleDiff(fqXdirE,fqYdirE,fqZdirE,geXdir,geYdir,geZdir)
                    h_derr_e_cut.Fill(derr_e_cut)


        # mu-like
        elif (fq1rnll[0][2] < fq1rnll[0][1] and isLepton(l_index, 1)):
            derr_mu = angleDiff(fqXdirMu,fqYdirMu,fqZdirMu,geXdir,geYdir,geZdir)
            h_derr_mu.Fill(derr_mu)

            # Selection Criteria (FV & SumQ Cut)
            if (fv_criteria(fqXvtxMu,fqYvtxMu,fqZvtxMu) and fv_criteria(wcXvtx,wcYvtx,wcZvtx) and
                sumq > 360 and fq1rmom[0][2] > 200 and Emu < 1330):
                    derr_mu_cut = angleDiff(fqXdirMu,fqYdirMu,fqZdirMu,geXdir,geYdir,geZdir)
                    h_derr_mu_cut.Fill(derr_mu_cut)



fout.cd()

print("Main loop completed in %f seconds." % (time()-tloop))



# ---------------------------------------------------------------------------------
# draw errors

tdrawerr = time()



# Function to calculate and return mean - 1 sigma
def get_mean_plus_sigma(hist):
    return hist.GetMean() + hist.GetStdDev()


canvas = ROOT.TCanvas("c1", "Direction Performance", 1200, 600)

canvas.Divide(2, 1)  # 2 columns, 1 row

h_derr_e.SetLineColor(ROOT.kBlue)
h_derr_e_cut.SetLineColor(ROOT.kRed)
h_derr_mu.SetLineColor(ROOT.kBlue)
h_derr_mu_cut.SetLineColor(ROOT.kRed)

h_derr_e.SetStats(0)
h_derr_e_cut.SetStats(0)
h_derr_mu.SetStats(0)
h_derr_mu_cut.SetStats(0)


# e pad
canvas.cd(1)  # Switch to the first pad
h_derr_e.Draw()
h_derr_e_cut.Draw("SAME")
value_e = get_mean_plus_sigma(h_derr_e)
value_e_cut = get_mean_plus_sigma(h_derr_e_cut)
label_e = ROOT.TLatex(0.45, 0.85, f"#mu + 1#sigma = {value_e:.2f}")
label_e_cut = ROOT.TLatex(0.45, 0.80, f"#mu + 1#sigma (cut) = {value_e_cut:.2f}")
label_e.SetTextSize(0.025)
label_e_cut.SetTextSize(0.025)
label_e.SetNDC()
label_e_cut.SetNDC()
label_e.Draw()
label_e_cut.Draw()
label_e2 = ROOT.TLatex(0.45, 0.75, f"Events: {h_derr_e.GetEntries():.0f} --> {h_derr_e_cut.GetEntries():.0f}")
label_e2.SetTextSize(0.025)
label_e2.SetNDC()
label_e2.Draw()

legend1 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend1.AddEntry(h_derr_e, "uncut", "l")
legend1.AddEntry(h_derr_e_cut, "cut", "l")
legend1.Draw()


# mu pad
canvas.cd(2)  # Switch to the second pad
h_derr_mu.Draw()
h_derr_mu_cut.Draw("SAME")
value_mu = get_mean_plus_sigma(h_derr_mu)
value_mu_cut = get_mean_plus_sigma(h_derr_mu_cut)
label_mu = ROOT.TLatex(0.45, 0.85, f"#mu + 1#sigma = {value_mu:.2f}")
label_mu_cut = ROOT.TLatex(0.45, 0.80, f"#mu + 1#sigma (cut) = {value_mu_cut:.2f}")
label_mu.SetTextSize(0.025)
label_mu_cut.SetTextSize(0.025)
label_mu.SetNDC()
label_mu_cut.SetNDC()
label_mu.Draw()
label_mu_cut.Draw()
label_mu2 = ROOT.TLatex(0.45, 0.75, f"Events: {h_derr_mu.GetEntries():.0f} --> {h_derr_mu_cut.GetEntries():.0f}")
label_mu2.SetTextSize(0.025)
label_mu2.SetNDC()
label_mu2.Draw()

# Add a legend for the second pad
legend2 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend2.AddEntry(h_derr_mu, "uncut", "l")
legend2.AddEntry(h_derr_mu_cut, "cut", "l")
legend2.Draw()

canvas.Update()
canvas.Write()


print("errors drawn in %f seconds." % (time()-tdrawerr))


fout.Close()


print("operation took %f seconds." % (time()-t))
