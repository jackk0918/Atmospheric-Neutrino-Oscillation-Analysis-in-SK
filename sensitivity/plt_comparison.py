import ROOT

# Open the ROOT files
file_n0 = ROOT.TFile.Open("0000.root")
file_n1 = ROOT.TFile.Open("0451.root")

# Retrieve histograms
n0e = file_n0.Get("0000e")
n0mu = file_n0.Get("0000mu")
n1e = file_n1.Get("0451e")
n1mu = file_n1.Get("0451mu")

# Remove default statistic boxes
n0e.SetStats(0)
n1e.SetStats(0)
n0mu.SetStats(0)
n1mu.SetStats(0)
ROOT.gStyle.SetOptStat(0)

# Create the first plot: n0e and n1e
canvas1 = ROOT.TCanvas("elike", "Sub-GeV e-like", 400, 300)
canvas1.SetTitle("Sub-GeV e-like")

# Adjust the Y-axis range to fit both histograms (minimum and maximum)
max_bin_content1 = max(n0e.GetMaximum(), n1e.GetMaximum())
n0e.SetMinimum(0)
n0e.SetMaximum(max_bin_content1 / 0.8)

n0e.SetLineColor(ROOT.kBlue)
n1e.SetLineColor(ROOT.kRed)
n0e.SetLineWidth(2)
n1e.SetLineWidth(2)

n0e.Draw()
n1e.Draw("SAME")

# Add legend
legend1 = ROOT.TLegend(0.7, 0.75, 0.9, 0.85)
legend1.AddEntry(n0e, "no oscillation", "l")
legend1.AddEntry(n1e, "sin^{2}#theta_{23} = 0.451", "l")
legend1.Draw()


# Create the second plot: n0mu and n1mu
canvas2 = ROOT.TCanvas("mulike", "Sub-GeV #mu-like", 400, 300)
canvas2.SetTitle("Sub-GeV #mu-like")

n0mu.SetMinimum(0)
n0mu.SetMaximum(n1mu.GetMaximum() / 0.8)

n0mu.SetLineColor(ROOT.kBlue)
n1mu.SetLineColor(ROOT.kRed)
n0mu.SetLineWidth(2)
n1mu.SetLineWidth(2)

n0mu.Draw()
n1mu.Draw("SAME")

# Add legend
legend2 = ROOT.TLegend(0.7, 0.75, 0.9, 0.85)  # Adjust position as needed
legend2.AddEntry(n0mu, "no oscillation", "l")
legend2.AddEntry(n1mu, "sin^{2}#theta_{23} = 0.451", "l")
legend2.Draw()

output = ROOT.TFile("plt_comparison.root", "RECREATE")
canvas1.Write("e-like_distribution")
canvas2.Write("mu-like_distribution")
output.Close()

# Close the files
file_n0.Close()
file_n1.Close()
