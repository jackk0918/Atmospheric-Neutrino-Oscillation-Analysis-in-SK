# Plot of the ratio of weighted and unweighted histograms

import ROOT

fout = ROOT.TFile("h_ratios.root", "RECREATE")

file1 = ROOT.TFile.Open("h_recon.root", "READ")
h_e_weighted = file1.Get("e_reduced")
h_mu_weighted = file1.Get("mu_reduced")

file2 = ROOT.TFile.Open("h_unweighted.root", "READ")
h_e_unweighted = file2.Get("e_recon")
h_mu_unweighted = file2.Get("mu_recon")

h_e_ratio = h_e_weighted.Clone("e_ratio")
h_e_ratio.SetTitle("e with / witout oscillation")
h_e_ratio.Divide(h_e_unweighted)

h_mu_ratio = h_mu_weighted.Clone("mu_ratio")
h_mu_ratio.SetTitle("mu with / without oscillation")
h_mu_ratio.Divide(h_mu_unweighted)

fout.cd()

h_e_ratio.Write()
h_mu_ratio.Write()


fout.Close()