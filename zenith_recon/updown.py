# Calculate the up-down asymmetry.

import ROOT


def updown(hist):
    up_sum = 0
    down_sum = 0
    nbins = hist.GetNbinsX()

    for i in range(1, nbins + 1):
        # Get bin center and content
        bin_center = hist.GetBinCenter(i)
        bin_content = hist.GetBinContent(i)

        # Sum up content based on bin center
        if bin_center < 0:
            up_sum += bin_content
        elif bin_center > 0:
            down_sum += bin_content
        # Bins exactly at 0 are ignored

    if down_sum == 0:
        print("Warning: Division by zero. Cannot calculate up-down ratio.")
        return None

    ratio = up_sum / down_sum
    return ratio


file = ROOT.TFile.Open("0000.root")
h_e = file.Get("0000e")
h_mu = file.Get("0000mu")

print("Up/down e = ", updown(h_e))
print("Up/down mu = ", updown(h_mu))