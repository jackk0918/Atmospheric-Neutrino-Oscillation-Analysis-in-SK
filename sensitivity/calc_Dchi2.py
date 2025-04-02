import ROOT
import numpy as np

def delta_chi2(h0, h1):
    """
    Calculate the Delta Chi-Square value between two histograms.

    Args: TH1D
    Returns: float (Delta Chi-Square value)
    """

    Dchi2 = 0.0
    nbins = h0.GetNbinsX()

    for i in range(1, nbins + 1):
        N0 = h0.GetBinContent(i)  # Null hypothesis event rate
        N1 = h1.GetBinContent(i)   # Alternative hypothesis event rate

        if N0 > 0 and N1 > 0:  # Avoid division by zero or log of zero
            Dchi2 += 2 * (N0 - N1 + N1 * np.log(N1 / N0))

    return Dchi2


def main():
    # Load ROOT files
    file0 = ROOT.TFile.Open("0000.root")
    file1 = ROOT.TFile.Open("0451.root")

    # Load histograms
    n0e = file0.Get("0000e")
    n0mu = file0.Get("0000mu")
    n1e = file1.Get("0451e")
    n1mu = file1.Get("0451mu")

    if not n0e or not n1e or not n0mu or not n1mu:
        print("Error: One or more histograms could not be found in the file.")
        return

    # Calculate Delta Chi-Square for each category
    Dchi2e = delta_chi2(n0e, n1e)
    Dchi2mu = delta_chi2(n0mu, n1mu)

    # Calculate sensitivities
    sensitivity_e = np.sqrt(Dchi2e)
    sensitivity_mu = np.sqrt(Dchi2mu)

    # Print results
    print(f"Electron-like sensitivity (sigma): {sensitivity_e:.2f}")
    print(f"Muon-like sensitivity (sigma): {sensitivity_mu:.2f}")

    file0.Close()
    file1.Close()

if __name__ == "__main__":
    main()
