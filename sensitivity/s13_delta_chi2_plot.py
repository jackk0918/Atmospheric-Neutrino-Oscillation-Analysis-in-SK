import ROOT
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from calc_Dchi2 import delta_chi2

def calculate_chi2_values(s13sq_values, root_files):
    """Calculate chi-squared values from ROOT files"""
    delta_chi2_e = np.zeros(len(s13sq_values))
    delta_chi2_mu = np.zeros(len(s13sq_values))

    null_file = ROOT.TFile.Open("histograms/histograms_theta13/020.root", "READ")
    h0_e = null_file.Get("020e")
    h0_mu = null_file.Get("020mu")

    for i, (value, filename) in enumerate(zip(s13sq_values, root_files)):
        f = ROOT.TFile.Open(f"histograms/histograms_theta13/{filename}.root", "READ")
        h1_e = f.Get(f"{filename}e")
        h1_mu = f.Get(f"{filename}mu")

        delta_chi2_e[i] = delta_chi2(h0_e, h1_e)
        delta_chi2_mu[i] = delta_chi2(h0_mu, h1_mu)
        f.Close()

    null_file.Close()
    return delta_chi2_e, delta_chi2_mu

def find_sigma_errors(null_value, s13sq_values, chi2_values):
    """
    Compute 1σ intervals and return the errors relative to the null value.
    """
    crossings = []
    for i in range(len(s13sq_values) - 1):
        if (chi2_values[i] - 1) * (chi2_values[i + 1] - 1) <= 0:
            # Linear interpolation to find the crossing point
            x1, x2 = s13sq_values[i], s13sq_values[i + 1]
            y1, y2 = chi2_values[i], chi2_values[i + 1]
            if y1 != y2:
                crossing = x1 + (x2 - x1) * (1 - y1) / (y2 - y1)
                crossings.append(crossing)

    if len(crossings) >= 2:
        low_error = null_value - crossings[0]
        high_error = crossings[-1] - null_value
        return low_error, high_error
    else:
        return None, None  # In case no valid interval is found


def plot_chi2_combined(s13sq_values, delta_chi2_sum):
    """Create a single plot for both e-like and mu-like channels with enhanced text scaling"""
    fig, ax = plt.subplots(figsize=(8, 7))

    # Colors
    teal_color = '#2AA198'
    orange_color = '#CB4B16'

    # Plot
    ax.plot(s13sq_values, delta_chi2_sum, '-', color=teal_color, linewidth=4)

    # Add 1σ line
    ax.axhline(y=1, color='black', linestyle=':', linewidth=1.5, alpha=0.9)
    ax.text(0.001, 1.5, '1σ', color='black', fontsize=16)

    # Customize plot
    ax.set_xlabel(r'$\sin^2\theta_{13}$', fontsize=20, labelpad=10)
    ax.set_ylabel(r'$\Delta\chi^2$', fontsize=20, labelpad=10)
    ax.set_xlim(0.0, 0.075)
    ax.set_ylim(0, 16)
    ax.grid(False)
    ax.tick_params(axis='both', direction='in', labelsize=16, length=6, width=2)
    ax.set_xticks(np.arange(0, 0.076, 0.01))
    ax.set_xticks(np.arange(0, 0.076, 0.002), minor=True)
    ax.set_yticks(np.arange(0, 16, 0.5), minor=True)
    ax.tick_params(axis='both', which='minor', length=3)

    # Thicken axes
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)

    plt.tight_layout()
    return fig, ax


def main():
    # Define input values
    s13sq_values = [0.000,0.005,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075]
    root_files = ['000', '005', '010', '015', '020', '025', '030', '035', '040', '045', '050', '055', '060', '065', '070', '075']
    null_value = 0.020

    # Calculate chi-squared values
    delta_chi2_e, delta_chi2_mu = calculate_chi2_values(s13sq_values, root_files)
    delta_chi2_sum = delta_chi2_e + delta_chi2_mu


    # Find 1σ errors relative to the null value
    low, high = find_sigma_errors(null_value, s13sq_values, delta_chi2_sum)

    # Create combined plot
    fig, ax = plot_chi2_combined(s13sq_values, delta_chi2_sum)

    # Save plot
    plt.savefig("s13_delta_chi2_plot.png", dpi=300, bbox_inches='tight')
    print("Plot saved as 's13_delta_chi2_plot.png'")


    # Print results
    print("\n1σ errors:")
    if low is not None and high is not None:
        print(f"-{low:.3f}, +{high:.3f}")
    else:
        print("No complete 1σ interval found")



if __name__ == "__main__":
    main()