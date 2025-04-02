import ROOT
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from calc_Dchi2 import delta_chi2

def calculate_chi2_values(dmsq_values, root_files):
    """Calculate chi-squared values from ROOT files"""
    delta_chi2_e = np.zeros(len(dmsq_values))
    delta_chi2_mu = np.zeros(len(dmsq_values))

    null_file = ROOT.TFile.Open("histograms/histograms_dmsq/24.root", "READ")
    h0_e = null_file.Get("e")
    h0_mu = null_file.Get("mu")

    for i, (value, filename) in enumerate(zip(dmsq_values, root_files)):
        f = ROOT.TFile.Open(f"histograms/histograms_dmsq/{filename}.root", "READ")
        h1_e = f.Get(f"e")
        h1_mu = f.Get(f"mu")

        delta_chi2_e[i] = delta_chi2(h0_e, h1_e)
        delta_chi2_mu[i] = delta_chi2(h0_mu, h1_mu)
        f.Close()

    null_file.Close()
    return delta_chi2_e, delta_chi2_mu

def find_sigma_errors(null_value, dmsq_values, chi2_values):
    """
    Compute 1σ intervals and return the errors relative to the null value.
    """
    crossings = []
    for i in range(len(dmsq_values) - 1):
        if (chi2_values[i] - 1) * (chi2_values[i + 1] - 1) <= 0:
            # Linear interpolation to find the crossing point
            x1, x2 = dmsq_values[i], dmsq_values[i + 1]
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


def plot_chi2_combined(dmsq_values, delta_chi2_sum):
    """Create a single plot for both e-like and mu-like channels with enhanced text scaling"""
    fig, ax = plt.subplots(figsize=(8, 7))

    # Colors
    teal_color = '#2AA198'
    orange_color = '#CB4B16'

    # Plot
    ax.plot(dmsq_values, delta_chi2_sum, '-', color=teal_color, linewidth=4)

    # Add 1σ line
    ax.axhline(y=1, color='black', linestyle=':', linewidth=1.5, alpha=0.9)
    ax.text(1.21, 1.5, '1σ', color='black', fontsize=16)

    # Customize plot
    ax.set_xlabel(r'$\Delta m^2_{32,31} (10^{-3} eV^2)$', fontsize=20, labelpad=10)
    ax.set_ylabel(r'$\Delta\chi^2$', fontsize=20, labelpad=10)
    ax.set_xlim(1.2, 3.6)
    ax.set_ylim(0, 16)
    ax.grid(False)
    ax.tick_params(axis='both', direction='in', labelsize=16, length=6, width=2)
    ax.set_xticks(np.arange(1.5, 3.6, 0.5))
    ax.set_xticks(np.arange(1.2, 3.7, 0.1), minor=True)
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
    dmsq_values = [1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6]
    root_files = ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36']
    null_value = 2.4

    # Calculate chi-squared values
    delta_chi2_e, delta_chi2_mu = calculate_chi2_values(dmsq_values, root_files)
    delta_chi2_sum = delta_chi2_e + delta_chi2_mu


    # Find 1σ errors relative to the null value
    low, high = find_sigma_errors(null_value, dmsq_values, delta_chi2_sum)

    # Create combined plot
    fig, ax = plot_chi2_combined(dmsq_values, delta_chi2_sum)

    # Save plot
    plt.savefig("dmsq_delta_chi2_plot.png", dpi=300, bbox_inches='tight')
    print("Plot saved as 'dmsq_delta_chi2_plot.png'")


    # Print results
    print("\n1σ errors:")
    if low is not None and high is not None:
        print(f"-{low:.3f}, +{high:.3f}")
    else:
        print("No complete 1σ interval found")



if __name__ == "__main__":
    main()