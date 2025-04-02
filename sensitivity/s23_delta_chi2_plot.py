import ROOT
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from chisq import delta_chi2

def calculate_chi2_values(s23sq_values, root_files):
    """Calculate chi-squared values from ROOT files"""
    delta_chi2_e = np.zeros(len(s23sq_values))
    delta_chi2_mu = np.zeros(len(s23sq_values))

    null_file = ROOT.TFile.Open("histograms/0451.root", "READ")
    h0_e = null_file.Get("0451e")
    h0_mu = null_file.Get("0451mu")

    for i, (value, filename) in enumerate(zip(s23sq_values, root_files)):
        f = ROOT.TFile.Open(f"histograms/{filename}.root", "READ")
        h1_e = f.Get(f"{filename}e")
        h1_mu = f.Get(f"{filename}mu")

        delta_chi2_e[i] = delta_chi2(h0_e, h1_e)
        delta_chi2_mu[i] = delta_chi2(h0_mu, h1_mu)
        f.Close()

    null_file.Close()
    return delta_chi2_e, delta_chi2_mu

def find_sigma_errors(null_value, s23sq_values, chi2_values):
    """
    Compute 1σ intervals and return the errors relative to the null value.
    """
    crossings = []
    for i in range(len(s23sq_values) - 1):
        if (chi2_values[i] - 1) * (chi2_values[i + 1] - 1) <= 0:
            # Linear interpolation to find the crossing point
            x1, x2 = s23sq_values[i], s23sq_values[i + 1]
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


def plot_chi2_combined(s23sq_values, delta_chi2_sum):
    """Create a single plot for both e-like and mu-like channels with enhanced text scaling"""
    fig, ax = plt.subplots(figsize=(8, 7))  # Adjusted aspect ratio

    # Colors
    teal_color = '#2AA198'
    orange_color = '#CB4B16'

    # Plot
    ax.plot(s23sq_values, delta_chi2_sum, '-', color=teal_color, linewidth=4)

    # Add 1σ line
    ax.axhline(y=1, color='black', linestyle=':', linewidth=1.5, alpha=0.9)
    ax.text(0.31, 1.5, '1σ', color='black', fontsize=16)

    # Customize plot
    ax.set_xlabel(r'$\sin^2\theta_{23}$', fontsize=20, labelpad=10)
    ax.set_ylabel(r'$\Delta\chi^2$', fontsize=20, labelpad=10)
    ax.set_xlim(0.3, 0.775)
    ax.set_ylim(0, 16)
    ax.grid(False)
    ax.tick_params(axis='both', direction='in', labelsize=16, length=6, width=2)
    ax.set_xticks(np.arange(0.3, 0.776, 0.05))
    ax.set_xticks(np.arange(0.3, 0.776, 0.01), minor=True)
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
    #s23sq_values = [0.300, 0.353, 0.406, 0.451, 0.511, 0.564, 0.617, 0.669, 0.722, 0.775]
    #root_files = ['0300', '0353', '0406', '0451', '0511', '0564', '0617', '0669', '0722', '0775']
    #s23sq_values = [0.300, 0.323, 0.346, 0.368, 0.391, 0.414, 0.437, 0.451, 0.482, 0.505, 0.528, 0.551, 0.574, 0.596, 0.619, 0.642, 0.665, 0.688, 0.710, 0.775]
    #root_files = ['0300', '0323', '0346', '0368', '0391', '0414', '0437', '0451', '0482', '0505', '0528', '0551', '0574', '0596', '0619', '0642', '0665', '0688', '0710', '0775']
    s23sq_values = [0.300, 0.326, 0.353, 0.379, 0.406, 0.430, 0.451, 0.479, 0.511, 0.538, 0.564, 0.589, 0.617, 0.643, 0.669, 0.695, 0.722, 0.748, 0.775]
    root_files = ['0300', '0326', '0353', '0379', '0406', '0430', '0451', '0479', '0511', '0538', '0564', '0589', '0617', '0643', '0669', '0695', '0722', '0748', '0775']
    null_value = 0.451

    # Calculate chi-squared values
    delta_chi2_e, delta_chi2_mu = calculate_chi2_values(s23sq_values, root_files)
    delta_chi2_sum = delta_chi2_e + delta_chi2_mu


    # Find 1σ errors relative to the null value
    low, high = find_sigma_errors(null_value, s23sq_values, delta_chi2_sum)

    # Create combined plot
    fig, ax = plot_chi2_combined(s23sq_values, delta_chi2_sum)

    # Save plot
    plt.savefig("NEWdelta_chi2_allE_plot_20.png", dpi=300, bbox_inches='tight')
    print("Plot saved as 'NEWdelta_chi2_allE_plot_20.png'")


    # Print results
    print("\n1σ errors:")
    if low is not None and high is not None:
        print(f"-{low:.3f}, +{high:.3f}")
    else:
        print("No complete 1σ interval found")



if __name__ == "__main__":
    main()