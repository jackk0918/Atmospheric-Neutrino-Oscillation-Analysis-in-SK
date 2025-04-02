# this one calculates the oscillation probability weight for single zenith and energy value.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

#Constants
radius = 6371 # km
height = 15 # km

# NuFast constants
eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4
YerhoE2a = 1.52e-4
Ye = 0.5
N_Newton = 0

# ---------------------------------------------------------------
# OSCILLATION PARAMETERS (PDG 2023)
s12sq = 0.307
s13sq = 0.020           # 0.020
s23sq = 0.451           # 0.451 normal hierarchy <-------------------------------
delta = 1.23 * np.pi
Dmsq21 = 7.53e-5        # eV^2
Dmsq31 = 2.40e-3        # 2.40e-3 eV^2
# ---------------------------------------------------------------

densities = []
totlengths = []

# Total path length calculation
def length(z):
    return -radius * np.cos(z) + np.sqrt(radius**2 * np.cos(z)**2 + height**2 + 2 * radius * height)

# The density model
def density(z):

    regions = [
        {"name": "inner", "r_min": 0, "r_max": 1220, "density": 13},
        {"name": "outer", "r_min": 1220, "r_max": 3480, "density": 11.3},
        {"name": "mantle", "r_min": 3480, "r_max": 5701, "density": 5},
        {"name": "crust", "r_min": 5701, "r_max": 6371, "density": 3}
    ]

    # Calculate the path length through each region
    path_lengths = []
    for region in regions:
        #print()
        #print(f"region: {region['name']}")
        #print(f"r sin z: {round(radius * np.sin(z), 1)}, r_max: {region['r_max']}")
        if radius * np.sin(z) >= region["r_max"] or z <= np.pi/2:   # condition: does the path NOT cross the region OR cosine angle <= pi/2?
            path_lengths.append(0)  # if yes, then set path_length = 0
            #print("path_length: 0")
            continue

        outer_path_length = 2 * np.sqrt(region["r_max"]**2 - (radius * np.sin(z))**2)   # chord length formula a = 2√r^2 - d^2
        #print(f"outer_path_length: {round(outer_path_length, 0)}")

        if region["r_min"] > abs(radius * np.sin(z)):
            inner_path_length = 2 * np.sqrt(region["r_min"]**2 - (radius * np.sin(z))**2)
            path_length = outer_path_length - inner_path_length
            #print(f"inner_path_length: {round(inner_path_length, 0)}")
        else:
            path_length = outer_path_length
            #print(f"inner_path_length: 0")

        path_lengths.append(path_length)
        #print(f"path_length: {round(path_length, 0)}")

    # Calculate the total path length
    total_path_length = sum(path_lengths)
    totlengths.append(round(total_path_length, 0))

    # Calculate the weighted average density
    if total_path_length == 0:
        density = 0
    else:
        density = sum(region["density"] * path_length for region, path_length in zip(regions, path_lengths)) / total_path_length
    densities.append(round(density, 2))
    return density

# Calculate the transition probability (NuFast algorithm)
def Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton):
    # --------------------------------------------------------------------- #
    # First calculate useful simple functions of the oscillation parameters #
    # --------------------------------------------------------------------- #
    c13sq = 1 - s13sq

    # Ueisq's
    Ue2sq = c13sq * s12sq
    Ue3sq = s13sq

    # Umisq's, Utisq's and Jvac
    Um3sq = c13sq * s23sq
    # Um2sq and Ut2sq are used here as temporary variables, will be properly defined later
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)

    Jrr = np.sqrt(Um2sq * Ut2sq)
    sind = np.sin(delta)
    cosd = np.cos(delta)

    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jmatter = 8 * Jrr * c13sq * sind
    Amatter = Ye * rho * E * YerhoE2a
    Dmsqee = Dmsq31 - s12sq * Dmsq21

    # calculate A, B, C, See, Tee, and part of Tmm
    A = Dmsq21 + Dmsq31 # temporary variable
    See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq
    Tmm = Dmsq21 * Dmsq31 # using Tmm as a temporary variable
    Tee = Tmm * (1 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A + Amatter

    # ---------------------------------- #
    # Get lambda3 from lambda+ of MP/DMP #
    # ---------------------------------- #
    xmat = Amatter / Dmsqee
    tmp = 1 - xmat
    lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + np.sqrt(tmp * tmp + 4 * s13sq * xmat))

    # ---------------------------------------------------------------------------- #
    # Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) #
    # ---------------------------------------------------------------------------- #
    B = Tmm + Amatter * See # B is only needed for N_Newton >= 1
    for i in range(N_Newton):
        lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B) # this strange form prefers additions to multiplications

    # ------------------- #
    # Get  Delta lambda's #
    # ------------------- #
    tmp = A - lambda3

    Dlambda21 = np.sqrt(tmp * tmp - 4 * C / lambda3)
    lambda2 = 0.5 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21
    # ----------------------- #
    # Use Rosetta for Veisq's #
    # ----------------------- #
    # denominators
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31

    # numerators
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2

    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq
    Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A)

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2

    # ------------- #
    # Use NHS for J #
    # ------------- #
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv

    # ----------------------- #
    # Get all elements of Usq #
    # ----------------------- #
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq

    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq

    # ----------------------- #
    # Get the kinematic terms #
    # ----------------------- #
    Lover4E = eVsqkm_to_GeV_over4 * L / E

    D21 = Dlambda21 * Lover4E
    D32 = Dlambda32 * Lover4E

    sinD21 = np.sin(D21)
    sinD31 = np.sin(D32 + D21)
    sinD32 = np.sin(D32)
    triple_sin = sinD21 * sinD31 * sinD32

    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32

    # ------------------------------------------------------------------- #
    # Calculate the three necessary probabilities, separating CPC and CPV #
    # ------------------------------------------------------------------- #
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 \
                + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 \
                + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    Pme_CPV = -Jmatter * triple_sin

    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 \
                 + Um3sq * Um1sq * sinsqD31_2 \
                 + Um3sq * Um2sq * sinsqD32_2)

    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 \
                 + Ue3sq * Ue1sq * sinsqD31_2 \
                 + Ue3sq * Ue2sq * sinsqD32_2)

    probs_me = Pme_CPC + Pme_CPV
    probs_mm = Pmm
    probs_ee = Pee
    probs_em = Pme_CPC - Pme_CPV

    return [probs_me, probs_mm, probs_ee, probs_em]

# Check function (answer should be around 0.5):
#test = Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, length(2.5), 6, density(2.5), Ye, N_Newton)
#print("this is a test:")
#print(test)
#exit()

def calcNuFast(zenith, energy, particle):

    # Values
    z_value = np.arccos(zenith)

    # Calculate p_matter for each combination of z and energy
    if particle == 12:
        p_value1 = Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, length(z_value), energy, density(z_value), Ye, N_Newton)[0]
        p_value2 = Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, length(z_value), energy, density(z_value), Ye, N_Newton)[2]
    elif particle == 14:
        p_value1 = Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, length(z_value), energy, density(z_value), Ye, N_Newton)[1]
        p_value2 = Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, length(z_value), energy, density(z_value), Ye, N_Newton)[3]

    transitioned_nu = p_value1 * 2
    survived_nu = p_value2
    weight = transitioned_nu + survived_nu

    # Normalise
    #max_val = np.max(weight)
    #weight /= max_val

    return weight