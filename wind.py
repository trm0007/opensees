import numpy as np
import math
# Table 6.2.12: Wind Directionality Factor, Kd
def wind_directionality_factor(structure_type):
    factors = {
        "Buildings": {
            "Main Wind Force Resisting System": 0.85,
            "Components and Cladding": {
                "Arched Roofs": 0.85,
                "Chimneys, Tanks, and Similar Structures": {
                    "Square": 0.90,
                    "Hexagonal": 0.95,
                    "Round": 0.95
                }
            }
        },
        "Solid Signs": 0.85,
        "Open Signs and Lattice Framework": 0.85,
        "Trussed Towers": {
            "Triangular, square, rectangular": 0.85,
            "All other cross section": 0.95
        }
    }

    def get_factor(data, structure_type):
        for key, value in data.items():
            if key == structure_type:
                return value
            elif isinstance(value, dict):
                result = get_factor(value, structure_type)
                if result is not None:
                    return result
        return None

    factor = get_factor(factors, structure_type)
    if factor is None:
        print("Structure type not found.")
        return None
    else:
        return factor
# Table 6.2.12: Wind Directionality Factor, Kd
# Example usage:
structure_type = "Main Wind Force Resisting System"
kd = wind_directionality_factor(structure_type)
if kd is not None:
    print(f"The wind directionality factor for {structure_type} is {kd}.")

# 2.4.7.2
def calculate_topographic_factor(k1,k2,k3):
    """
    Calculate the topographic factor for wind speed-up effect.

    Parameters:
        k1 (float): Value from Figure 6.2.4
        k2 (float): Value from Figure 6.2.4
        k3 (float): Value from Figure 6.2.4

    Returns:
        float: Topographic factor
        # 2.4.7.2
    """
    kzt = (1 + k1*k2*k3)**2
    return kzt
# 2.4.7.2
# Function call with sample values (replace with actual values from Figure 6.2.4)
kzt = calculate_topographic_factor(0.1, 0.2, 0.3)
print("kzt:", kzt)



# Table 6.2.8: Basic Wind Speeds, V, for Selected Locations in Bangladesh
def basic_wind_speed(location):
    # Table 6.2.8: Basic Wind Speeds, V, for Selected Locations in Bangladesh
    wind_speeds = {
        "Angarpota": 47.8, "Lalmonirhat": 63.7,
        "Bagerhat": 77.5, "Madaripur": 68.1,
        "Bandarban": 62.5, "Magura": 65.0,
        "Barguna": 80.0, "Manikganj": 58.2,
        "Barisal": 78.7, "Meherpur": 58.2,
        "Bhola": 69.5, "Maheshkhali": 80.0,
        "Bogra": 61.9, "Moulvibazar": 53.0,
        "Brahmanbaria": 56.7, "Munshiganj": 57.1,
        "Chandpur": 50.6, "Mymensingh": 67.4,
        "Chapai Nawabganj": 41.4, "Naogaon": 55.2,
        "Chittagong": 80.0, "Narail": 68.6,
        "Chuadanga": 61.9, "Narayanganj": 61.1,
        "Comilla": 61.4, "Narsinghdi": 59.7,
        "Cox’s Bazar": 80.0, "Natore": 61.9,
        "Dahagram": 47.8, "Netrokona": 65.6,
        "Dhaka": 65.7, "Nilphamari": 44.7,
        "Dinajpur": 41.4, "Noakhali": 57.1,
        "Faridpur": 63.1, "Pabna": 63.1,
        "Feni": 64.1, "Panchagarh": 41.4,
        "Gaibandha": 65.6, "Patuakhali": 80.0,
        "Gazipur": 66.5, "Pirojpur": 80.0,
        "Gopalganj": 74.5, "Rajbari": 59.1,
        "Habiganj": 54.2, "Rajshahi": 49.2,
        "Hatiya": 80.0, "Rangamati": 56.7,
        "Ishurdi": 69.5, "Rangpur": 65.3,
        "Joypurhat": 56.7, "Satkhira": 57.6,
        "Jamalpur": 56.7, "Shariatpur": 61.9,
        "Jessore": 64.1, "Sherpur": 62.5,
        "Jhalakati": 80.0, "Sirajganj": 50.6,
        "Jhenaidah": 65.0, "Srimangal": 50.6,
        "Khagrachhari": 56.7, "St. Martin’s Island": 80.0,
        "Khulna": 73.3, "Sunamganj": 61.1,
        "Kutubdia": 80.0, "Sylhet": 61.1,
        "Kishoreganj": 64.7, "Sandwip": 80.0,
        "Kurigram": 65.6, "Tangail": 50.6,
        "Kushtia": 66.9, "Teknaf": 80.0,
        "Lakshmipur": 51.2, "Thakurgaon": 41.4
    }

    # Check if the location is in the wind_speeds dictionary
    if location in wind_speeds:
        return wind_speeds[location]
    else:
        print("Location not found.")
        return None
# Table 6.2.8: Basic Wind Speeds, V, for Selected Locations in Bangladesh
# Example usage:
location = "Dhaka"
wind_speed = basic_wind_speed(location)
if wind_speed is not None:
    print(f"The basic wind speed for {location} is {wind_speed} m/s.")

# Table 6.2.9: Importance Factor, I (Wind Loads)
def importance_factor(occupancy_category, basic_wind_speed):
    # Table 6.2.9: Importance Factor, I (Wind Loads)
    # Importance factors for different occupancy categories based on basic wind speed ranges
    importance_factors = {
        "I": {"Non-Cyclone Prone Regions": 0.87, "Cyclone Prone Regions (V 5 38-44 m/s)": 0.87,
              "Cyclone Prone Regions (V > 44 m/s)": 0.77},
        "II": {"Non-Cyclone Prone Regions": 1.0, "Cyclone Prone Regions (V 5 38-44 m/s)": 1.0,
               "Cyclone Prone Regions (V > 44 m/s)": 1.0},
        "III": {"Non-Cyclone Prone Regions": 1.15, "Cyclone Prone Regions (V 5 38-44 m/s)": 1.15,
                "Cyclone Prone Regions (V > 44 m/s)": 1.15},
        "IV": {"Non-Cyclone Prone Regions": 1.15, "Cyclone Prone Regions (V 5 38-44 m/s)": 1.15,
               "Cyclone Prone Regions (V > 44 m/s)": 1.15}
    }

    # Check if the occupancy category is in the importance_factors dictionary
    if occupancy_category in importance_factors:
        if basic_wind_speed >= 44:
            cyclone_prone_region_key = "Cyclone Prone Regions (V > 44 m/s)"
        elif 38 <= basic_wind_speed <= 44:
            cyclone_prone_region_key = "Cyclone Prone Regions (V 5 38-44 m/s)"
        else:
            cyclone_prone_region_key = "Non-Cyclone Prone Regions"

        return importance_factors[occupancy_category][cyclone_prone_region_key]
    else:
        print( "Occupancy category not found." )
        return None

# Table 6.2.9: Importance Factor, I (Wind Loads)
# Example usage:
occupancy_category = "III"
basic_wind_speed = 50  # example basic wind speed in m/s
importance_factor = importance_factor( occupancy_category, basic_wind_speed )
if importance_factor is not None:
    print(
        f"The importance factor for Occupancy Category {occupancy_category} and basic wind speed {basic_wind_speed} m/s is {importance_factor}." )
# Table 6.2.10: Terrain Exposure Constants

M = 0.000613 * kzt * kd * wind_speed**2 * importance_factor
print("M=", {M})
def terrain_exposure_constants(exposure_category):
    """
    Retrieve terrain exposure constants for the specified exposure category.

    Args:
    exposure_category (str): The exposure category (A, B, or C).

    Returns:
    dict: Dictionary containing the terrain exposure constants.
    """
    # Define the table as a dictionary
    # Table 6.2.10: Terrain Exposure Constants
    constants = {
        "A": {
            "exposure": "A",
            "alpha": 7.0,           # ¯ˆ˘ (m): Height of the terrain (m)
            "zg": 365.76,           # ¸˝ ˛_: Terrain roughness (m)
            "a": 1/7,               # ¸_: Constant 'a' for terrain exposure
            "b": 0.84,              # c  (m): Constant 'b' for terrain exposure
            "a_var": 1/4.0,         # _ ˆ (m): Constant 'a_var' for terrain exposure
            "b_var": 0.45,          # ˆ: Constant 'b_var' for terrain exposure
            "c": 0.30,              # ˆ: Constant 'c' for terrain exposure
            "L": 97.54,             # L: Maximum fetch (m)
            "epsilon": 1/3.0,       # 1/α: Ratio of terrain roughness to structure height
            "z_min": 9.14           # z_min: Minimum height used to ensure that the equivalent height z is greater of 0.6*h or z_min
        },
        "B": {
            "exposure": "B",
            "alpha": 9.5,
            "zg": 274.32,
            "a": 1/9.5,
            "b": 1.00,
            "a_var": 1/6.5,
            "b_var": 0.65,
            "c": 0.20,
            "L": 152.4,
            "epsilon": 1/5.0,
            "z_min": 4.57
        },
        "C": {
            "exposure": "C",
            "alpha": 11.5,
            "zg": 213.36,
            "a": 1/11.5,
            "b": 1.07,
            "a_var": 1/9.0,
            "b_var": 0.80,
            "c": 0.15,
            "L": 198.12,
            "epsilon": 1/8.0,
            "z_min": 2.13
        }
    }

    # Check if the exposure category is in the constants dictionary
    if exposure_category in constants:
        return constants[exposure_category]
    else:
        print("Exposure category not found.")
        return None
# Table 6.2.10: Terrain Exposure Constants
# Example usage:
exposure_category = "A"
constants = terrain_exposure_constants(exposure_category)
# Assuming constants is a dictionary containing the provided key-value pairs

# Assigning float values to variables using constants.get() method
exposure = (constants.get("exposure"))
alpha = float(constants.get("alpha"))
zg = float(constants.get("zg"))
a = float(constants.get("a"))
b = float(constants.get("b"))
a_var = float(constants.get("a_var"))
b_var = float(constants.get("b_var"))
c = float(constants.get("c"))
l = float(constants.get("L"))
epsilon = float(constants.get("epsilon"))
z_min = float(constants.get("z_min"))
# exposure,alpha,zg,a,b,a_var,b_var,c,L,epsilon,z_min

if constants is not None:
    print(f"Terrain exposure constants for category {exposure_category}:")
    for key, value in constants.items():
        print(f"{key}: {value}")

def velocity_pressure_coefficient(z, zg, alpha):
    """
    Calculate the velocity pressure exposure coefficient Kz.

    Args:
    z (float): Height above ground level (m).
    zg (float): Height of the structure above ground level (m).
    alpha (float): Coefficient depending on exposure category.

    Returns:
    float: Velocity pressure exposure coefficient Kz.
    """
    if z >= 4.57:  # For z >= 4.57 m
        kz = 2.01 * (z / zg) ** (2 / alpha)
        return kz
    elif z < 4.57 and z >= 9.1:  # For z < 4.57 m and z >= 9.1 m
        kz = 2.01 * (4.57 / zg) ** (2 / alpha)
        return kz
    else:  # For z < 9.1 m
        print("Warning: z should not be less than 9.1 m for Case 1 in exposure A.")
        return None

# Example usage:
z = 5  # Height above ground level in meters
Kz = velocity_pressure_coefficient(z, zg, alpha)
if Kz is not None:
    print(f"The velocity pressure exposure coefficient Kz is {Kz}.")

# Compute the Time Period of the Structure, T from Teismic Analysis
T = 1.23
frequency = 1/T
# exposure,alpha,zg,a,b,a_var,b_var,c,L,epsilon,z_min

def Compute_Gust_factor(T,z_min,c,l,epsilon,wind_speed,a_var,b_var,h,L,B):
    """
    Computes various parameters for flexible buildings or structures based on provided inputs.

    Parameters:
    n (float): Fundamental frequency of the structure.
    Vkph (float): Velocity in kilometers per hour.

    Returns:
    tuple: A tuple containing computed values for V2, R, RRB, RL, R9, I, Q, and G.
    """

    # Constants and initial values
    n1 = 1/T # n1 = building natural frequency
    B_damping_percentage = 0.01  # 1% for steel and 2% for concrete
    V = 0.2778 * wind_speed
    z_var = np.maximum(0.6*h, z_min)
    print(f'z_var={z_var}')

    # 1. Compute V2
    Vz = b_var * (z_var / 10)**a_var * V
    print("Vz:", Vz)

    # 2. Compute Rh, Rb, and Rl
    ah = (4.6 * n1 * h / Vz)
    ab = (4.6 * n1 * B / Vz)
    al = (15.4 * n1 * L / Vz)
    nh = (ah) ** -1
    nb = (ab) ** -1
    nl = (al) ** -1
    Rh = np.maximum( nh - 0.5 * nh * nh * (1 - np.exp( -2 * ah )), 0 )
    Rb = np.maximum( nb - 0.5 * nb * nb * (1 - np.exp( -2 * ab )), 0 )
    Rl = np.maximum( nl - 0.5 * nl * nl * (1 - np.exp( -2 * al )), 0 )

    print("Rh:", Rh)
    print("Rb:", Rb)
    print("Rl:", Rl)

    # 3. Compute Resonant Response Factor (Rn)
    # Assuming a value for z, change according to requirement
    Lz = l * (z_var / 10)**epsilon
    print( "Lz:", Lz )
    N1 = n1 * Lz / Vz
    print( "N1:", N1 )
    Rn = (7.47 * N1) / (1+10.3 * N1)**(5/3)
    print("Rn:", Rn)

    # 4. Compute Resonant Response Factor (R)
    R = math.sqrt((1/B_damping_percentage)*Rn*Rh*Rb*(0.53+0.47*Rl))
    print("R:", R)

    # 5. Compute gR
    gR = np.sqrt(2*np.log(3600*n1)) + 0.577/np.sqrt(2*np.log(3600*n1))
    print("gR:", gR)
    gQ = 3.4
    gv = 3.4

    # 6. Compute I and Q:
    Iz = c*(10/z)**(1/6)
    print( "Iz:", Iz )
    Q = np.sqrt(1/(1+0.63*((B+h)/Lz)**0.63))
    print("Q:", Q)

    # 6. Gust Factor, Gf:
    if T > 1:
        Gf = 0.925 * ((1+1.7*Iz*np.sqrt(gQ**2+gR**2*R**2))/(1+1.7*gv*Iz))
        print( "Gust Factor, Gf:", Gf )
    # Return computed values
        return Gf
    else:
        asr = 1 + 1.7 * gQ * Iz * Q
        bsr = 1 + 1.7 * gv * Iz
        Gf = np.maximum(0.925 * (asr / bsr),0.85)
        print( "Gust Factor, Gf:", Gf )
        # Return computed values
        return Gf


print("successful")
h,L,B = 40, 20, 50
Compute_Gust_factor(T,z_min,c,l,epsilon,wind_speed,a_var,b_var,h,L,B)


def calculate_wall_pressure_coefficient(surface, L_over_B):
    """
    Calculate the wall pressure coefficient (Cp) based on surface type and dimensions.

    Parameters:
        surface (str): Type of wall surface (e.g., 'Windward Wall', 'Leeward Wall', 'Side Wall')
        L_over_B (float): Ratio of length to breadth
        q (float): Velocity pressure (q)

    Returns:
        float: Wall pressure coefficient (Cp)
    """
    
    if surface == 'Windward Wall':
        Cp = 0.8
    elif surface == 'Leeward Wall':
        if L_over_B <= 1:
            Cp = -0.5
        elif L_over_B <= 2:
            Cp = -0.3
        else:
            # Perform linear interpolation for L_over_B between 2 and infinity
            Cp = -0.3 + (-0.2 + (-0.3)) * (L_over_B - 2) / (L_over_B - 2 + 1)

    elif surface == 'Side Wall':
        Cp = -0.7
    else:
        Cp = None  # Unknown surface type

    return Cp


# Function call with sample values
surface = 'Leeward Wall'
L_over_B = 0.5  # Example value for L/B ratio
q = 10  # Sample velocity pressure (q)

Cp = calculate_wall_pressure_coefficient( surface, L_over_B)
if Cp is not None:
    print( f"Wall pressure coefficient (Cp) for {surface} (L/B={L_over_B}): {Cp}" )
else:
    print( "Invalid surface type. Please specify 'Windward Wall', 'Leeward Wall', or 'Side Wall'." )
