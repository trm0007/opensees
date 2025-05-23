import numpy as np


def calculate_N(d, N):
    total = 0
    a = 0
    b = 0
    for i in range(len(d)):
        a += (d[i] / N[i])
        b += d[i]
        total = b/a
    return total


# Table 6.2.20: Values for Coefficients to Estimate Approximate Period
def calculate_base_shear_and_overturning_moment(type):
    """
    Calculate design base shear and overturning moment for a building.

    Table 6.2.20: Values for Coefficients to Estimate Approximate Period
    Structure type       Ct       m
    Concrete moment-resisting frames    0.0466    0.9
    Steel moment-resisting frames       0.0724    0.8
    Eccentrically braced steel frame    0.0731    0.75
    All other structural systems       0.0488    0.75

    Args:
    - num_stories (int): Number of stories in the building.
    - floor_height (float): Height of each floor in meters.

    Returns:
    - base_shear (float): Design base shear in kN.
    - overturning_moment (float): Overturning moment in kNm.
    """

    # Table 6.2.20: Values for Coefficients
    Ct_values = {
        "Concrete moment-resisting frames": 0.0466,
        "Steel moment-resisting frames": 0.0724,
        "Eccentrically braced steel frame": 0.0731,
        "All other structural systems": 0.0488
    }
    m_values = {
        "Concrete moment-resisting frames": 0.9,
        "Steel moment-resisting frames": 0.8,
        "Eccentrically braced steel frame": 0.75,
        "All other structural systems": 0.75
    }

    # Using Concrete moment-resisting frames values
    Ct = Ct_values[type]
    m = m_values[type]


    return Ct, m
# Table 6.2.14: Description of Seismic Zones
# def get_seismic_zone_coefficient(location):
#     """
#     Get the seismic zone coefficient for a given location.
#
#     Table 6.2.14: Description of Seismic Zones
#     Seismic Zone   Location                                      Seismic Intensity    Seismic Zone Coefficient, Z
#     1              Southwestern part including Barisal, Khulna, Low                 0.12
#                    Jessore, Rajshahi
#     2              Lower Central and Northwestern part including Moderate             0.20
#                    Noakhali, Dhaka, Pabna, Dinajpur, as well as
#                    Southwestern corner including Sundarbans
#     3              Upper Central and Northwestern part including Severe              0.28
#                    Brahmanbaria, Sirajganj, Rangpur
#     4              Northeastern part including Sylhet,          Very Severe         0.36
#                    Mymensingh, Kurigram
#
#     Args:
#     - location (str): Location for which seismic zone coefficient is required.
#
#     Returns:
#     - zone_coefficient (float): Seismic zone coefficient (Z).
#     """
#     # Table 6.2.14: Description of Seismic Zones
#     zone_coefficients = {
#         1: 0.12,
#         2: 0.20,
#         3: 0.28,
#         4: 0.36
#     }
#
#     # # Find the seismic zone coefficient for the given location
#     # if "Barisal" in location or "Khulna" in location or "Jessore" in location or "Rajshahi" in location:
#     #     return zone_coefficients[1]
#     # elif "Noakhali" in location or "Dhaka" in location or "Pabna" in location or "Dinajpur" in location or "Sundarbans" in location:
#     #     return zone_coefficients[2]
#     # elif "Brahmanbaria" in location or "Sirajganj" in location or "Rangpur" in location:
#     #     return zone_coefficients[3]
#     # elif "Sylhet" in location or "Mymensingh" in location or "Kurigram" in location:
#     #     return zone_coefficients[4]
#     # else:
#     #     return None  # Location not found in the table
#
#     # Find the seismic zone coefficient for the given location
#     for zone, locations in {
#         1: ["Barisal", "Khulna", "Jessore", "Rajshahi"],
#         2: ["Noakhali", "Dhaka", "Pabna", "Dinajpur", "Sundarbans"],
#         3: ["Brahmanbaria", "Sirajganj", "Rangpur"],
#         4: ["Sylhet", "Mymensingh", "Kurigram"]
#     }.items():
#         if any( loc in location for loc in locations ):
#             return zone
#
#     return None  # Location not found in the table

# Table 6.2.15: Seismic Zone Coefficient Z for Some Important Towns of Bangladesh
def get_seismic_zone_coefficient_town(town):
    """
    Get the seismic zone coefficient for a given town.

    Table 6.2.15: Seismic Zone Coefficient Z for Some Important Towns of Bangladesh
    Town         Z     Town         Z     Town         Z     Town         Z
    Bagerhat   0.12   Gaibandha  0.28   Magura     0.12   Patuakhali 0.12
    Bandarban  0.28   Gazipur    0.20   Manikganj  0.20   Pirojpur    0.12
    Barguna    0.12   Gopalganj  0.12   Maulvibazar0.36   Rajbari     0.20
    Barisal    0.12   Habiganj   0.36   Meherpur   0.12   Rajshahi    0.12
    Bhola      0.12   Jaipurhat  0.20   Mongla     0.12   Rangamati   0.28
    Bogra      0.28   Jamalpur   0.36   Munshiganj 0.20   Rangpur     0.28
    Brahmanbaria0.28   Jessore    0.12   Mymensingh 0.36   Satkhira    0.12
    Chandpur   0.20   Jhalokati  0.12   Narail     0.12   Shariatpur  0.20
    Chapainababganj 0.12   Jhenaidah 0.12  Narayanganj0.20   Sherpur     0.36
    Chittagong 0.28   Khagrachari0.28   Narsingdi  0.28   Sirajganj   0.28
    Chuadanga  0.12   Khulna     0.12   Natore     0.20   Srimangal   0.36
    Comilla    0.20   Kishoreganj0.36   Naogaon    0.20   Sunamganj   0.36
    Cox's Bazar0.28   Kurigram   0.36   Netrakona  0.36   Sylhet      0.36
    Dhaka      0.20   Kushtia    0.20   Nilphamari 0.12   Tangail     0.28
    Dinajpur   0.20   Lakshmipur 0.20   Noakhali   0.20   Thakurgaon  0.20
    Faridpur   0.20   Lalmanirhat0.28   Pabna      0.20
    Feni       0.20   Madaripur  0.20   Panchagarh 0.20

    Args:
    - town (str): Town for which seismic zone coefficient is required.

    Returns:
    - zone_coefficient (float): Seismic zone coefficient (Z).
    """
    # Table 6.2.15: Seismic Zone Coefficient Z for Some Important Towns of Bangladesh
    town_coefficients = {
        "Bagerhat": 0.12, "Gaibandha": 0.28, "Magura": 0.12, "Patuakhali": 0.12,
        "Bandarban": 0.28, "Gazipur": 0.20, "Manikganj": 0.20, "Pirojpur": 0.12,
        "Barguna": 0.12, "Gopalganj": 0.12, "Maulvibazar": 0.36, "Rajbari": 0.20,
        "Barisal": 0.12, "Habiganj": 0.36, "Meherpur": 0.12, "Rajshahi": 0.12,
        "Bhola": 0.12, "Jaipurhat": 0.20, "Mongla": 0.12, "Rangamati": 0.28,
        "Bogra": 0.28, "Jamalpur": 0.36, "Munshiganj": 0.20, "Rangpur": 0.28,
        "Brahmanbaria": 0.28, "Jessore": 0.12, "Mymensingh": 0.36, "Satkhira": 0.12,
        "Chandpur": 0.20, "Jhalokati": 0.12, "Narail": 0.12, "Shariatpur": 0.20,
        "Chapainababganj": 0.12, "Jhenaidah": 0.12, "Narayanganj": 0.20, "Sherpur": 0.36,
        "Chittagong": 0.28, "Khagrachari": 0.28, "Narsingdi": 0.28, "Sirajganj": 0.28,
        "Chuadanga": 0.12, "Khulna": 0.12, "Natore": 0.20, "Srimangal": 0.36,
        "Comilla": 0.20, "Kishoreganj": 0.36, "Naogaon": 0.20, "Sunamganj": 0.36,
        "Cox's Bazar": 0.28, "Kurigram": 0.36, "Netrakona": 0.36, "Sylhet": 0.36,
        "Dhaka": 0.20, "Kushtia": 0.20, "Nilphamari": 0.12, "Tangail": 0.28,
        "Dinajpur": 0.20, "Lakshmipur": 0.20, "Noakhali": 0.20, "Thakurgaon": 0.20,
        "Faridpur": 0.20, "Lalmanirhat": 0.28, "Pabna": 0.20,
        "Feni": 0.20, "Madaripur": 0.20, "Panchagarh": 0.20
    }

    # Find the seismic zone coefficient for the given town
    if town in town_coefficients:
        return town_coefficients[town]
    else:
        return None  # Town not found in the table
# Table 6.2.17: Importance Factors for Buildings and Structures for Earthquake Design
def get_importance_factor(occupancy_category):
    """
    Get the importance factor for a given occupancy category.

    Table 6.2.17: Importance Factors for Buildings and Structures for Earthquake Design
    Occupancy Category   Importance Factor (I)
    I, II                1.00
    III                  1.25
    IV                   1.50

    Args:
    - occupancy_category (str): Occupancy category for which the importance factor is required.

    Returns:
    - importance_factor (float): Importance factor (I).
    """
    # Table 6.2.17: Importance Factors for Buildings and Structures for Earthquake Design
    importance_factors = {
        "I": 1.00,
        "II": 1.00,
        "III": 1.25,
        "IV": 1.50
    }

    # Find the importance factor for the given occupancy category
    for category, factor in importance_factors.items():
        if category == occupancy_category:
            return factor

    return None  # Occupancy category not found in the table
# Table 6.2.18: Seismic Design Category of Buildings
def get_seismic_design_category(site_class, occupancy_category, seismic_zone):
    """
    Get the seismic design category for a given site class, occupancy category, and seismic zone.

    Table 6.2.18: Seismic Design Category of Buildings
    Site Class      Occupancy Category I, II and III   Occupancy Category IV
                    Zone 1  Zone 2  Zone 3  Zone 4    Zone 1  Zone 2  Zone 3  Zone 4
    SA              B       C       C       D         C       D       D       D
    SB              B       C       D       D         C       D       D       D
    SC              B       C       D       D         C       D       D       D
    SD              C       D       D       D         D       D       D       D
    SE, S1, S2      D       D       D       D         D       D       D       D

    Args:
    - site_class (str): Site class (e.g., SA, SB, SC, SD, SE, S1, S2).
    - occupancy_category (str): Occupancy category (I, II, III, IV).
    - seismic_zone (int): Seismic zone (1, 2, 3, 4).

    Returns:
    - seismic_design_category (str): Seismic design category (A, B, C, D).
    """
    # Table 6.2.18: Seismic Design Category of Buildings
    seismic_design_categories = {
        "SA": {"I": {"1": "B", "2": "C", "3": "C", "4": "D"}, "IV": {"1": "C", "2": "D", "3": "D", "4": "D"}},
        "SB": {"I": {"1": "B", "2": "C", "3": "D", "4": "D"}, "IV": {"1": "C", "2": "D", "3": "D", "4": "D"}},
        "SC": {"I": {"1": "B", "2": "C", "3": "D", "4": "D"}, "IV": {"1": "C", "2": "D", "3": "D", "4": "D"}},
        "SD": {"I": {"1": "C", "2": "D", "3": "D", "4": "D"}, "IV": {"1": "D", "2": "D", "3": "D", "4": "D"}},
        "SE": {"I": {"1": "D", "2": "D", "3": "D", "4": "D"}, "IV": {"1": "D", "2": "D", "3": "D", "4": "D"}},
        "S1": {"I": {"1": "D", "2": "D", "3": "D", "4": "D"}, "IV": {"1": "D", "2": "D", "3": "D", "4": "D"}},
        "S2": {"I": {"1": "D", "2": "D", "3": "D", "4": "D"}, "IV": {"1": "D", "2": "D", "3": "D", "4": "D"}}
    }

    # Find the seismic design category for the given parameters
    try:
        return seismic_design_categories[site_class][occupancy_category][str(seismic_zone)]
    except KeyError:
        return None  # Invalid input combination

# Table 6.1.1: Occupancy Category of Buildings and other Structures for Flood, Surge, Wind and Earthquake Loads.

def get_occupancy_category(nature_of_occupancy):
    """
    Get the occupancy category based on the nature of occupancy.

    Table 6.1.1: Occupancy Category of Buildings and other Structures for Flood, Surge, Wind and Earthquake Loads.

    Args:
    - nature_of_occupancy (str): Nature of occupancy for which the occupancy category is required.

    Returns:
    - occupancy_category (str): Occupancy category (I, II, III, IV).
    """
    if "Agricultural facilities" in nature_of_occupancy \
            or "Certain temporary facilities" in nature_of_occupancy \
            or "Minor storage facilities" in nature_of_occupancy:
        return "I"
    elif "More than 300 people congregate in one area" in nature_of_occupancy \
            or "Day care facilities with a capacity greater than 150" in nature_of_occupancy \
            or "Elementary school or secondary school facilities with a capacity greater than 250" in nature_of_occupancy \
            or "Capacity greater than 500 for colleges or adult education facilities" in nature_of_occupancy \
            or "Healthcare facilities with a capacity of 50 or more resident patients" in nature_of_occupancy \
            or "Jails and detention facilities" in nature_of_occupancy:
        return "II"
    elif "Power generating stations" in nature_of_occupancy \
            or "Water treatment facilities" in nature_of_occupancy \
            or "Sewage treatment facilities" in nature_of_occupancy \
            or "Telecommunication centers" in nature_of_occupancy \
            or "Manufacture, process, handle, store, use, or dispose of hazardous substances" in nature_of_occupancy:
        return "III"
    elif "Hospitals and other healthcare facilities having surgery or emergency treatment facilities" in nature_of_occupancy \
            or "Fire, rescue, ambulance, and police stations and emergency vehicle garages" in nature_of_occupancy \
            or "Designated earthquake, hurricane, or other emergency shelters" in nature_of_occupancy \
            or "Designated emergency preparedness, communication, and operation centers" in nature_of_occupancy \
            or "Power generating stations and other public utility facilities required in an emergency" in nature_of_occupancy \
            or "Ancillary structures required for operation of Occupancy Category IV structures during an emergency" in nature_of_occupancy \
            or "Aviation control towers, air traffic control centers, and emergency aircraft hangars" in nature_of_occupancy \
            or "Community water storage facilities and pump structures required to maintain water pressure for fire suppression" in nature_of_occupancy \
            or "Buildings and other structures having critical national defense functions" in nature_of_occupancy \
            or "Facilities containing highly toxic substances" in nature_of_occupancy:
        return "IV"
    else:
        return None  # Nature of occupancy not found in the table

# 2.5.7.2 Building period
def calculate_building_period(height, m, x, structural_type="concrete"):
    """
    Calculate the fundamental period of a building based on the given guidelines.

    Args:
    - height (float): Height of building in meters from foundation or top of rigid basement.
    - m (float): Parameter obtained from Table 6.2.20.
    - x (int): Number of shear walls in the building.
    - structural_type (str): Type of structure (default is "concrete").

    Returns:
    - T (float): Fundamental period of the building in seconds.
    """
    if structural_type == "concrete":
        T = height / (m * 0.06) ** 0.5
    elif structural_type == "masonry":
        T = 0.0062 * (sum( [(0.83 * (Di ** 2) / hi) for Di, hi in zip( D, h )] ) + (100 / AB))
    else:
        raise ValueError( "Invalid structural type. Supported types are 'concrete' and 'masonry'." )

    return T
# Table 6.2.13: Site Classification Based on Soil Properties
def classify_spt_value(spt_value):
    if spt_value > 50:
        return "SB"
    elif 15 <= spt_value <= 50:
        return "SC"
    elif spt_value < 15:
        return "SD"
    elif 10 <= spt_value <= 20:
        return "S1"
    else:
        return "Unknown"
# Table 6.2.13: Site Classification Based on Soil Properties
def get_site_classification(shear_wave_velocity, spt_value, undrained_shear_strength):
    """
    Determine the site classification based on shear wave velocity (Vs), SPT value, and undrained shear strength.

    Table 6.2.13: Site Classification Based on Soil Properties

    Args:
    - shear_wave_velocity (float): Shear wave velocity (Vs) in m/s.
    - spt_value (int): SPT value (blows/30cm).
    - undrained_shear_strength (int): Undrained shear strength in kPa.

    Returns:
    - site_class (str): Site classification (SA, SB, SC, SD, SE, S1, S2).
    """
    if shear_wave_velocity > 800:
        return "SA"
    elif 360 <= shear_wave_velocity <= 800 or spt_value > 50 or undrained_shear_strength > 250:
        return "SB"
    elif 180 <= shear_wave_velocity < 360 or 15 <= spt_value <= 50 or 70 <= undrained_shear_strength <= 250:
        return "SC"
    elif shear_wave_velocity < 180 or spt_value < 15 or undrained_shear_strength < 70:
        return "SD"
    elif shear_wave_velocity < 100 or (spt_value >= 10 and spt_value <= 20):
        return "S1"
    else:
        return "S2"

# Table 6.2.18: Seismic Design Category of Buildings
# def get_seismic_design_category(site_class, occupancy_category, seismic_zone):
#     categories = {
#         "SA": {"1": "B", "2": "C", "3": "C", "4": "D"},
#         "SB": {"1": "B", "2": "C", "3": "D", "4": "D"},
#         "SC": {"1": "B", "2": "C", "3": "D", "4": "D"},
#         "SD": {"1": "C", "2": "D", "3": "D", "4": "D"},
#         "SE": {"1": "D", "2": "D", "3": "D", "4": "D"},
#         "S1": {"1": "D", "2": "D", "3": "D", "4": "D"},
#         "S2": {"1": "D", "2": "D", "3": "D", "4": "D"}
#     }
#
#     if occupancy_category in ['I', 'II', 'III']:
#         if seismic_zone == 1:
#             return categories[site_class]["1"]
#         elif seismic_zone == 2:
#             return categories[site_class]["2"]
#         elif seismic_zone == 3 or seismic_zone == 4:
#             return categories[site_class]["3"]  # For seismic zones 3 and 4, use the same value as for zone 3 in the table
#         else:
#             return "Unknown Seismic Zone"
#     elif occupancy_category == 'IV':
#         if seismic_zone == 1:
#             return categories[site_class]["C"]
#         elif seismic_zone == 2:
#             return categories[site_class]["D"]
#         elif seismic_zone == 3 or seismic_zone == 4:
#             return categories[site_class]["D"]
#         else:
#             return "Unknown Seismic Zone"
#     else:
#         return "Unknown Occupancy Category"

# Table 6.2.16: Site Dependent Soil Factor and Other Parameters Defining Elastic Response Spectrum
def get_soil_parameters(soil_type):
    """
    Determine the site-dependent soil factor and other parameters defining the elastic response spectrum.

    Table 6.2.16: Site Dependent Soil Factor and Other Parameters Defining Elastic Response Spectrum

    Args:
    - soil_type (str): Soil type (SA, SB, SC, SD, SE).

    Returns:
    - soil_factor (float): Site-dependent soil factor (S).
    - tb (float): TB parameter.
    - tc (float): TC parameter.
    - td (float): TD parameter.
    """
    soil_parameters = {
        "SA": (1.0, 0.15, 0.40, 2.0),
        "SB": (1.2, 0.15, 0.50, 2.0),
        "SC": (1.15, 0.20, 0.60, 2.0),
        "SD": (1.35, 0.20, 0.80, 2.0),
        "SE": (1.4, 0.15, 0.50, 2.0)
    }

    try:
        return soil_parameters[soil_type]
    except KeyError:
        return None

# Calculate Cs:
def calculate_normalized_acceleration_spectrum(S, TB, TC, TD):
    """
    Calculate the normalized acceleration response spectrum based on structure period and soil type.

    Args:
    - T (float): Structure (building) period in seconds.
    - soil_type (str): Soil type (SA, SB, SC, SD, SE).

    Returns:
    - Cs (float): Normalized acceleration response spectrum.
    η 5 Damping correction factor as a function of damping with a
    reference value of η51 for 5% viscous damping. It is given by the
    following expression:
      10 /(5   )  0.55
    """
    # Constants
    # η = Damping correction factor as a function of damping with a  reference
    #     value of η51 for 5 % viscous damping.It is given by the following expression:
    viscous_damping = 0.05
    Cs = 0
    eta = 0
    # Calculate eta
    eta = np.sqrt( 10 / (5 + viscous_damping) )
    # Ensure eta is not less than 0.55
    eta = max( eta, 0.55 )

    # Calculate Cs based on T
    if 0 <= T < TB:
        Cs = S * (1 + (T / TB) * (2.5 * eta - 1))
    elif TB <= T <= TC:
        Cs = 2.5 * S * eta
    elif TC <= T <= TD:
        Cs = 2.5 * S * eta * (TC / T)
    elif TD <= T <= 4:
        Cs = 2.5 * S * eta * (TC * TD / T ** 2)
    else:
        Cs = None  # Handle invalid range of T

    return Cs
# Table 6.2.19: Response Reduction Factor, Deflection Amplification Factor and Height Limitations for Different Structural Systems
def get_system_info(system_type, system_number):
    '''
    Table 6.2.19: Response Reduction Factor, Deflection Amplification Factor and Height
    Limitations for Different Structural Systems
    Notes:
    1. Seismic design category, NL 5 No height restriction, NP 5 Not permitted.
    Number represents maximum allowable height (m).
    2.Dual Systems include buildings which consist of both moment resisting frame
    and shear walls (or braced frame) where both systems resist the total design
    forces in proportion to their lateral stiffness.
    3. See Sec. 10.20 of Chapter 10 of this Part for additional values of R and " and
    height limits for some other types of steel structures not covered in this Table.
    4. Where data specific to a structure type is not available in this Table, reference
    may be made to Table 12.2-1 of ASCE 7-05.

    R = Response Reduction Factor,
    omega = System Overstrength Factor, omega
    Cd = Deflection Amplification Factor,
    B = Seismic Design Category B
    C = Seismic Design Category C
    D = Seismic Design Category D
    '''
    system_info = {
        "A. BEARING WALL SYSTEMS": {
            1: {"R": 5, "omega": 2.5, "Cd": 5, "B": "NL", "C": "NL", "D": 50},
            2: {"R": 4, "omega": 2.5, "Cd": 4, "B": "NL", "C": "NL", "D": "NP"},
            3: {"R": 2, "omega": 2.5, "Cd": 1.75, "B": "NL", "C": 50, "D": "NP"},
            4: {"R": 1.5, "omega": 2.5, "Cd": 1.25, "B": 18, "C": "NP", "D": "NP"}
        },
        "B.BUILDING FRAME SYSTEMS (with bracing or shear wall)": {
            1: {"R": 8, "omega": 2, "Cd": 4, "B": "NL", "C": "NL", "D": 50},
            2: {"R": 7, "omega": 2, "Cd": 4, "B": "NL", "C": "NL", "D": 50},
            3: {"R": 6, "omega": 2, "Cd": 5, "B": "NL", "C": "NL", "D": 50},
            4: {"R": 3.25, "omega": 2, "Cd": 3.25, "B": "NL", "C": "NL", "D": 11},
            5: {"R": 6, "omega": 2.5, "Cd": 5, "B": "NL", "C": "NL", "D": 50},
            6: {"R": 5, "omega": 2.5, "Cd": 4.25, "B": "NL", "C": "NL", "D": "NP"},
            7: {"R": 2, "omega": 2.5, "Cd": 2, "B": "NL", "C": 50, "D": "NP"},
            8: {"R": 1.5, "omega": 2.5, "Cd": 1.25, "B": 18, "C": "NP", "D": "NP"}
        },
        "C.MOMENT RESISTING FRAME SYSTEMS (no shear wall)": {
            1: {"R": 8, "omega": 3, "Cd": 5.5, "B": "NL", "C": "NL", "D": "NL"},
            2: {"R": 4.5, "omega": 3, "Cd": 4, "B": "NL", "C": "NL", "D": 35},
            3: {"R": 3.5, "omega": 3, "Cd": 3, "B": "NL", "C": "NL", "D": "NP"},
            4: {"R": 8, "omega": 3, "Cd": 5.5, "B": "NL", "C": "NL", "D": "NL"},
            5: {"R": 5, "omega": 3, "Cd": 4.5, "B": "NL", "C": "NL", "D": "NP"},
            6: {"R": 3, "omega": 3, "Cd": 2.5, "B": "NL", "C": "NP", "D": "NP"}
        },
        "D. DUAL SYSTEMS: SPECIAL MOMENT FRAMES CAPABLE OF RESISTING AT LEAST 25% OF PRESCRIBED SEISMIC FORCES (with bracing or shear wall)": {
            1: {"R": 8, "omega": 2.5, "Cd": 4, "B": "NL", "C": "NL", "D": "NL"},
            2: {"R": 7, "omega": 2.5, "Cd": 5.5, "B": "NL", "C": "NL", "D": "NL"},
            3: {"R": 7, "omega": 2.5, "Cd": 5.5, "B": "NL", "C": "NL", "D": "NL"},
            4: {"R": 6, "omega": 2.5, "Cd": 5, "B": "NL", "C": "NL", "D": "NP"}
        },
        "E. DUAL SYSTEMS: INTERMEDIATE MOMENT FRAMES CAPABLE OF RESISTING AT LEAST 25% OF PRESCRIBED SEISMIC FORCES (with bracing or shear wall) ": {
            1: {"R": 6, "omega": 2.5, "Cd": 5, "B": "NL", "C": "NL", "D": 11},
            2: {"R": 6.5, "omega": 2.5, "Cd": 5, "B": "NL", "C": "NL", "D": 50},
            3: {"R": 3, "omega": 3, "Cd": 3, "B": "NL", "C": 50, "D": "NP"},
            4: {"R": 5.5, "omega": 2.5, "Cd": 4.5, "B": "NL", "C": "NL", "D": "NP"}
        },
        "F. DUAL SHEAR WALL FRAME SYSTEM: ORDINARY REINFORCED CONCRETE MOMENT FRAMES AND ORDINARY REINFORCED CONCRETE SHEAR WALLS": {
            1: {"R": 4.5, "omega": 2.5, "Cd": 4, "B": "NL", "C": "NP", "D": "NP"}
        },
        "G. STEEL SYSTEMS NOT SPECIFICALLY DETAILED FOR SEISMIC RESISTANCE": {
            1: {"R": 3, "omega": 3, "Cd": 3, "B": "NL", "C": "NL", "D": "NP"}
        }
    }

    return system_info.get( system_type, {} ).get( system_number, None )












# Example
d = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
N = [2, 2, 2, 2, 12, 9, 18, 7, 14, 10, 20, 6, 16, 17, 13, 23, 10, 29]

d_float = [float(num) for num in d]
N_float = [float(num) for num in N]

print("d_float:", d_float)
print("N_float:", N_float)

corrected_SPT_value = calculate_N(d_float, N_float)
print("corrected_SPT_value:", corrected_SPT_value)

FloorHeight =3.048
BuildingHeight=30.48
'''"
Concrete moment-resisting frames": 
"Steel moment-resisting frames": 
"Eccentrically braced steel frame": 
"All other structural systems":
'''
type = "Concrete moment-resisting frames"
# MomentResistingFrameSystem  # (From Table 6.2.20)


ct, m = calculate_base_shear_and_overturning_moment(type)

# Output results
print("ct:", round(ct, 2), "kN")
print("m:", round(m, 2), "kNm")



# Table 6.2.14: and Table 6.2.15::Description of Seismic Zones
# Test the function
location = "Madaripur"  # Example location
# zone = get_seismic_zone_coefficient(location)
zone_coefficient = get_seismic_zone_coefficient_town(location)

# Output result
if zone_coefficient is not None:
    print("Seismic Zone ", location + ":", zone_coefficient)
    print("Seismic Zone Coefficient for", location + ":", zone_coefficient)
else:
    print("Seismic zone information not available for", location)

width = 19
length = 25

# Table 6.1.1: Occupancy Category of Buildings and other Structures for Flood, Surge, Wind and Earthquake Loads.

# Call the function with a sample nature_of_occupancy
nature_of_occupancy = "Hospitals and other healthcare facilities having surgery or emergency treatment facilities"

occupancy_category = get_occupancy_category(nature_of_occupancy)

if occupancy_category is not None:
    print(f"The occupancy category for the given nature of occupancy is: {occupancy_category}")
else:
    print("Nature of occupancy not found in the table.")

# Table 6.2.18: Seismic Design Category of Buildings
# Call the function with sample inputs
site_class = get_site_classification(10, corrected_SPT_value, 10)
print(f"site_class={site_class}")
seismic_zone = 2

seismic_design_category = get_seismic_design_category(site_class, occupancy_category, seismic_zone)

if seismic_design_category is not None:
    print(f"The seismic design category for site class {site_class}, occupancy category {occupancy_category}, and seismic zone {seismic_zone} is: {seismic_design_category}")
else:
    print("Invalid input combination.")

# Table 6.2.17: Importance Factors for Buildings and Structures for Earthquake Design
# Call the function with a sample occupancy category
occupancy_category = "III"

importance_factor = get_importance_factor(occupancy_category)

if importance_factor is not None:
    print(f"The importance factor for occupancy category {occupancy_category} is: {importance_factor}")
else:
    print("Occupancy category not found in the table.")
# 2.5.7.2 Building period
# Call the function with sample parameters
height = 50  # Height of building in meters
m = 2.5  # Value of m obtained from Table 6.2.20
x = 4  # Number of shear walls
structural_type = "concrete"

T = calculate_building_period(height, m, x, structural_type)
print(f"The fundamental period of the building is approximately {T:.2f} seconds.")
# Table 6.2.13: Site Classification Based on Soil Properties
spt_value = 30
print(classify_spt_value(spt_value))
# Table 6.2.13: Site Classification Based on Soil Properties
# Call the function with sample parameters
shear_wave_velocity = 00  # m/s
spt_value = 30  # blows/30cm
undrained_shear_strength = 00  # kPa

site_class = get_site_classification(shear_wave_velocity, spt_value, undrained_shear_strength)
print(f"The site classification based on the given soil properties is: {site_class}")


# Table 6.2.18: Seismic Design Category of Buildings
site_class = "SB"
occupancy_category = "II"
seismic_zone = 2

print("Inputs:")
print("Site Class:", site_class)
print("Occupancy Category:", occupancy_category)
print("Seismic Zone:", seismic_zone)

seismic_design_category = get_seismic_design_category(site_class, occupancy_category, seismic_zone)

if seismic_design_category != "Unknown":
    print(f"The seismic design category is: {seismic_design_category}")
else:
    print("Invalid input combination.")


# Table 6.2.16: Site Dependent Soil Factor and Other Parameters Defining Elastic Response Spectrum
# Call the function with sample inputs
soil_type = "SB"

soil_factor, TB, TC, TD = get_soil_parameters(soil_type)

if soil_factor is not None:
    print(f"For soil type {soil_type}:")
    print(f"Site-dependent soil factor (S): {soil_factor}")
    print(f"TB parameter: {TB}")
    print(f"TC parameter: {TC}")
    print(f"TD parameter: {TD}")
else:
    print("Invalid soil type.")


# Calculate Cs:
# Example usage:
T = 1.5  # Example structure period in seconds
soil_type = "SC"  # Example soil type

Cs = calculate_normalized_acceleration_spectrum(soil_factor, TB, TC, TD)
print(f"The normalized acceleration response spectrum (Cs) for T={T} s and soil type {soil_type} is Cs: {Cs}")

# Sa = Design spectral acceleration (in units of Q) which shall not be less than 0.67 * beta * Z * I * S
# c 5 Coefficient used to calculate lower bound for C'. Recommended
# value for c is 0.11

# Calculate Sa_min = Design spectral acceleration
beta = 0.11
Z = 0.2
I = 1
S = 1.15
Sa_min = 0.67 * beta * Z * I * S



# Table 6.2.19: Response Reduction Factor, Deflection Amplification Factor and Height Limitations for Different Structural Systems
# Example usage:
system_type = "C.MOMENT RESISTING FRAME SYSTEMS (no shear wall)"
system_number = 1
system_info = get_system_info(system_type, system_number)

R = None
omega = None
Cd = None
B = None
C = None
D = None

if system_info:
    R = system_info['R']
    omega = system_info['omega']
    Cd = system_info['Cd']
    B = system_info['B']
    C = system_info['C']
    D = system_info['D']

    print(f"Response Reduction Factor, R: {R}")
    print(f"System Overstrength Factor, omega: {omega}")
    print(f"Deflection Amplification Factor, Cd: {Cd}")
    print(f"Seismic Design Category B: {B}")
    print(f"Seismic Design Category C: {C}")
    print(f"Seismic Design Category D: {D}")
else:
    print("System information not found.")

# Example calculation of Sa, replace Z, I, and Cs with appropriate values
# Z = 1.5  # Example value for Z
# I = 2.0  # Example value for I
# Cs = 1.2  # Example value for Cs

if R is not None:
    Sa = (2 / 3) * Z * I * Cs / R
    print(f"Sa: {Sa}")
else:
    print("Cannot calculate Sa because R is not provided.")

Weight_of_building_W = 5110

# V = Sa * W
def calculate_base_shear(n, wx, hx, w, h, T):
    """
    Calculate the base shear force induced at different floor levels.

    Args:
    - n (int): Number of stories.
    - wx (list of floats): List of horizontal seismic weights at each floor level.
    - hx (list of floats): List of heights from the base to each floor level.
    - w (float): Total effective seismic weight of the structure.
    - h (float): Height from the base to the top of the structure.
    - T (float): Structure period in seconds.

    Returns:
    - Fx (list of floats): Base shear force induced at each floor level.
    """
    m = 1 if T <= 0.5 else 2 if T >= 2.5 else 1 + (T - 0.5) / (2.5 - 0.5)  # Interpolation for other periods

    Fx = [0] * n  # Initialize base shear forces at each floor level

    for i in range(n):
        Fx[i] = wx[i] * hx[i] / h * w * m

    return Fx

def calculate_storey_shear(n, Fx):
    """
    Calculate the design storey shear at each storey.

    Args:
    - n (int): Number of stories.
    - Fx (list of floats): Base shear force induced at each floor level.

    Returns:
    - Vx (list of floats): Design storey shear at each storey.
    """
    Vx = [0] * n  # Initialize storey shears

    for i in range(n):
        Vx[i] = sum(Fx[i:])  # Sum of forces from this storey and all stories above

    return Vx

# Example usage:
n = 5  # Number of stories
wx = [200, 180, 160, 140, 120]  # Horizontal seismic weights at each floor level (kN)
hx = [4, 8, 12, 16, 20]  # Heights from the base to each floor level (m)
w = sum(wx)  # Total effective seismic weight of the structure (kN)
h = sum(hx)  # Height from the base to the top of the structure (m)
T = 1.8  # Structure period (s)

# Calculate base shear force induced at different floor levels
Fx = calculate_base_shear(n, wx, hx, w, h, T)
print("Base shear forces induced at each floor level:", Fx)

# Calculate design storey shear at each storey
Vx = calculate_storey_shear(n, Fx)
print("Design storey shear at each storey:", Vx)

T = .0466*(22.7)**0.9
print(f"T = .0466*(22.7)**0.9={T}")
