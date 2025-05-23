import json

# Given data as a list of dictionaries
elements = [
    {"name": "A", "Dx": 500, "Dy": 1750, "x": 0, "y": 22.895, "weight":10},
    {"name": "B", "Dx": 500, "Dy": 500, "x": 6, "y": 23.52, "weight":10},
    {"name": "C", "Dx": 500, "Dy": 500, "x": 15, "y": 23.52, "weight":10},
    {"name": "D", "Dx": 500, "Dy": 500, "x": 21, "y": 23.52, "weight":10},
    {"name": "E", "Dx": 300, "Dy": 700, "x": 0, "y": 15.35, "weight":10},
    {"name": "F", "Dx": 700, "Dy": 300, "x": 6, "y": 15.35, "weight":10},
    {"name": "G", "Dx": 700, "Dy": 300, "x": 15, "y": 15.35, "weight":10},
    {"name": "H", "Dx": 300, "Dy": 700, "x": 21, "y": 15.35, "weight":10},
    {"name": "I", "Dx": 500, "Dy": 500, "x": 0, "y": 7.18, "weight":10},
    {"name": "J", "Dx": 500, "Dy": 500, "x": 6, "y": 7.18, "weight":10},
    {"name": "K", "Dx": 500, "Dy": 500, "x": 15, "y": 7.18, "weight":10},
    {"name": "L", "Dx": 500, "Dy": 500, "x": 21, "y": 7.18, "weight":10},
    {"name": "M", "Dx": 500, "Dy": 500, "x": 0, "y": 0, "weight":10},
    {"name": "N", "Dx": 500, "Dy": 500, "x": 6, "y": 0, "weight":10},
    {"name": "O", "Dx": 500, "Dy": 500, "x": 15, "y": 0, "weight":10},
    {"name": "P", "Dx": 1750, "Dy": 500, "x": 20.375, "y": 0, "weight":10}
]

# Calculate kx and ky for each element
for element in elements:
    Dx = element["Dx"]  # in mm
    Dy = element["Dy"]  # in mm
    x = element["x"]*1000    # in m
    y = element["y"]*1000    # in m
    
    # Calculate moments of inertia (convert mm to m for consistent units)
    element["kx"] = (1/12) * (Dy) * ((Dx) ** 3)  # kx = (1/12) * a * b³
    element["ky"] = (1/12) * (Dx) * ((Dy) ** 3)  # ky = (1/12) * b * a³
    
    # Calculate products for center of rigidity
    element["kx_y"] = element["kx"] * y
    element["ky_x"] = element["ky"] * x

    # Calculate area
    element["A"] = Dx * Dy
    
    # Calculate weighted positions
    element["x_A"] = x * element["A"] * element["weight"]
    element["y_A"] = y * element["A"] * element["weight"]

# Sum all kx, ky, kx_y, and ky_x
sum_kx = sum(element["kx"] for element in elements)
sum_ky = sum(element["ky"] for element in elements)
sum_kx_y = sum(element["kx_y"] for element in elements)
sum_ky_x = sum(element["ky_x"] for element in elements)
sum_A = sum(element["A"] for element in elements)
sum_x_A = sum(element["x_A"] for element in elements)
sum_y_A = sum(element["y_A"] for element in elements)

# Calculate center of rigidity (CR)
CR_y = sum_kx_y / sum_ky if sum_ky != 0 else 0
CR_x = sum_ky_x / sum_kx if sum_kx != 0 else 0

# Calculate center of mass (CM)
CM_x = sum_x_A / sum_A if sum_A != 0 else 0
CM_y = sum_y_A / sum_A if sum_A != 0 else 0

enx = CR_x - CM_x
eny = CM_y - CR_y

Wx = 21.4*1000
Wy = 24.02*1000

ex = enx-0.05*Wx if enx<0 else enx+0.05*Wx
ey = eny-0.05*Wy if eny<0 else eny+0.05*Wy

Vy = -34.163 
Vx = 0.0
Ty = Vx * ey
Tx = Vy * ex
T = (Tx + Ty)/1000

# NEW PART: Calculate additional parameters for each element
for element in elements:
    x = element["x"]*1000    # in mm
    y = element["y"]*1000    # in mm
    
    # Calculate relative positions from center of rigidity
    element["xr"] = CR_x - x
    element["yr"] = y - CR_y
    
    # Calculate ky * xr² and kx * yr² (convert xr and yr to meters for consistent units)
    element["ky_xr2"] = element["ky"] * ((element["xr"]/1000) ** 2)
    element["kx_yr2"] = element["kx"] * ((element["yr"]/1000) ** 2)

# Sum the new parameters
sum_ky_xr2 = sum(element["ky_xr2"] for element in elements)
sum_kx_yr2 = sum(element["kx_yr2"] for element in elements)

# Calculate J (polar moment of inertia)
J = sum_ky_xr2 + sum_kx_yr2

# Calculate Fx and Fy for each element
for element in elements:
    # Fy = (Vy * ky)/Σky + (T * ky * xr)/J
    element["Fy"] = (Vy * element["ky"])/sum_ky + (T * element["ky"] * (element["xr"]/1000))/J
    
    # Fx = (Vx * kx)/Σkx + (T * kx * yr)/J
    element["Fx"] = (Vx * element["kx"])/sum_kx + (T * element["kx"] * (element["yr"]/1000))/J

# Prepare results
results = {
    "elements": elements,
    "sum_kx": sum_kx,
    "sum_ky": sum_ky,
    "sum_kx_y": sum_kx_y,
    "sum_ky_x": sum_ky_x,
    "center_of_rigidity": {"x": CR_x, "y": CR_y},
    "center_of_mass": {"x": CM_x, "y": CM_y},
    "enx" : enx,
    "eny" : eny,
    "ex" : ex,
    "ey" : ey,
    "Tx" : Tx,
    "Ty" : Ty,
    "T" : T,
    "sum_ky_xr2": sum_ky_xr2,
    "sum_kx_yr2": sum_kx_yr2,
    "J": J,
}

# Convert to JSON
json_data = json.dumps(results, indent=4)

# Print results
print(json_data)

# Print Fx and Fy for each element
print("\n" + "="*50)
print("Forces for each element:")
print("="*50)
for element in elements:
    print(f"Element {element['name']}: Fx = {element['Fx']:.3f}, Fy = {element['Fy']:.3f}")

# # Optionally, save to a file
# with open("center_of_rigidity_results.json", "w") as f:
#     f.write(json_data)



import math

# ========================================
# WALL RIGIDITY CALCULATIONS
# ========================================

def calculate_wall_rigidity(h, L):
    """
    Calculate wall rigidity using the given formula
    R = 1 / (0.4 * (h/d)^3 + 0.3 * (h/d))
    Where h = height, d = length (L)
    """
    h_over_d = h / L
    denominator = 0.4 * (h_over_d)**3 + 0.3 * h_over_d
    R = 1 / denominator
    return R

# ========================================
# COLUMN RIGIDITY CALCULATIONS
# ========================================

def calculate_column_rigidity(E, I, L):
    """
    Calculate column rigidity using k = 12EI/L³
    
    Parameters:
    E: Modulus of elasticity (psi or Pa)
    I: Moment of inertia (in⁴ or m⁴)
    L: Column length/height (in or m)
    
    Returns:
    k: Column rigidity
    """
    k = (12 * E * I) / (L**3)
    return k

def calculate_moment_of_inertia_rectangular(b, h):
    """Calculate moment of inertia for rectangular cross-section: I = bh³/12"""
    I = (b * h**3) / 12
    return I

def calculate_moment_of_inertia_circular(d):
    """Calculate moment of inertia for circular cross-section: I = πd⁴/64"""
    I = (math.pi * d**4) / 64
    return I

# ========================================
# MAIN CALCULATIONS
# ========================================

def main():
    print("=" * 60)
    print("COMBINED WALL AND COLUMN RIGIDITY CALCULATOR")
    print("=" * 60)
    
    # ========================================
    # WALL RIGIDITY CALCULATIONS
    # ========================================
    print("\nWALL RIGIDITY CALCULATIONS")
    print("-" * 40)
    print("Formula: R = 1 / (0.4 * (h/d)³ + 0.3 * (h/d))")
    print("Where h = height, d = length")
    print()

    # Wall data
    walls = [
        {"name": "10ft Wall", "h": 9, "L": 10},
        {"name": "20ft Wall", "h": 9, "L": 20}
    ]

    wall_results = []
    for wall in walls:
        h = wall["h"]
        L = wall["L"]
        
        # Calculate h/d ratio
        h_over_d = h / L
        
        # Calculate denominator components
        cubic_term = 0.4 * (h_over_d)**3
        linear_term = 0.3 * h_over_d
        denominator = cubic_term + linear_term
        
        # Calculate rigidity
        R = calculate_wall_rigidity(h, L)
        wall_results.append({"name": wall["name"], "R": R})
        
        print(f"{wall['name']}:")
        print(f"  h = {h} ft, L = {L} ft")
        print(f"  h/d = {h}/{L} = {h_over_d:.3f}")
        print(f"  0.4 * (h/d)³ = 0.4 * ({h_over_d:.3f})³ = {cubic_term:.6f}")
        print(f"  0.3 * (h/d) = 0.3 * {h_over_d:.3f} = {linear_term:.6f}")
        print(f"  Denominator = {cubic_term:.6f} + {linear_term:.6f} = {denominator:.6f}")
        print(f"  R = 1 / {denominator:.6f} = {R:.2f}")
        print()

    # ========================================
    # COLUMN RIGIDITY CALCULATIONS
    # ========================================
    print("\nCOLUMN RIGIDITY CALCULATIONS")
    print("-" * 40)
    print("Formula: k = 12EI/L³")
    print()

    # Example 1: Concrete column
    print("Example 1: Concrete Column")
    E_concrete = 3000000  # psi (typical concrete E)
    b = 12  # width in inches
    h = 12  # height in inches
    L = 120  # column length in inches (10 ft)

    I_concrete = calculate_moment_of_inertia_rectangular(b, h)
    k_concrete = calculate_column_rigidity(E_concrete, I_concrete, L)

    print(f"  Material: Concrete")
    print(f"  E = {E_concrete:,} psi")
    print(f"  Cross-section: {b}\" × {h}\"")
    print(f"  I = bh³/12 = {b} × {h}³/12 = {I_concrete:.1f} in⁴")
    print(f"  L = {L} in")
    print(f"  k = 12 × {E_concrete:,} × {I_concrete:.1f} / {L}³")
    print(f"  k = {k_concrete:,.0f} lb/in")
    print()

    # Example 2: Steel column
    print("Example 2: Steel Column")
    E_steel = 29000000  # psi (steel E)
    d = 10  # diameter in inches
    L = 120  # column length in inches

    I_steel = calculate_moment_of_inertia_circular(d)
    k_steel = calculate_column_rigidity(E_steel, I_steel, L)

    print(f"  Material: Steel")
    print(f"  E = {E_steel:,} psi")
    print(f"  Cross-section: Circular, d = {d}\"")
    print(f"  I = πd⁴/64 = π × {d}⁴/64 = {I_steel:.1f} in⁴")
    print(f"  L = {L} in")
    print(f"  k = 12 × {E_steel:,} × {I_steel:.1f} / {L}³")
    print(f"  k = {k_steel:,.0f} lb/in")
    print()

    # ========================================
    # SUMMARY
    # ========================================
    print("\nSUMMARY OF RESULTS")
    print("-" * 40)
    print("Wall Rigidities:")
    for result in wall_results:
        print(f"  {result['name']}: R = {result['R']:.2f}")
    
    print(f"\nColumn Rigidities:")
    print(f"  Concrete Column (12\"×12\"): k = {k_concrete:,.0f} lb/in")
    print(f"  Steel Column (d=10\"): k = {k_steel:,.0f} lb/in")

# ========================================
# CUSTOM CALCULATION FUNCTIONS
# ========================================

def calculate_custom_wall(h, L):
    """Calculate wall rigidity with custom dimensions"""
    R = calculate_wall_rigidity(h, L)
    h_over_d = h / L
    print(f"Custom Wall: h={h}, L={L}")
    print(f"  h/d = {h_over_d:.3f}")
    print(f"  R = {R:.2f}")
    return R

def calculate_custom_column(E, cross_section_type, dimensions, L):
    """Calculate column rigidity with custom inputs"""
    if cross_section_type.lower() == "rectangular":
        b, h = dimensions
        I = calculate_moment_of_inertia_rectangular(b, h)
        print(f"Custom Column - Rectangular: {b}\" × {h}\"")
    elif cross_section_type.lower() == "circular":
        d = dimensions[0]
        I = calculate_moment_of_inertia_circular(d)
        print(f"Custom Column - Circular: d = {d}\"")
    else:
        I = dimensions[0]  # Assume I is given directly
        print(f"Custom Column - Given I = {I} in⁴")
    
    k = calculate_column_rigidity(E, I, L)
    
    print(f"  E = {E:,} psi")
    print(f"  I = {I:.1f} in⁴")
    print(f"  L = {L} in")
    print(f"  k = {k:,.0f} lb/in")
    return k

# ========================================
# USAGE EXAMPLES
# ========================================

print("\nUSAGE EXAMPLES:")
print("-" * 40)
print("# For wall rigidity:")
print("R = calculate_wall_rigidity(height, length)")
print("Example: R = calculate_wall_rigidity(9, 15)")
print(f"Result: {calculate_wall_rigidity(9, 15):.2f}")
print()

print("# For column rigidity:")
print("k = calculate_column_rigidity(E, I, L)")
print("Example: k = calculate_column_rigidity(3000000, 1728, 120)")
print(f"Result: {calculate_column_rigidity(3000000, 1728, 120):,.0f} lb/in")

# Run main calculations
if __name__ == "__main__":
    main()
