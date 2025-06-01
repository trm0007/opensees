import math

def structural_analysis(Lx, Ly,E, poisson_ratio, fc, Vx, Vy, components, elements):
    """Perform all structural calculations in one function."""
    
    # 1. Calculate center of mass
    # print("\n1. CENTER OF MASS CALCULATION:")
    total_wt = sum(comp["wt"] for comp in components)
    sum_wt_xi = sum(comp["wt"] * comp["xi"] for comp in components)
    sum_wt_yi = sum(comp["wt"] * comp["yi"] for comp in components)
    Xm = sum_wt_xi / total_wt
    Ym = sum_wt_yi / total_wt
    COM =(Xm, Ym)
    total_Kx = 0
    total_Ky = 0
    sum_Kx_x = 0
    sum_Ky_y = 0
    element_details = []
    
    for element in elements:
        Dx = element["Dx"]
        Dy = element["Dy"]
        x = element["x"]
        y = element["y"]
        height = element["height"]
        element_type = element.get("type")
        
        # Calculate stiffness in both directions
        Ix = (Dy * Dx**3) / 12
        Iy = (Dx * Dy**3) / 12
        
        
        if element_type == "column":
            Kx = (12 * E * Ix) / (height**3)
            Ky = (12 * E * Iy) / (height**3)
        else:
            denom_x = 1 + 0.6 * (1 + poisson_ratio) * (Dx**2 / height**2)
            denom_y = 1 + 0.6 * (1 + poisson_ratio) * (Dy**2 / height**2)
            Kx = (3 * E * Ix) / (height**3 * denom_x)
            Ky = (3 * E * Iy) / (height**3 * denom_y)
        
        element_details.append({
            "name": element.get("name", "Unknown"),
            "Kx": Kx,
            "Ky": Ky,
            "x": x,
            "y": y
        })
        
        total_Kx += Kx
        total_Ky += Ky
        sum_Kx_x += Kx * x
        sum_Ky_y += Ky * y
    
    Xr = sum_Kx_x / total_Kx if total_Kx > 0 else 0
    Yr = sum_Ky_y / total_Ky if total_Ky > 0 else 0
    COR = (Xr, Yr)


    # 3. Calculate eccentricity
    # print("\n3. ECCENTRICITY CALCULATION:")
    ex = Xm - Xr
    ey = Ym - Yr
    min_ex = 0.05 * Lx
    min_ey = 0.05 * Ly
    ex_total = ex + (-min_ex if ex < 0 else min_ex)
    ey_total = ey + (-min_ey if ey < 0 else min_ey)
    eccentricity = (ex, ey)
    eccentricity_total = (ex_total, ey_total)

    # 4. Calculate torsional moment
    T = Vy * ex_total + Vx * ey_total

    # 5. Calculate torsional rigidity
    # print("\n4. TORSIONAL RIGIDITY CALCULATION:")
    Jr_total = 0
    Jr_X = 0
    Jr_Y = 0
    

    for element in element_details:
        dx = element["x"] - Xr
        dy = element["y"] - Yr
        Kx = element["Kx"]
        Ky = element["Ky"]

        Jr_contribution_x = Kx * dy**2
        Jr_contribution_y = Ky * dx**2
        
        Jr_X += Jr_contribution_x
        Jr_Y += Jr_contribution_y

        # print(f"{element['name']:<8} {Kx/1e6:<12.2f} {Ky/1e6:<12.2f} {dx:<8.2f} {dy:<8.2f} {Jr_contribution_x/1e6:<12.2f} {Jr_contribution_y/1e6:<12.2f}")

    Jr_total = Jr_X + Jr_Y
    torsional_rigidity = (Jr_X, Jr_Y, Jr_total)

    # 6. Calculate floor shear distribution
    print("\n5. FLOOR SHEAR DISTRIBUTION:")
    results = []
    total_Fx = 0
    total_Fy = 0
    
    for element in element_details:
        dx = element['x'] - Xr
        dy = element['y'] - Yr
        Kx = element['Kx']
        Ky = element['Ky']
        
        Fy = (Vy * Ky / total_Ky) + (T * Kx * dx / Jr_total)
        Fx = (Vx * Kx / total_Kx) - (T * Ky * dy / Jr_total)
        
        results.append({
            'name': element['name'],
            'Fx': Fx,
            'Fy': Fy,
            'dx': dx,
            'dy': dy
        })
        total_Fx += Fx
        total_Fy += Fy
        
    return {
        "center_of_mass": COM,
        "center_of_rigidity": COR,
        "total_stiffness": (total_Kx, total_Ky),
        "element_details": element_details,
        "eccentricity": eccentricity,
        "eccentricity_total": eccentricity_total,
        "torsional_moment": T,
        "torsional_rigidity": torsional_rigidity,
        "shear_distribution": results,
        "total_shear": (total_Fx, total_Fy)
    }

# Example usage with the same input data
components = [
    {"wt": 2800, "xi": 8, "yi": 6.25},
    {"wt": -84, "xi": 14, "yi": 11.75},
    {"wt": 115.2, "xi": 0.15, "yi": 6},
    {"wt": 115.2, "xi": 8, "yi": 6},
    {"wt": 115.2, "xi": 15.85, "yi": 6},
    {"wt": 115.2, "xi": 2, "yi": 0.15},
    {"wt": 115.2, "xi": 14, "yi": 0.15},
    {"wt": 115.2, "xi": 2, "yi": 12.35},
    {"wt": 115.2, "xi": 14, "yi": 10.85},
]

elements = [
    {"name": "W2-1", "type": "wall", "Dx": 0.3, "Dy": 4.0, "height": 4.0, "x": 2.0, "y": 0.15},
    {"name": "W2-2", "type": "wall", "Dx": 0.3, "Dy": 4.0, "height": 4.0, "x": 14.0, "y": 0.15},
    {"name": "W2-3", "type": "wall", "Dx": 0.3, "Dy": 4.0, "height": 4.0, "x": 2.0, "y": 12.35},
    {"name": "W2-4", "type": "wall", "Dx": 0.3, "Dy": 4.0, "height": 4.0, "x": 14.0, "y": 10.85},
    {"name": "W1-1", "type": "wall", "Dx": 4.0, "Dy": 0.3, "height": 4.0, "x": 0.15, "y": 6.0},
    {"name": "W1-2", "type": "wall", "Dx": 4.0, "Dy": 0.3, "height": 4.0, "x": 8.0, "y": 6.0},
    {"name": "W1-3", "type": "wall", "Dx": 4.0, "Dy": 0.3, "height": 4.0, "x": 15.85, "y": 6.0},
]

Lx = 16.0
Ly = 12.5
fc = 30
E = 4700 * 10**6 * math.sqrt(fc)
poisson_ratio = 0.25
fc = 30
Vx = 0.0
Vy = -34.163

results = structural_analysis(Lx, Ly, E, poisson_ratio, fc, Vx, Vy, components, elements)

# Print all results in a clean, organized format
print("\n" + "="*50)
print("STRUCTURAL ANALYSIS RESULTS".center(50))
print("="*50)

# 1. Center of Mass and Rigidity
print("\n1. CENTER OF MASS & RIGIDITY")
print("-"*50)
print(f"{'Center of Mass (COM)':<30}: X = {results['center_of_mass'][0]:.3f} m, Y = {results['center_of_mass'][1]:.3f} m")
print(f"{'Center of Rigidity (COR)':<30}: X = {results['center_of_rigidity'][0]:.3f} m, Y = {results['center_of_rigidity'][1]:.3f} m")

# 2. Stiffness Properties
print("\n2. STIFFNESS PROPERTIES")
print("-"*50)
print(f"{'Total Stiffness in X-direction':<30}: {results['total_stiffness'][0]/1e6:.3f} x10^6 kN/m")
print(f"{'Total Stiffness in Y-direction':<30}: {results['total_stiffness'][1]/1e6:.3f} x10^6 kN/m")

# 3. Eccentricity
print("\n3. ECCENTRICITY")
print("-"*50)
print(f"{'Natural Eccentricity (ex, ey)':<30}: ex = {results['eccentricity'][0]:.3f} m, ey = {results['eccentricity'][1]:.3f} m")
print(f"{'Design Eccentricity (ex_total, ey_total)':<30}: ex = {results['eccentricity_total'][0]:.3f} m, ey = {results['eccentricity_total'][1]:.3f} m")

# 4. Torsional Properties
print("\n4. TORSIONAL PROPERTIES")
print("-"*50)
print(f"{'Torsional Moment (T)':<30}: {results['torsional_moment']/1e3:.3f} kN·m")
print(f"{'Torsional Rigidity (Jr_X, Jr_Y)':<30}: {results['torsional_rigidity'][0]/1e6:.3f}, {results['torsional_rigidity'][1]/1e6:.3f} x10^6 kN·m²")
print(f"{'Total Torsional Rigidity (Jr_total)':<30}: {results['torsional_rigidity'][2]/1e6:.3f} x10^6 kN·m²")

# 5. Element Stiffness Details
print("\n5. ELEMENT STIFFNESS DETAILS")
print("-"*50)
print(f"{'Element':<10}{'Kx (kN/m)':<15}{'Ky (kN/m)':<15}{'X (m)':<10}{'Y (m)':<10}")
print("-"*50)
for element in results['element_details']:
    print(f"{element['name']:<10}{element['Kx']/1e3:<15.3f}{element['Ky']/1e3:<15.3f}{element['x']:<10.3f}{element['y']:<10.3f}")

# 6. Shear Distribution
print("\n6. SHEAR DISTRIBUTION")
print("-"*50)
print(f"{'Element':<10}{'Fx (kN)':<15}{'Fy (kN)':<15}{'dx (m)':<10}{'dy (m)':<10}")
print("-"*50)
for result in results['shear_distribution']:
    print(f"{result['name']:<10}{result['Fx']:<15.3f}{result['Fy']:<15.3f}{result['dx']:<10.3f}{result['dy']:<10.3f}")

# Print totals
print("-"*50)
print(f"{'TOTAL':<10}{results['total_shear'][0]:<15.3f}{results['total_shear'][1]:<15.3f}")
print("="*50)
