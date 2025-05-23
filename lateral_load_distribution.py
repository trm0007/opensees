import json

# Given data as a list of dictionaries
elements = [
    {"name": "A", "Dx": 500, "Dy": 1750, "x": 0, "y": 22.895},
    {"name": "B", "Dx": 500, "Dy": 500, "x": 6, "y": 23.52},
    {"name": "C", "Dx": 500, "Dy": 500, "x": 15, "y": 23.52},
    {"name": "D", "Dx": 500, "Dy": 500, "x": 21, "y": 23.52},
    {"name": "E", "Dx": 300, "Dy": 700, "x": 0, "y": 15.35},
    {"name": "F", "Dx": 700, "Dy": 300, "x": 6, "y": 15.35},
    {"name": "G", "Dx": 700, "Dy": 300, "x": 15, "y": 15.35},
    {"name": "H", "Dx": 300, "Dy": 700, "x": 21, "y": 15.35},
    {"name": "I", "Dx": 500, "Dy": 500, "x": 0, "y": 7.18},
    {"name": "J", "Dx": 500, "Dy": 500, "x": 6, "y": 7.18},
    {"name": "K", "Dx": 500, "Dy": 500, "x": 15, "y": 7.18},
    {"name": "L", "Dx": 500, "Dy": 500, "x": 21, "y": 7.18},
    {"name": "M", "Dx": 500, "Dy": 500, "x": 0, "y": 0},
    {"name": "N", "Dx": 500, "Dy": 500, "x": 6, "y": 0},
    {"name": "O", "Dx": 500, "Dy": 500, "x": 15, "y": 0},
    {"name": "P", "Dx": 1750, "Dy": 500, "x": 20.375, "y": 0}
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
    element["x_A"] = x * element["A"]
    element["y_A"] = y * element["A"]

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
