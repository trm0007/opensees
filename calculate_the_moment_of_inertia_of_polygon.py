import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


# METHOD 2: Improved method from the second part of your code
def calculate_properties_method2(vertices):
    """Method 2: More detailed implementation with centroid calculation"""
    # Extract x and y coordinates (ignore z for 2D calculations)
    x = vertices[:, 0]
    y = vertices[:, 1]
    
    # Close the polygon by adding the first point at the end
    x_closed = np.append(x, x[0])
    y_closed = np.append(y, y[0])
    
    # Calculate area using shoelace formula
    area = 0.5 * abs(sum(x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i] 
                         for i in range(len(x_closed)-1)))
    
    # Calculate centroid
    if area > 0:
        cx = sum((x_closed[i] + x_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i]) 
                 for i in range(len(x_closed)-1)) / (6 * area)
        cy = sum((y_closed[i] + y_closed[i+1]) * (x_closed[i] * y_closed[i+1] - x_closed[i+1] * y_closed[i]) 
                 for i in range(len(y_closed)-1)) / (6 * area)
    else:
        cx = cy = 0
    
    # Calculate second moments of area (moments of inertia)
    Ix = 0
    Iy = 0
    
    for i in range(len(x_closed)-1):
        x1, y1 = x_closed[i], y_closed[i]
        x2, y2 = x_closed[i+1], y_closed[i+1]
        
        # Contribution from this edge to Ix and Iy
        cross_product = x1 * y2 - x2 * y1
        
        Ix += cross_product * (y1**2 + y1*y2 + y2**2)
        Iy += cross_product * (x1**2 + x1*x2 + x2**2)
    
    Ix = abs(Ix) / 12
    Iy = abs(Iy) / 12
    
    return area, Ix, Iy, cx, cy

# METHOD 3: Analytical verification for simple shapes
def calculate_triangle_analytical(vertices):
    """Analytical calculation for triangle verification"""
    # For a triangle with vertices (x1,y1), (x2,y2), (x3,y3)
    x1, y1 = vertices[0, 0], vertices[0, 1]
    x2, y2 = vertices[1, 0], vertices[1, 1]
    x3, y3 = vertices[2, 0], vertices[2, 1]
    
    # Area using cross product
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    
    # Centroid
    cx = (x1 + x2 + x3) / 3
    cy = (y1 + y2 + y3) / 3
    
    # Second moments about centroidal axes (approximate for verification)
    # Using the fact that for a triangle, the moments can be calculated analytically
    # but the exact formula is complex, so we'll use the numerical methods
    
    return area, cx, cy

# Test all methods
print("="*80)
print("COMPARISON OF METHODS FOR CALCULATING POLYGON PROPERTIES")
print("="*80)

# Define the same shapes for comparison
shape1 = np.array([
    [0, 0, 0],
    [2, 0.5, 0],
    [2.5, 3, 0],
    [0.5, 2.5, 0]
])

shape2 = np.array([
    [3, 1, 0],
    [4.5, 0.5, 0],
    [4, 2.5, 0]
])

area1_m2, Ix1_m2, Iy1_m2, cx1_m2, cy1_m2 = calculate_properties_method2(shape1)
area2_m2, Ix2_m2, Iy2_m2, cx2_m2, cy2_m2 = calculate_properties_method2(shape2)
area2_analytical, cx2_analytical, cy2_analytical = calculate_triangle_analytical(shape2)

print(f"\nSHAPE 1 - QUADRILATERAL:")


print(f"Method 2 Results:")
print(f"  Area: {area1_m2:.6f}")
print(f"  Ix: {Ix1_m2:.6f}")
print(f"  Iy: {Iy1_m2:.6f}")
print(f"  Centroid: ({cx1_m2:.6f}, {cy1_m2:.6f})")



# Test Shape 2 (Triangle)



print(f"\nSHAPE 2 - TRIANGLE:")

print(f"  Area: {area2_m2:.6f}")
print(f"  Ix: {Ix2_m2:.6f}")
print(f"  Iy: {Iy2_m2:.6f}")
print(f"  Centroid: ({cx2_m2:.6f}, {cy2_m2:.6f})")

print(f"Analytical verification:")
print(f"  Area: {area2_analytical:.6f}")
print(f"  Centroid: ({cx2_analytical:.6f}, {cy2_analytical:.6f})")



# Function to create mesh elements (each with 4 nodes) in a rectangular grid using x, y, z (y=0)



from shapely.geometry import Polygon


def calculate_moments_of_inertia(polygon):
    # Get centroid
    centroid = polygon.centroid
    cx, cy = centroid.x, centroid.y
    
    # Translate coordinates to centroidal system
    exterior_coords = list(polygon.exterior.coords[:-1])  # Remove duplicate
    translated_coords = [(x - cx, y - cy) for x, y in exterior_coords]
    
    n = len(translated_coords)
    Ix = Iy = 0.0
    
    for i in range(n):
        x_i, y_i = translated_coords[i]
        x_j, y_j = translated_coords[(i + 1) % n]
        
        cross_term = x_i * y_j - x_j * y_i
        Ix += (y_i**2 + y_i * y_j + y_j**2) * cross_term
        Iy += (x_i**2 + x_i * x_j + x_j**2) * cross_term
    
    return abs(Ix / 12.0), abs(Iy / 12.0)


polygon = Polygon([(0, 0), (4, 0), (4, 3), (0, 3)])
Ix, Iy = calculate_moments_of_inertia(polygon)
print(f"Numerical (Polygon) Ix: {Ix}, Iy: {Iy}")

# Expected analytical solution
width = 4  # x-direction
height = 3  # y-direction
Ix_expected = width * height**3 / 12  # 9.0
Iy_expected = height * width**3 / 12  # 16.0
print(f"Analytical (Rectangle) Ix: {Ix_expected}, Iy: {Iy_expected}")
