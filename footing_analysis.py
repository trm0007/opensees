
import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops

# Define constants
E = 2.1e11       # Elastic modulus of the foundation material (Pa)
H = 0.5          # Height (thickness) of the foundation (m)
B = 5.0          # Width of the foundation (m)
L = 10.0         # Length of the foundation (m)
K = 1e8          # Modulus of subgrade (soil spring stiffness) reaction (Pa/m)
P = -1e6         # Compression force applied at the center of the foundation (N)

# Derived properties
nu = 0.3         # Poisson's ratio for the foundation material
rho = 2500       # Density of the foundation material (kg/m^3)

t = H            # Thickness of the shell elements

# Derived properties
I = B * H**3 / 12  # Moment of inertia for the rectangular cross-section (m^4)

# Discretization
nx = 50  # Number of grid points along length
ny = 25  # Number of grid points along width
dx = L / (nx - 1)
dy = B / (ny - 1)

nodes = nx * ny  # Total number of nodes
center_x, center_y = nx // 2, ny // 2
############################
#   FINITE ELEMENT METHOD  #
############################

# OpenSees Modeling
ops.wipe()
ops.model("BasicBuilder", "-ndm", 3, "-ndf", 6)

ops.uniaxialMaterial("ENT", 1, K)            # No Tension Soil
#ops.uniaxialMaterial("Elastic", 1, K)       # Soil
# Create nodes
node_tags = {}
spring_tags = {}
tag = 1
spring_tag = 1
NZ = 200000 # soil spring base node started 
for i in range(nx):
    for j in range(ny):
        x = i * dx
        y = j * dy
        node_tags[(i, j)] = tag
        ops.node(tag, x, y, 0.0)
        # Create a spring node directly below the foundation node
        spring_node = tag + nodes + NZ
        spring_tags[(i, j)] = spring_node
        ops.node(spring_node + NZ, x, y, 0.0)       # define base soil springs nodes
        ops.fix(spring_node + NZ, 1, 1, 1, 1, 1, 1) # fix base soil springs nodes
        # Create zero-length element for vertical spring
        ops.element("zeroLength", spring_node, spring_node + NZ, tag, "-mat", 1, "-dir", 3)
        spring_tag += 1
        tag += 1

# Define material and section
ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
ops.section("PlateFiber", 1, 1, t)

# Create ShellMITC4 elements
tag = 1
element_tags = {}
for i in range(nx - 1):
    for j in range(ny - 1):
        node1 = node_tags[(i, j)]
        node2 = node_tags[(i + 1, j)]
        node3 = node_tags[(i + 1, j + 1)]
        node4 = node_tags[(i, j + 1)]
        ops.element("ShellMITC4", tag, node1, node2, node3, node4, 1)
        element_tags[(i, j)] = tag
        tag += 1

# Apply distributed load
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

center_tag = node_tags[(center_x, center_y)]
ops.load(center_tag, 0.0, 0.0, P, 0.0, 0.0, 0.0)

# Run analysis
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system("BandGeneral")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)

# Extract results
displacements = np.zeros((nx, ny))
reactions_os = np.zeros((nx, ny))
for i in range(nx):
    for j in range(ny):
        disp = ops.nodeDisp(node_tags[(i, j)], 3)  # Z-direction displacement
        displacements[i, j] = disp
        reactions_os[i, j] = -K * disp
print(reactions_os)

################
# Plot Results #
################

import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops

import numpy as np
import openseespy.opensees as ops

# Initialize arrays for nodal results
nodal_stresses = np.zeros((nx, ny, 3))  # Sxx, Syy, Sxy
nodal_moments = np.zeros((nx, ny, 3))   # Mxx, Myy, Mxy
nodal_displacements = np.zeros((nx, ny, 3))  # Ux, Uy, Uz
nodal_reactions = np.zeros((nx, ny, 3))  # Rx, Ry, Rz

# Calculate reactions before extraction
ops.reactions()

# Extract nodal results
for i in range(nx):
    for j in range(ny):
        node_tag = node_tags[(i, j)]
        
        # --- Displacements ---
        disp = ops.nodeDisp(node_tag)
        nodal_displacements[i, j, :] = disp[:3]  # Store Ux, Uy, Uz
        
        # --- Reactions ---
        spring_node = spring_tags.get((i, j))
        if spring_node is not None:
            react = ops.nodeReaction(spring_node + NZ) if (spring_node + NZ) in ops.getNodeTags() else ops.nodeReaction(spring_node)
            if isinstance(react, (list, np.ndarray)) and len(react) >= 3:
                nodal_reactions[i, j, :] = react[:3]
            else:
                print(f"Warning: No reaction for node {spring_node}")
                nodal_reactions[i, j, :] = 0.0
        
        # --- Stresses & Moments (from connected elements) ---
        stresses = []
        moments = []
        
        # Find connected elements (4-noded quad: check 4 surrounding elements)
        for di, dj in [(0,0), (1,0), (0,1), (1,1)]:
            ii, jj = i - di, j - dj  # Check left/bottom elements
            if 0 <= ii < nx-1 and 0 <= jj < ny-1:
                elem_tag = element_tags.get((ii, jj))
                if elem_tag:
                    # --- Get stresses ---
                    stress = ops.eleResponse(elem_tag, 'stress')  # Format depends on element type!
                    if stress:
                        try:
                            stress = np.array(stress).reshape(-1, 3)[:, :3]  # Reshape to (nIP, 3)
                            stresses.append(np.mean(stress, axis=0))  # Avg integration points
                        except:
                            print(f"Warning: Stress format mismatch for element {elem_tag}")
                    
                    # --- Get moments ---
                    moment = ops.eleResponse(elem_tag, 'moment')  # For shells/plates
                    if moment:
                        try:
                            moment = np.array(moment).reshape(-1, 3)[:, :3]  # Reshape to (nIP, 3)
                            moments.append(np.mean(moment, axis=0))
                        except:
                            print(f"Warning: Moment format mismatch for element {elem_tag}")
        
        # Average results from all connected elements
        if stresses:
            nodal_stresses[i, j, :] = np.mean(stresses, axis=0)
        if moments:
            nodal_moments[i, j, :] = np.mean(moments, axis=0)

# --- Verification ---
total_reaction = np.sum(nodal_reactions[:, :, 2])  # Sum Rz components
# print(f"\nTotal vertical reaction: {total_reaction:.1f} N")
# print(f"Applied load: {P:.1f} N")
# print(f"Difference: {total_reaction - P:.1f} N ({100*(total_reaction - P)/P:.2f}%)")

# print("Nodal Stresses (Sxx, Syy, Sxy):\n", nodal_stresses)
# print("\nNodal Moments (Mxx, Myy, Mxy):\n", nodal_moments)
# print("\nNodal Displacements (Ux, Uy, Uz):\n", nodal_displacements)
print("\nNodal Reactions (Rx, Ry, Rz):\n", nodal_reactions)

# - plot model
# plt.figure()
# opsv.plot_model()
# Create grid for plotting
X, Y = np.meshgrid(np.linspace(0, L, nx), np.linspace(0, B, ny))

# Create a 3x3 grid of plots
plt.figure(figsize=(18, 15))

# 1. Stress Sxx
plt.subplot(3, 3, 1)
contour = plt.contourf(X, Y, nodal_stresses[:,:,0].T, levels=20, cmap='jet')
plt.colorbar(contour, label='Stress (Pa)')
plt.title('Stress Sxx')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 2. Stress Syy
plt.subplot(3, 3, 2)
contour = plt.contourf(X, Y, nodal_stresses[:,:,1].T, levels=20, cmap='jet')
plt.colorbar(contour, label='Stress (Pa)')
plt.title('Stress Syy')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 3. Stress Sxy
plt.subplot(3, 3, 3)
contour = plt.contourf(X, Y, nodal_stresses[:,:,2].T, levels=20, cmap='jet')
plt.colorbar(contour, label='Stress (Pa)')
plt.title('Stress Sxy')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 4. Moment Mxx
plt.subplot(3, 3, 4)
contour = plt.contourf(X, Y, nodal_moments[:,:,0].T, levels=20, cmap='viridis')
plt.colorbar(contour, label='Moment (N·m/m)')
plt.title('Moment Mxx')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 5. Moment Myy
plt.subplot(3, 3, 5)
contour = plt.contourf(X, Y, nodal_moments[:,:,1].T, levels=20, cmap='viridis')
plt.colorbar(contour, label='Moment (N·m/m)')
plt.title('Moment Myy')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 6. Moment Mxy
plt.subplot(3, 3, 6)
contour = plt.contourf(X, Y, nodal_moments[:,:,2].T, levels=20, cmap='viridis')
plt.colorbar(contour, label='Moment (N·m/m)')
plt.title('Moment Mxy')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 7. Displacement Uz
plt.subplot(3, 3, 7)
contour = plt.contourf(X, Y, nodal_displacements[:,:,2].T, levels=20, cmap='coolwarm')
plt.colorbar(contour, label='Displacement (m)')
plt.title('Vertical Displacement Uz')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 8. Reaction Rz
plt.subplot(3, 3, 8)
contour = plt.contourf(X, Y, nodal_reactions[:,:,2].T, levels=20, cmap='plasma')
plt.colorbar(contour, label='Reaction (N)')
plt.title('Vertical Reaction Rz')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

# 9. Von Mises Stress (example of derived quantity)
plt.subplot(3, 3, 9)
sxx = nodal_stresses[:,:,0]
syy = nodal_stresses[:,:,1]
sxy = nodal_stresses[:,:,2]
von_mises = np.sqrt(sxx**2 - sxx*syy + syy**2 + 3*sxy**2)
contour = plt.contourf(X, Y, von_mises.T, levels=20, cmap='hot')
plt.colorbar(contour, label='Stress (Pa)')
plt.title('Von Mises Stress')
plt.xlabel('Length (m)')
plt.ylabel('Width (m)')

plt.tight_layout()
plt.show()

# Define missing parameters for concrete foundation
fc = 25e6          # Concrete compressive strength (Pa) - typical for 25MPa concrete
column_size = 0.5   # Assume square column 0.5m x 0.5m
column_area = column_size ** 2  # 0.25 m²
cover = 0.05       # Concrete cover (50mm)
effective_depth = H - cover  # 0.45m (total thickness minus cover)

# Now we can perform the punching and beam shear checks properly
def check_punching_shear_nodal(reactions, column_area, d, fc):
    """Check punching shear using nodal reactions"""
    # Assuming square column
    col_side = np.sqrt(column_area)
    # Critical perimeter (code typically uses d/2 from column face)
    perimeter = 4 * (col_side + d)  # Simplified for square column
    
    # Total reaction in column area (using Rz)
    center_i, center_j = nx//2, ny//2
    # Determine nodes within column area
    col_nodes_x = int(np.ceil(col_side / dx))
    col_nodes_y = int(np.ceil(col_side / dy))
    
    total_shear = np.sum(nodal_reactions[
        center_i-col_nodes_x//2:center_i+col_nodes_x//2+1,
        center_j-col_nodes_y//2:center_j+col_nodes_y//2+1,
        2  # Rz component
    ])
    
    # Punching shear stress (force divided by perimeter × depth)
    v_punch = total_shear / (perimeter * d)
    
    # Allowable shear (ACI 318 simplified method)
    v_allow = 0.33 * np.sqrt(fc)  # For non-prestressed concrete
    
    print(f"\nPunching Shear Check:")
    print(f"Column size: {col_side:.2f}m x {col_side:.2f}m")
    print(f"Critical perimeter: {perimeter:.2f}m")
    print(f"Effective depth: {d:.3f}m")
    print(f"Total shear force: {total_shear/1e3:.1f} kN")
    print(f"Punching shear stress: {v_punch/1e6:.2f} MPa")
    print(f"Allowable shear stress: {v_allow/1e6:.2f} MPa")
    
    safety_factor = v_allow / v_punch
    if v_punch < v_allow:
        print(f"PASS - Safety factor: {safety_factor:.2f}")
        return True
    else:
        print(f"FAIL - {safety_factor:.2f}x over capacity")
        return False

def check_beam_shear_nodal(reactions, width, d, fc):
    """Check beam shear using nodal reactions"""
    # Critical section at distance d from column face
    col_side = np.sqrt(column_area)
    critical_distance = col_side/2 + d
    
    # Find nodes at critical section
    critical_node = int(critical_distance / dx)
    center_i = nx//2
    
    # Total shear in strip (using Rz)
    total_shear = np.sum(nodal_reactions[
        center_i-critical_node:center_i+critical_node+1,
        :,
        2  # Rz component
    ])
    
    # Shear stress (force divided by width × depth)
    v_beam = total_shear / (width * d)
    
    # Allowable shear (ACI 318 simplified method)
    v_allow = 0.17 * np.sqrt(fc)  # For non-prestressed concrete
    
    print(f"\nBeam Shear Check:")
    print(f"Critical section at: {critical_distance:.2f}m from center")
    print(f"Effective depth: {d:.3f}m")
    print(f"Total shear force: {total_shear/1e3:.1f} kN")
    print(f"Beam shear stress: {v_beam/1e6:.2f} MPa")
    print(f"Allowable shear stress: {v_allow/1e6:.2f} MPa")
    
    safety_factor = v_allow / v_beam
    if v_beam < v_allow:
        print(f"PASS - Safety factor: {safety_factor:.2f}")
        return True
    else:
        print(f"FAIL - {safety_factor:.2f}x over capacity")
        return False

# Perform the checks
check_punching_shear_nodal(nodal_reactions, column_area, effective_depth, fc)
check_beam_shear_nodal(nodal_reactions, B, effective_depth, fc)
