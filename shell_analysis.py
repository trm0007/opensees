
from matplotlib.pylab import tri
import numpy as np
import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt

ops.wipe()
# Change model to 3D with 6 DOFs
ops.model('basic', '-ndm', 3, '-ndf', 6)

# Define nodes with z-coordinate = 0.0
ops.node(1, 0., 0., 0.)
ops.node(2, 0., 1., 0.)
ops.node(3, 0., 2., 0.)
ops.node(4, 0., 3., 0.)
ops.node(5, 0., 4., 0.)
ops.node(6, 1., 0., 0.)
ops.node(7, 1., 1., 0.)
ops.node(8, 1., 2., 0.)
ops.node(9, 1., 3., 0.)
ops.node(10, 1., 4., 0.)
ops.node(11, 2., 0., 0.)
ops.node(12, 2., 1., 0.)
ops.node(13, 2., 2., 0.)
ops.node(14, 2., 3., 0.)
ops.node(15, 2., 4., 0.)
ops.node(16, 3., 0., 0.)
ops.node(17, 3., 1., 0.)
ops.node(18, 3., 2., 0.)
ops.node(19, 3., 3., 0.)
ops.node(20, 3., 4., 0.)
ops.node(21, 4., 0., 0.)
ops.node(22, 4., 1., 0.)
ops.node(23, 4., 2., 0.)
ops.node(24, 4., 3., 0.)
ops.node(25, 4., 4., 0.)

# Define material and shell section
ops.nDMaterial('ElasticIsotropic', 1, 1000, 0.3)
thickness = 0.1  # Example thickness
ops.section('PlateFiber', 2, 1, thickness)

# Replace quad elements with ShellMITC4
ops.element('ShellMITC4', 1, 1, 6, 7, 2, 2)
ops.element('ShellMITC4', 2, 2, 7, 8, 3, 2)
ops.element('ShellMITC4', 3, 3, 8, 9, 4, 2)
ops.element('ShellMITC4', 4, 4, 9, 10, 5, 2)
ops.element('ShellMITC4', 5, 6, 11, 12, 7, 2)
ops.element('ShellMITC4', 6, 7, 12, 13, 8, 2)
ops.element('ShellMITC4', 7, 8, 13, 14, 9, 2)
ops.element('ShellMITC4', 8, 9, 14, 15, 10, 2)
ops.element('ShellMITC4', 9, 11, 16, 17, 12, 2)
ops.element('ShellMITC4', 10, 12, 17, 18, 13, 2)
ops.element('ShellMITC4', 11, 13, 18, 19, 14, 2)
ops.element('ShellMITC4', 12, 14, 19, 20, 15, 2)
ops.element('ShellMITC4', 13, 16, 21, 22, 17, 2)
ops.element('ShellMITC4', 14, 17, 22, 23, 18, 2)
ops.element('ShellMITC4', 15, 18, 23, 24, 19, 2)
ops.element('ShellMITC4', 16, 19, 24, 25, 20, 2)

# Boundary conditions
all_nodes = range(1, 26)
bottom_nodes = [1, 6, 11, 16, 21]

# Fix only non-bottom nodes in Z-direction (UZ)
for node in all_nodes:
    if node not in bottom_nodes:
        ops.fix(node, 0, 0, 1, 0, 0, 0)  # Free X,Y; Fix Z translation

# Fully fix bottom nodes (X,Y,Z translations)
for node in bottom_nodes:
    ops.fix(node, 1, 1, 1, 0, 0, 0)  # Fix X,Y,Z translations but free rotations

# Equal DOF commands (couple X and Y directions)
# ops.equalDOF(2, 22, 1, 2)  # UX and UY
# ops.equalDOF(3, 23, 1, 2)
# ops.equalDOF(4, 24, 1, 2)
# ops.equalDOF(5, 25, 1, 2)

# Define load (applied to Y-direction)
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(15, 0., -1., 0., 0., 0., 0.)  # Load node 15 in Y-direction

# Add these analysis configuration commands before running the analysis
ops.system("ProfileSPD")  # Matches the default warning
ops.numberer("RCM")       # Reverse Cuthill-McKee ordering
ops.constraints("Plain")  # Plain constraint handler
ops.algorithm("Newton")   # Newton-Raphson algorithm
ops.integrator("LoadControl", 1.0)  # Load-controlled integration
ops.analysis("Static")    # Static analysis type

# Run analysis
ops.analyze(1)

# Modified plotting section (remove arrow_scale parameter)
# plt.figure()
opsv.plot_model("3D")
plt.title('3D Model Visualization')
plt.axis('equal')

# Plot loads without arrow_scale
# plt.figure()
opsv.plot_load()  # Remove arrow_scale argument
# plt.title('Load Distribution')

# Plot deformed shape with adjusted scale
# plt.figure()
opsv.plot_defo()
# plt.title('Deformed Shape')
# plt.axis('equal')

# Add these at the top of your file
import numpy as np
import matplotlib.tri as mtri  # Correct triangulation import
import matplotlib.pyplot as plt

# Modified stress plotting function
def plot_shell_stress(jstr):
    plt.figure()
    stresses = []
    coords = []
    
    # Collect element data
    for ele_tag in range(1, 17):
        # Get element response (forces in N/m)
        response = ops.eleResponse(ele_tag, 'forces')
        membrane_stress = response[:3]  # Nxx, Nyy, Nxy
        
        # Calculate centroid coordinates
        nodes = ops.eleNodes(ele_tag)
        node_coords = np.array([ops.nodeCoord(n)[:2] for n in nodes])
        centroid = np.mean(node_coords, axis=0)
        
        stresses.append(membrane_stress)
        coords.append(centroid)
    
    # Convert to numpy arrays
    coords = np.array(coords)
    stresses = np.array(stresses)
    
    # Calculate stress component values
    if jstr == 'sxx':
        values = stresses[:, 0]
    elif jstr == 'syy':
        values = stresses[:, 1]
    elif jstr == 'sxy':
        values = stresses[:, 2]
    elif jstr == 'vmis':
        sxx = stresses[:, 0]
        syy = stresses[:, 1]
        sxy = stresses[:, 2]
        values = np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*sxy**2)
    
    # Create triangulation using proper module
    triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
    
    # Create contour plot
    plt.tricontourf(triangulation, values, levels=20, cmap='jet')
    plt.colorbar()
    plt.title(f'{jstr} Stress Distribution')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.axis('equal')

# Plot stresses
stress_components = ['sxx', 'syy', 'sxy', 'vmis']
for comp in stress_components:
    plot_shell_stress(comp)

def plot_shell_moment(jstr):
    plt.figure()
    moments = []
    coords = []
    
    # Collect element data
    for ele_tag in range(1, 17):
        # Get element forces response (contains moments at indices 3-5)
        response = ops.eleResponse(ele_tag, 'forces')
        element_moments = response[3:6]  # Mxx, Myy, Mxy
        
        # Calculate element centroid
        nodes = ops.eleNodes(ele_tag)
        node_coords = np.array([ops.nodeCoord(n)[:2] for n in nodes])
        centroid = np.mean(node_coords, axis=0)
        
        moments.append(element_moments)
        coords.append(centroid)
    
    # Convert to numpy arrays
    coords = np.array(coords)
    moments = np.array(moments)
    
    # Select moment component
    if jstr == 'Mxx':
        values = moments[:, 0]
    elif jstr == 'Myy':
        values = moments[:, 1]
    elif jstr == 'Mxy':
        values = moments[:, 2]
    else:
        raise ValueError(f"Invalid moment component: {jstr}")
    
    # Create triangulation
    triangulation = mtri.Triangulation(coords[:, 0], coords[:, 1])
    
    # Create contour plot
    plt.tricontourf(triangulation, values, levels=20, cmap='jet')
    plt.colorbar()
    plt.title(f'{jstr} Moment Distribution [N·m/m]')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.axis('equal')

# Plot moments
moment_components = ['Mxx', 'Myy', 'Mxy']
for comp in moment_components:
    plot_shell_moment(comp)

plt.show()
