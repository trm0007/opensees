import os
import time
import numpy as np
import openseespy.opensees as ops
from collections import defaultdict

# Constants
KN = 1.0  # Assuming kN as base unit
KN_m = 1.0  # kN/m

# *******---------********-------xxxxxxxxx--------*********---------xxxxxxx
# 1. Load Case Definition - Improved with better structure
class LoadCase:
    def __init__(self, name, case_type, factor=1.0):
        self.name = name
        self.type = case_type  # 'gravity' or 'lateral'
        self.factor = factor
        self.nodal_loads = defaultdict(list)  # {node_tag: [Fx, Fy, Fz, Mx, My, Mz]}
        self.member_loads = []  # [{'eleTags': [], 'type': 'beamUniform', 'values': [Wy, Wz, Wx]}]
        self.surface_loads = []  # [{'eleTags': [], 'face': face_id, 'pressure': p}]
        self.direction = None  # For lateral loads
        self.mass_loads = defaultdict(list)  # {node_tag: [massX, massY, massZ, massRotX, massRotY, massRotZ]}
        
    def add_nodal_load(self, node_tag, loads):
        self.nodal_loads[node_tag] = loads
        
    def add_member_load(self, ele_tags, load_type, values):
        self.member_loads.append({
            'eleTags': ele_tags,
            'type': load_type,
            'values': values
        })
        
    def add_surface_load(self, ele_tags, face, pressure):
        self.surface_loads.append({
            'eleTags': ele_tags,
            'face': face,  # Face ID depends on element type
            'pressure': pressure
        })
        
    def add_mass_load(self, node_tag, mass_values):
        """Add mass to a node.
        
        Args:
            node_tag: Node identifier
            mass_values: List of mass values [mx, my, mz, mRotX, mRotY, mRotZ]
                        (typically rotational masses are zero)
        """
        self.mass_loads[node_tag] = mass_values


# 2. Applying Loads - Improved with better surface load and mass handling
def apply_loads(load_case, pattern_tag=1):
    """Apply nodal, member, surface loads, and masses for a given load case."""
    ops.pattern('Plain', pattern_tag, 1)
    
    # Apply nodal loads
    for node_tag, load_values in load_case.nodal_loads.items():
        ops.load(node_tag, *load_values)
    
    # Apply member loads
    for member_load in load_case.member_loads:
        ele_tags = member_load['eleTags']
        load_type = member_load['type']
        values = member_load['values']
        
        if 'beamUniform' in load_type:
            ops.eleLoad('-ele', *ele_tags, '-type', '-beamUniform', *values)
    
    # Apply surface loads
    for surface_load in load_case.surface_loads:
        ele_tags = surface_load['eleTags']
        pressure = surface_load['pressure']
        face = surface_load.get('face', 1)  # Default to face 1 if not specified
        
        for ele_tag in ele_tags:
            face_nodes = get_face_nodes(ele_tag, face)
            if face_nodes:
                surf_tag = 1000 + ele_tag  
                ops.element('SurfaceLoad', surf_tag, *face_nodes, pressure)
            else:
                print(f"Warning: Could not determine face nodes for element {ele_tag}")
    
    # Apply mass loads (NEW)
    for node_tag, mass_values in load_case.mass_loads.items():
        ops.mass(node_tag, *mass_values)
    
    return pattern_tag

def get_face_nodes(ele_tag, face):
    """Example implementation for getting face nodes for different element types"""
    # For brick elements (8 nodes), face node ordering varies by face
    if ops.eleType(ele_tag) == 'stdBrick':
        all_nodes = ops.eleNodes(ele_tag)
        if face == 1:   # Front face
            return [all_nodes[0], all_nodes[1], all_nodes[2], all_nodes[3]]
        elif face == 2: # Right face
            return [all_nodes[1], all_nodes[5], all_nodes[6], all_nodes[2]]
        # Add other faces as needed
    
    # For shell elements (4 nodes), same as element nodes
    elif ops.eleType(ele_tag) == 'ShellMITC4':
        return ops.eleNodes(ele_tag)
    
    return None


# 3. Running Analysis for Combinations - Improved with better error handling
def run_analysis_for_combinations(combinations_dict, load_cases_dict, analysis_type="gravity"):
    """
    Run analysis for each load combination with specified analysis type.
    
    Args:
        combinations_dict: Dictionary of load combinations
        load_cases_dict: Dictionary of load cases
        analysis_type: One of ["gravity", "modal", "pushover", "time_history"]
    """
    results = {}
    
    for comb_name, case_factors in combinations_dict.items():
        print(f"\nRunning {analysis_type} analysis for combination: {comb_name}")
        
        try:
            # Reset model
            ops.wipe()
            build_model()  # Your model building function
            
            # Apply all loads in combination (scaled by factors)
            pattern_tags = []
            for case_name, factor in case_factors:
                pattern_tag = len(pattern_tags) + 1
                scaled_case = scale_load_case(load_cases_dict[case_name], factor)
                apply_loads(scaled_case, pattern_tag)
                pattern_tags.append(pattern_tag)
            
            # Run specified analysis type
            if analysis_type == "gravity":
                success = run_gravity(comb_name=comb_name)
            elif analysis_type == "modal":
                success = run_modal(comb_name=comb_name)
            elif analysis_type == "pushover":
                # Determine direction from lateral load cases
                lateral_case = next((c for c,_ in case_factors 
                                  if load_cases_dict[c].type == 'lateral'), None)
                direction = load_cases_dict[lateral_case].direction if lateral_case else 'X'
                success = run_pushover(m_1, direction=direction, comb_name=comb_name)
            elif analysis_type == "time_history":
                lateral_case = next((c for c,_ in case_factors 
                                   if 'Seismic' in c), None)
                direction = load_cases_dict[lateral_case].direction if lateral_case else 'X'
                success = run_time_history(direction=direction, comb_name=comb_name)
            else:
                raise ValueError(f"Unknown analysis type: {analysis_type}")
            
            # Store results
            results[comb_name] = {
                'status': 'success' if success else 'failed',
                'pattern_tags': pattern_tags,
                'analysis_type': analysis_type
            }
            
        except Exception as e:
            print(f"Error in combination {comb_name}: {str(e)}")
            results[comb_name] = {
                'status': 'failed',
                'error': str(e)
            }
            
        finally:
            reset_analysis()
    
    return results


def scale_load_case(load_case, factor):
    """Create a scaled copy of a load case.
    
    Args:
        load_case: Original LoadCase object
        factor: Scaling factor
        
    Returns:
        LoadCase: New scaled load case
    """
    scaled_case = LoadCase(load_case.name, load_case.type, factor)
    scaled_case.direction = load_case.direction
    
    # Scale nodal loads
    for node, loads in load_case.nodal_loads.items():
        scaled_case.nodal_loads[node] = [v*factor for v in loads]
    
    # Scale member loads
    for load in load_case.member_loads:
        scaled_case.member_loads.append({
            'eleTags': load['eleTags'],
            'type': load['type'],
            'values': [v*factor for v in load['values']]
        })
    
    # Scale surface loads
    for load in load_case.surface_loads:
        scaled_case.surface_loads.append({
            'eleTags': load['eleTags'],
            'face': load['face'],
            'pressure': load['pressure']*factor
        })
    
    # Scale mass loads (NEW)
    for node, masses in load_case.mass_loads.items():
        scaled_case.mass_loads[node] = [m*factor for m in masses]
    
    return scaled_case

def apply_masses(loaded_nodes, m_1, mass_distribution=None):
    """
    Apply nodal masses correctly.
    
    Args:
        loaded_nodes: List of node tags to apply mass to
        m_1: Total mass to distribute
        mass_distribution: Optional list of factors for each node 
                         (default: equal distribution)
    """
    if mass_distribution is None:
        # Equal distribution by default
        mass_distribution = [1.0/len(loaded_nodes)] * len(loaded_nodes)
    
    if abs(sum(mass_distribution) - 1.0) > 1e-6:
        raise ValueError("Mass distribution factors must sum to 1.0")
    
    for n, factor in zip(loaded_nodes, mass_distribution):
        # Apply mass in X, Y directions (Z typically for vertical structures)
        # Rotational masses typically zero for lumped mass model
        ops.mass(n, m_1*factor, m_1*factor, 0, 0, 0, 0)


# Analysis Functions - Improved with better recorder handling
def ensure_output_dir():
    """Ensure output directory exists."""
    os.makedirs('FGU_RC3DF_files', exist_ok=True)

def run_gravity(steps=10, comb_name=None):
    """Run gravity analysis with dynamic recorder naming."""
    ensure_output_dir()
    
    reaction_file = f"Gravity_Reactions_{comb_name}.out" if comb_name else "Gravity_Reactions.out"
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    
    ops.recorder('Node', '-file', reaction_path,
                '-time', '-node', *list(range(1,5)), '-dof', *list(range(1,7)), 'reaction')

    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1/steps)
    ops.analysis('Static')
    ops.record()
    
    ok = ops.analyze(steps)
    
    if ok == 0:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} completed successfully")
        return True
    else:
        print(f"Gravity analysis {'for ' + comb_name if comb_name else ''} failed")
        return False

def run_modal(n_evs=3, comb_name=None):
    """
    Runs Modal analysis with support for load combinations.
    
    Args:
        n_evs (int): Number of eigenvalues to compute
        comb_name (str): Name of the load combination (for output files)
        
    Returns:
        np.array: Array of eigenvalues
    """
    ensure_output_dir()
    
    # Create unique recorder names if combination is specified
    if comb_name:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec1_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec2_{comb_name}.out'),
            os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVec3_{comb_name}.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', f'ModalAnalysis_EigenVal_{comb_name}.out')
    else:
        eigen_files = [
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec1.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec2.out'),
            os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVec3.out')
        ]
        eigenval_file = os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out')

    # Set up recorders for each mode
    for i, file in enumerate(eigen_files[:n_evs], 1):
        ops.recorder('Node', '-file', file,
                    '-node', *list(range(5,9)), '-dof', 1, 2, f'eigen {i}')

    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-12, 25, 0, 2)
    ops.algorithm('Newton')
    ops.analysis('Transient')  # Needed for eigen commands
    
    # Compute eigenvalues
    lamda = np.array(ops.eigen(n_evs))
    
    # Write eigenvalues to file
    with open(eigenval_file, "w") as eig_file:
        eig_file.write("lambda omega period frequency\n")
        for l in lamda:
            omega = l**0.5
            period = 2*np.pi/omega
            freq = omega/(2*np.pi)
            eig_file.write(f"{l:2.6e} {omega:2.6e} {period:2.6e} {freq:2.6e}\n")

    ops.record()
    print(f"Modal analysis {'for ' + comb_name if comb_name else ''} completed")
    return lamda

def run_pushover(m_1, steps=10000, direction='X', comb_name=None):
    """Run pushover analysis with dynamic recorder naming."""
    ensure_output_dir()
    
    reaction_file = f"Pushover_Horizontal_Reactions{direction}_{comb_name}.out" if comb_name else f"Pushover_Horizontal_Reactions{direction}.out"
    disp_file = f"Pushover_Story_Displacement{direction}_{comb_name}.out" if comb_name else f"Pushover_Story_Displacement{direction}.out"
    
    reaction_path = os.path.join('FGU_RC3DF_files', reaction_file)
    disp_path = os.path.join('FGU_RC3DF_files', disp_file)

    d_o_f = 1 if direction == 'X' else 2
    phi = 1.0
    
    ops.recorder('Node', '-file', reaction_path,
                '-time', '-node', *list(range(1,5)), '-dof', d_o_f, 'reaction')
    ops.recorder('Node', '-file', disp_path,
                '-time', '-node', *list(range(5,9)), '-dof', d_o_f, 'disp')

    # Create lateral load pattern
    pattern_tag = 100 if comb_name else 2  # Use high tag for combination cases
    ops.pattern('Plain', pattern_tag, 1)
    for node in range(5, 9):
        if direction == 'X':
            ops.load(node, m_1*phi, 0, 0, 0, 0, 0)
        else:
            ops.load(node, 0, m_1*phi, 0, 0, 0, 0)

    step = 1.0e-05
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 0.000001, 100)
    ops.algorithm('NewtonLineSearch', True, 0.8, 1000, 0.1, 10.0)
    ops.integrator('DisplacementControl', 5, d_o_f, step)
    ops.analysis('Static')
    ops.record()

    ok = ops.analyze(steps)
    
    if ok == 0:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} completed successfully')
        return True
    else:
        print(f'Pushover Analysis in {direction} {"for " + comb_name if comb_name else ""} failed')
        return False

def run_time_history(direction='X', g_motion_id=1, scaling_id=1,
                    lamda=1.0, acc_file='FGU_RC3DF_files/acc_1.txt',
                    comb_name=None, analysis_duration_ratio=0.29):
    """
    Runs Time history analysis with support for load combinations.
    
    Args:
        direction (str): Direction of excitation ('X' or 'Y')
        g_motion_id (int): Ground motion identifier
        scaling_id (int): Scaling factor identifier
        lamda (float): Scaling factor for ground motion
        acc_file (str): Path to acceleration file
        comb_name (str): Name of the load combination
        analysis_duration_ratio (float): Fraction of full duration to analyze (0-1)
        
    Returns:
        bool: True if analysis succeeded, False otherwise
    """
    ensure_output_dir()
    
    # Create output file names
    if comb_name:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                   f'TimeHistory_Reactions_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                 f'TimeHistory_Acceleration_{direction}_{g_motion_id}_{scaling_id}_{comb_name}.out')
    else:
        reaction_file = os.path.join('FGU_RC3DF_files', 
                                    f'TimeHistory_Reactions_{direction}_{g_motion_id}.out')
        disp_file = os.path.join('FGU_RC3DF_files', 
                               f'TimeHistory_Displacement_{direction}_{g_motion_id}.out')
        accel_file = os.path.join('FGU_RC3DF_files', 
                                f'TimeHistory_Acceleration_{direction}_{g_motion_id}.out')

    # Determine DOF and damping parameters
    dof = 1 if direction == 'X' else 2
    
    # Get modal properties (assuming first mode dominates)
    try:
        eigenvals = np.loadtxt(os.path.join('FGU_RC3DF_files', 'ModalAnalysis_EigenVal.out'), 
                             skiprows=1)
        omega = eigenvals[0,1] if direction == 'X' else eigenvals[2,1]
    except:
        print("Warning: Could not read eigenvalue file, using default omega=2π")
        omega = 2*np.pi  # Fallback value
    
    xi = 0.05  # Damping ratio (5%)
    alpha_M = 0.0       # Mass proportional damping
    beta_K = 2*xi/omega # Stiffness proportional damping
    
    # Set up recorders
    ops.recorder('Node', '-file', reaction_file,
                '-time', '-node', *list(range(1,5)), '-dof', dof, 'reaction')
    ops.recorder('Node', '-file', disp_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'disp')
    ops.recorder('Node', '-file', accel_file,
                '-time', '-node', *list(range(5,9)), '-dof', dof, 'accel')

    # Load acceleration time history
    accelerogram = np.loadtxt(acc_file)
    dt = 0.02  # Time step in acceleration file
    n_steps = len(accelerogram)
    
    # Analysis parameters
    tol = 1.0e-6
    max_iter = 500
    analysis_dt = 0.01  # Analysis time step (should be ≤ dt/2 for accuracy)

    # Define time series and pattern
    ops.timeSeries('Path', 2, '-dt', dt, '-values', *accelerogram, '-factor', lamda)
    ops.pattern('UniformExcitation', 3, dof, '-accel', 2)
    
    # Analysis configuration
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', tol, max_iter, 0, 2)
    ops.algorithm('Newton')
    ops.integrator('Newmark', 0.5, 0.25)
    ops.rayleigh(alpha_M, beta_K, 0.0, 0.0)
    ops.analysis('Transient')

    # Run analysis
    print(f"Running Time-History analysis (λ={lamda}) {'for ' + comb_name if comb_name else ''}")
    start_time = time.time()
    
    ok = 0
    current_time = ops.getTime()
    final_time = n_steps * dt
    target_time = analysis_duration_ratio * final_time
    
    while ok == 0 and current_time < target_time:
        ok = ops.analyze(1, analysis_dt)
        current_time = ops.getTime()
        
        # Optional: Print progress
        if int(current_time/dt) % 100 == 0:
            print(f"Time: {current_time:.2f}s ({current_time/target_time:.1%})")
    
    elapsed_time = time.time() - start_time
    
    if ok == 0:
        print(f"Time-History completed in {elapsed_time:.2f}s")
        return True
    else:
        print(f"Time-History failed at {current_time:.2f}s")
        return False
    
def reset_analysis():
    """Reset the analysis state."""
    ops.setTime(0.0)
    ops.loadConst()
    ops.remove('recorders')
    ops.wipeAnalysis()


# Define load cases with mass information
load_cases = {
    "Dead": LoadCase("Dead", "gravity", 1.0),
    "Live": LoadCase("Live", "gravity", 1.0),
    "WindX": LoadCase("WindX", "lateral", 1.0, direction="X"),
    "WindY": LoadCase("WindY", "lateral", 1.0, direction="Y"),
    "SeismicX": LoadCase("SeismicX", "lateral", 1.0, direction="X"),
    "SeismicY": LoadCase("SeismicY", "lateral", 1.0, direction="Y"),
    "Mass": LoadCase("Mass", "gravity", 1.0)  # Special case for mass only
}

# Define load combinations (ASD or LRFD format) including mass combinations
load_combinations = {
    # Standard load combinations
    "Comb1": [("Dead", 1.2), ("Live", 1.6)],
    "Comb2": [("Dead", 1.2), ("Live", 1.0), ("WindX", 1.3)],
    "Comb3": [("Dead", 1.2), ("Live", 1.0), ("WindY", 1.3)],
    "Comb4": [("Dead", 1.2), ("Live", 0.5), ("SeismicX", 1.0)],
    "Comb5": [("Dead", 1.2), ("Live", 0.5), ("SeismicY", 1.0)],
    
    # Mass combinations for modal and dynamic analysis
    "Mass_Comb1": [("Mass", 1.0), ("Dead", 1.0)],  # For modal analysis
    "Mass_Comb2": [("Mass", 1.0), ("Dead", 1.0), ("Live", 0.3)],  # For seismic analysis
}

# Example Usage
# Define your specific load cases
load_cases["Dead"].add_nodal_load(5, [0, 0, -20*KN, 0, 0, 0])
load_cases["Dead"].add_nodal_load(6, [0, 0, -20*KN, 0, 0, 0])
load_cases["Dead"].add_nodal_load(7, [0, 0, -20*KN, 0, 0, 0])
load_cases["Dead"].add_nodal_load(8, [0, 0, -20*KN, 0, 0, 0])

load_cases["Dead"].add_member_load([1, 2, 3], 'beamUniform', [0, -5*KN_m, 0])

# Example surface load on elements 10-12, face 2 with 3 kN/m² pressure
load_cases["Dead"].add_surface_load([10, 11, 12], face=2, pressure=3*KN_m)

load_cases["Live"].add_nodal_load(5, [0, 0, -10*KN, 0, 0, 0])
load_cases["Live"].add_nodal_load(6, [0, 0, -10*KN, 0, 0, 0])
load_cases["Live"].add_nodal_load(7, [0, 0, -10*KN, 0, 0, 0])
load_cases["Live"].add_nodal_load(8, [0, 0, -10*KN, 0, 0, 0])

load_cases["WindX"].add_nodal_load(5, [5*KN, 0, 0, 0, 0, 0])
load_cases["WindX"].add_nodal_load(6, [5*KN, 0, 0, 0, 0, 0])
load_cases["WindX"].add_nodal_load(7, [5*KN, 0, 0, 0, 0, 0])
load_cases["WindX"].add_nodal_load(8, [5*KN, 0, 0, 0, 0, 0])

# Add mass to nodes (NEW)
mass_per_node = 1000.0  # kg
for node in [5, 6, 7, 8]:
    load_cases["Mass"].add_mass_load(node, [mass_per_node, mass_per_node, 0, 0, 0, 0])

# Run analyses for all combinations
results = run_analysis_for_combinations(load_combinations, load_cases)

# Run modal analysis with mass combinations
mass_combs = {k:v for k,v in load_combinations.items() if 'Mass' in k}
modal_results = run_analysis_for_combinations(mass_combs, load_cases, "modal")

# Run time history with mass combinations
th_results = run_analysis_for_combinations(
    {"Mass_Comb2": load_combinations["Mass_Comb2"]}, 
    load_cases,
    "time_history"
)
