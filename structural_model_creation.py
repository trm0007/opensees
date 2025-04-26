    
import json
import numpy as np
import math
import openseespy.opensees as ops

def create_structural_model(structure_data_path, member_section_mapping_path, numIntgrPts = 8):
    # Load structure data from JSON files
    with open(structure_data_path, 'r') as f:
        structure_data = json.load(f)
    
    with open(member_section_mapping_path, 'r') as f:
        member_section_mapping = json.load(f)
    
    # Extract nodes and members from structure data
    nodes_data = structure_data["nodes"]
    members_data = structure_data["members"]
    
    # Convert nodes to numpy array for easier access
    # Format: node_id (1-based) maps to index (0-based)
    nodes = np.zeros((len(nodes_data), 3))
    for node in nodes_data:
        node_id = node["id"]
        nodes[node_id-1] = [node["x"], node["y"], node["z"]]
    
    # Define section tags
    section_tags = {
        "section1": 1,
        "section2": 2,
        "section3": 3,
        "section4": 4,
    }
    

    
    # Integration parameters for nonlinearBeamColumn
      # Number of integration points
    
    # Lists to store results
    shapes = {}
    memberLengths = []
    
    # Process each member
    for n, mbr_data in enumerate(members_data):
        mbr_id = mbr_data["id"]
        mbr_name = mbr_data["name"]
        node_i_id = mbr_data["start_node_id"] 
        node_j_id = mbr_data["end_node_id"]
        
        # Create a list of the node IDs for this member
        mbr = [node_i_id, node_j_id]
        
        # Determine the transformation type based on member name
        if mbr_name.startswith("cz"):  # Column
            transType = 'PDelta'
        else:  # Beam
            transType = 'Linear'
        
        # Element and transformation tags
        transTag = mbr_id
        eleTag = mbr_id
        
        # Get the section tag for this member
        section_name = member_section_mapping[mbr_name]
        secTag = section_tags[section_name]
        
        # Get node coordinates (0-based indexing for nodes array)
        node_i = mbr[0] - 1  # Convert to 0-based indexing for numpy array
        node_j = mbr[1] - 1
        
        ix = nodes[node_i, 0]
        iy = nodes[node_i, 1]
        iz = nodes[node_i, 2]
        jx = nodes[node_j, 0]
        jy = nodes[node_j, 1]
        jz = nodes[node_j, 2]
        
        # Calculate vector components along the member
        dx = jx - ix
        dy = jy - iy
        dz = jz - iz
        length = math.sqrt(dx**2 + dy**2 + dz**2)
        memberLengths.append(length)
        
        # Calculate local x unit vector
        if length > 0:
            local_x_unit = np.array([dx, dy, dz]) / length
        else:
            local_x_unit = np.array([1.0, 0.0, 0.0])  # Default if length is zero
        
        # Create a vector in the global XZ plane (typically [0,1,0] works as a reference)
        # This is a reference vector not parallel to local_x_unit
        reference_vector = np.array([0.0, 1.0, 0.0])
        
        # Check if local_x_unit is parallel to reference_vector
        if abs(np.dot(local_x_unit, reference_vector) / np.linalg.norm(reference_vector)) > 0.99:
            # If too parallel, choose a different reference vector
            reference_vector = np.array([1.0, 0.0, 0.0])
        
        # Calculate local y unit vector using cross product
        local_z_unit = np.cross(local_x_unit, reference_vector)
        local_z_unit = local_z_unit / np.linalg.norm(local_z_unit)
        
        # Ensure local z is perpendicular to local x
        local_z_unit = local_z_unit - np.dot(local_z_unit, local_x_unit) * local_x_unit
        local_z_unit = local_z_unit / np.linalg.norm(local_z_unit)
        
        # Finally, ensure y is perpendicular to both x and z
        local_y_unit = np.cross(local_z_unit, local_x_unit)
        
        # Define the geometric transformation
        # vecxz components define the local x-z plane
        vecxz = local_z_unit
        print(f'element:{eleTag}=[{float(vecxz[0])}, {float(vecxz[1])}, {float(vecxz[2])}]')
        
        # Define geometric transformation
        ops.geomTransf(transType, transTag, float(vecxz[0]), float(vecxz[1]), float(vecxz[2]))
        
        # Define the nonlinear beam column element
        ops.element('nonlinearBeamColumn', eleTag, node_i_id, node_j_id, numIntgrPts, 
                   secTag, transTag, '-integration', 'Lobatto')
    
    return shapes, memberLengths

# # Example usage:
# if __name__ == "__main__":
#     structure_data_path = "structure_data_json.json"
#     member_section_mapping_path = "member_section_mapping.json"
    
#     # Initialize OpenSees model
#     ops.wipe()
#     ops.model('basic', '-ndm', 3, '-ndf', 6)  # 3D model with 6 DOF per node
    
#     # Create nodes from structure data
#     with open(structure_data_path, 'r') as f:
#         structure_data = json.load(f)
    
#     for node in structure_data["nodes"]:
#         node_id = node["id"]
#         x, y, z = node["x"], node["y"], node["z"]
#         ops.node(node_id, x, y, z)
    
#     # Define materials and sections (simplified example)
#     # In a real application, you would define proper sections based on your requirements
#     E = 200000.0  # Young's modulus (MPa)
#     nu = 0.3      # Poisson's ratio
#     G = E / (2 * (1 + nu))  # Shear modulus
    
#     # Example section properties - would be different for columns and beams
#     # Section 1 (for columns)
#     ops.section('Elastic', 1, E, 0.01, 1.0e-4, 1.0e-4, G, 1.0e-6)
    
#     # Section 2 (for beams)
#     ops.section('Elastic', 2, E, 0.008, 8.0e-5, 8.0e-5, G, 8.0e-7)
    
#     # Create the structural model elements
#     shapes, memberLengths = create_structural_model(structure_data_path, member_section_mapping_path)
    
#     print("Model created successfully!")
#     print(f"Number of members: {len(memberLengths)}")
#     print(f"Member lengths: {memberLengths}")
