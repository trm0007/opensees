import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def create_dome_mesh(radius=1.0, height=0.0, num_segments=8):
    """Create a uniformly meshed dome and return mesh information."""
    # Create the dome vertices
    theta = np.linspace(0, 2*np.pi, num_segments, endpoint=False)
    phi = np.linspace(0, np.pi/2, num_segments//2)
    
    vertices = []
    for p in phi:
        for t in theta:
            x = float(radius * np.sin(p) * np.cos(t))
            y = float(radius * np.sin(p) * np.sin(t))
            z = float(height * np.cos(p))
            vertices.append([x, y, z])
    
    # Add the apex point
    vertices.append([0.0, 0.0, float(height)])
    apex_index = len(vertices) - 1
    
    # Create the base vertices
    for t in theta:
        x = float(radius * np.cos(t))
        y = float(radius * np.sin(t))
        vertices.append([x, y, 0.0])
    
    # Create the mesh elements in CCW order
    mesh_elements = {}
    all_nodes = {i: [float(v[0]), float(v[1]), float(v[2])] for i, v in enumerate(vertices)}
    
    # Generate quadrilateral elements
    quad_count = 0
    for i in range(num_segments//2 - 1):
        for j in range(num_segments):
            v1 = i * num_segments + j
            v2 = i * num_segments + (j + 1) % num_segments
            v3 = (i + 1) * num_segments + (j + 1) % num_segments
            v4 = (i + 1) * num_segments + j
            mesh_name = f"quad_{quad_count}"
            mesh_elements[mesh_name] = [v1, v2, v3, v4]
            quad_count += 1
    
    # Generate triangular elements
    tri_count = 0
    top_ring_start = (num_segments//2 - 1) * num_segments
    for j in range(num_segments):
        v1 = top_ring_start + j
        v2 = top_ring_start + (j + 1) % num_segments
        v3 = apex_index
        mesh_name = f"tri_{tri_count}"
        mesh_elements[mesh_name] = [v1, v2, v3]
        tri_count += 1
    
    # Generate base elements
    base_start = apex_index + 1
    for j in range(num_segments):
        v1 = base_start + j
        v2 = base_start + (j + 1) % num_segments
        v3 = top_ring_start + (j + 1) % num_segments
        v4 = top_ring_start + j
        mesh_name = f"base_quad_{j}"
        mesh_elements[mesh_name] = [v1, v2, v3, v4]
    
    return list(mesh_elements.keys()), all_nodes, mesh_elements

def plot_dome_with_labels(meshes, all_nodes, mesh_elements):
    """Plot the dome with node indices and mesh names labeled."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot each mesh element with center label
    for mesh_name, nodes in mesh_elements.items():
        verts = [all_nodes[n] for n in nodes]
        
        # Calculate center point for mesh label
        center = np.mean(verts, axis=0)
        
        if len(nodes) == 4:  # Quadrilateral
            poly = Poly3DCollection([verts], alpha=0.3, linewidths=1, edgecolor='k')
        else:  # Triangle
            poly = Poly3DCollection([verts], alpha=0.3, linewidths=1, edgecolor='k')
        
        poly.set_facecolor(np.random.rand(3))
        ax.add_collection3d(poly)
        
        # Add mesh name label at center
        ax.text(center[0], center[1], center[2], mesh_name, 
                color='black', fontsize=8, ha='center', va='center')
    
    # Plot nodes with indices
    xs, ys, zs = zip(*all_nodes.values())
    ax.scatter(xs, ys, zs, color='red', s=50)
    
    # Add node index labels
    for node_idx, (x, y, z) in all_nodes.items():
        ax.text(x, y, z, f'n{node_idx}', 
                color='blue', fontsize=8, ha='center', va='bottom')
    
    # Adjust view and labels
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.set_title('Dome Mesh with Node Indices and Element Labels', pad=20)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(0, 1.2)
    
    # Add legend
    ax.text2D(0.05, 0.95, "Blue: Node indices (n0, n1,...)\nBlack: Mesh names (quad_0, tri_0,...)", 
              transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.show()

# Create the dome mesh
meshes, all_nodes, mesh_elements = create_dome_mesh(num_segments=8)

# Print node and mesh information
print("All Nodes:")
for idx, coords in all_nodes.items():
    print(f"Node n{idx}: {coords}")

print("\nMesh Elements:")
for name, nodes in mesh_elements.items():
    print(f"{name}: Nodes {['n'+str(n) for n in nodes]}")

# Plot with labels
plot_dome_with_labels(meshes, all_nodes, mesh_elements)


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Polygon

class AdvancedStructuralMesher:
    def __init__(self):
        self.nodes = {}
        self.meshes = {}
        self.node_counter = 0
        self.mesh_counter = 0
    
    def add_node(self, x, y, z):
        """Add a node with uniform naming"""
        node_id = f"N{self.node_counter:04d}"
        self.nodes[node_id] = [float(x), float(y), float(z)]
        self.node_counter += 1
        return node_id
    
    def add_mesh(self, node_ids, mesh_type):
        """Add a mesh with uniform naming and CCW ordering"""
        # Ensure counter-clockwise ordering by checking normal direction
        if len(node_ids) >= 3:
            normal = self._calculate_normal([self.nodes[n] for n in node_ids])
            if normal[2] < 0:  # If normal points downward, reverse order
                node_ids = node_ids[::-1]
        
        mesh_id = f"M{self.mesh_counter:04d}_{mesh_type}"
        self.meshes[mesh_id] = {
            "nodes": node_ids,
            "type": mesh_type,
            "normal": normal if len(node_ids) >= 3 else [0,0,1]
        }
        self.mesh_counter += 1
        return mesh_id
    
    def _calculate_normal(self, vertices):
        """Calculate normal vector for a planar element"""
        if len(vertices) < 3:
            return [0, 0, 1]
        v1 = np.array(vertices[1]) - np.array(vertices[0])
        v2 = np.array(vertices[2]) - np.array(vertices[0])
        normal = np.cross(v1, v2)
        norm = np.linalg.norm(normal)
        return normal / norm if norm > 0 else normal
    
    def create_rectangular_surface(self, p1, p2, p3, p4, divisions_x, divisions_y, mesh_type):
        """Create uniformly meshed rectangular surface with CCW node ordering"""
        # Convert points to numpy arrays
        p1, p2, p3, p4 = [np.array(p) for p in [p1, p2, p3, p4]]
        
        # Create parametric grid
        u = np.linspace(0, 1, divisions_x + 1)
        v = np.linspace(0, 1, divisions_y + 1)
        
        # Generate nodes
        node_grid = []
        for i in range(divisions_x + 1):
            row = []
            for j in range(divisions_y + 1):
                # Bilinear interpolation
                point = (1-v[j])*((1-u[i])*p1 + u[i]*p2) + v[j]*((1-u[i])*p4 + u[i]*p3)
                row.append(self.add_node(*point))
            node_grid.append(row)
        
        # Create quadrilateral meshes
        for i in range(divisions_x):
            for j in range(divisions_y):
                node_ids = [
                    node_grid[i][j],
                    node_grid[i+1][j],
                    node_grid[i+1][j+1],
                    node_grid[i][j+1]
                ]
                self.add_mesh(node_ids, mesh_type)
    
    def create_inclined_wall(self, base_line, height, angle_deg, axis='x', divisions_x=5, divisions_y=5, mesh_type="inclined"):
        """
        Create an inclined wall at specific angle with respect to an axis
        
        Parameters:
            base_line: [(x1,y1,z1), (x2,y2,z2)] - base line of the wall
            height: height of the wall
            angle_deg: inclination angle in degrees
            axis: axis to incline relative to ('x', 'y', or 'z')
            divisions_x: number of divisions along length
            divisions_y: number of divisions along height
            mesh_type: type of mesh to create
        """
        start, end = base_line
        length = np.linalg.norm(np.array(end) - np.array(start))
        angle_rad = np.radians(angle_deg)
        
        # Calculate direction vectors
        length_dir = (np.array(end) - np.array(start)) / length
        if axis.lower() == 'x':
            height_dir = np.array([np.cos(angle_rad), 0, np.sin(angle_rad)])
        elif axis.lower() == 'y':
            height_dir = np.array([0, np.cos(angle_rad), np.sin(angle_rad)])
        elif axis.lower() == 'z':
            height_dir = np.array([0, 0, 1])  # Vertical by default
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")
        
        # Create nodes grid
        node_grid = []
        for i in range(divisions_x + 1):
            length_ratio = i / divisions_x
            length_pos = np.array(start) + length_ratio * (np.array(end) - np.array(start))
            row = []
            for j in range(divisions_y + 1):
                height_ratio = j / divisions_y
                point = length_pos + height * height_ratio * height_dir
                row.append(self.add_node(*point))
            node_grid.append(row)
        
        # Create meshes
        for i in range(divisions_x):
            for j in range(divisions_y):
                node_ids = [
                    node_grid[i][j],
                    node_grid[i+1][j],
                    node_grid[i+1][j+1],
                    node_grid[i][j+1]
                ]
                self.add_mesh(node_ids, mesh_type)
    
    def plot_structure(self):
        """Plot the structure with node and mesh labels"""
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111, projection='3d')
        
        # Color mapping for different mesh types
        type_colors = {
            "floor": [0.7, 0.7, 1.0],
            "wall": [1.0, 0.7, 0.7],
            "inclined": [0.7, 1.0, 0.7],
            "default": [0.8, 0.8, 0.8]
        }
        
        # Plot each mesh
        for mesh_id, mesh_data in self.meshes.items():
            verts = [self.nodes[n] for n in mesh_data["nodes"]]
            center = np.mean(verts, axis=0)
            
            # Get color based on mesh type
            mesh_type = mesh_data["type"]
            color = type_colors.get(mesh_type, type_colors["default"])
            
            # Create polygon
            poly = Poly3DCollection([verts], alpha=0.6)
            poly.set_facecolor(color)
            poly.set_edgecolor('k')
            poly.set_linewidth(0.5)
            ax.add_collection3d(poly)
            
            # Add mesh label at center (slightly offset for visibility)
            label_pos = center + np.array(mesh_data["normal"]) * 0.1
            ax.text(label_pos[0], label_pos[1], label_pos[2], mesh_id, 
                   color='black', fontsize=6, ha='center', va='center')
        
        # Plot all nodes with labels
        for node_id, coords in self.nodes.items():
            ax.scatter(*coords, color='red', s=20)
            ax.text(coords[0], coords[1], coords[2], node_id, 
                   color='blue', fontsize=6, ha='center', va='center')
        
        # Set view properties
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        ax.set_title('Structural Meshes with Precise Inclined Walls', pad=20)
        
        # Create legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', label='Nodes (red)',
                      markerfacecolor='red', markersize=8),
            plt.Rectangle((0,0),1,1, fc=type_colors["floor"], alpha=0.6, ec='k', label='Floor'),
            plt.Rectangle((0,0),1,1, fc=type_colors["wall"], alpha=0.6, ec='k', label='Wall'),
            plt.Rectangle((0,0),1,1, fc=type_colors["inclined"], alpha=0.6, ec='k', label='Inclined')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Equal aspect ratio
        ax.set_box_aspect([1, 1, 1])
        plt.tight_layout()
        plt.show()
    
    def print_mesh_data(self):
        """Print all meshes with CCW node ordering"""
        print("\nMesh Data (Counter-Clockwise Ordering):")
        print("=====================================")
        print(f"{'Mesh ID':<10} | {'Type':<8} | {'Normal Vector':<25} | Nodes")
        print("-"*80)
        
        for mesh_id, mesh_data in sorted(self.meshes.items()):
            normal_str = f"[{mesh_data['normal'][0]:.2f}, {mesh_data['normal'][1]:.2f}, {mesh_data['normal'][2]:.2f}]"
            nodes_str = ", ".join(mesh_data["nodes"])
            print(f"{mesh_id:<10} | {mesh_data['type']:<8} | {normal_str:<25} | {nodes_str}")
        
        print("\nSummary:")
        print(f"Total Nodes: {len(self.nodes)}")
        print(f"Total Meshes: {len(self.meshes)}")
        print(f"Mesh Types: {set(m['type'] for m in self.meshes.values())}")

# Example usage
if __name__ == "__main__":
    mesher = AdvancedStructuralMesher()
    
    # Create floor slab (10x6 units)
    p1 = [0, 0, 0]
    p2 = [10, 0, 0]
    p3 = [10, 6, 0]
    p4 = [0, 6, 0]
    mesher.create_rectangular_surface(p1, p2, p3, p4, 5, 3, "floor")
    
    # Create vertical wall (10 units long, 4 units high)
    mesher.create_inclined_wall(
        base_line=[[0, 0, 0], [10, 0, 0]],
        height=4,
        angle_deg=90,  # Vertical
        axis='y',
        divisions_x=5,
        divisions_y=4,
        mesh_type="wall"
    )
    
    # Create wall inclined 30 degrees from vertical (X-axis)
    mesher.create_inclined_wall(
        base_line=[[0, 6, 0], [6, 6, 0]],
        height=3,
        angle_deg=30,
        axis='x',
        divisions_x=4,
        divisions_y=3,
        mesh_type="inclined"
    )
    
    # Create wall inclined 45 degrees from vertical (Y-axis)
    mesher.create_inclined_wall(
        base_line=[[8, 0, 0], [8, 6, 0]],
        height=3,
        angle_deg=45,
        axis='y',
        divisions_x=3,
        divisions_y=3,
        mesh_type="inclined"
    )
    
    # Print and plot results
    mesher.print_mesh_data()
    mesher.plot_structure()
