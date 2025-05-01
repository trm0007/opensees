import math
import numpy as np

class SurfaceLoad6DOF:
    # Class variables
    oneOverRoot3 = 1.0 / math.sqrt(3.0)
    GsPts = np.array([
        [-oneOverRoot3, -oneOverRoot3],
        [oneOverRoot3, -oneOverRoot3],
        [oneOverRoot3, oneOverRoot3],
        [-oneOverRoot3, oneOverRoot3]
    ])
    
    # Constants
    SL_NUM_NODE = 4
    SL_NUM_NDF = 6  # Changed from 3 to 6 (3 translations + 3 rotations)
    SL_NUM_DOF = SL_NUM_NODE * SL_NUM_NDF
    
    def __init__(self, tag, nd1, nd2, nd3, nd4, pressure):
        self.tag = tag
        self.myExternalNodes = np.array([nd1, nd2, nd3, nd4])
        self.my_pressure = pressure
        self.mLoadFactor = 1.0
        
        # Initialize vectors and matrices
        self.g1 = np.zeros(3)  # Still 3D vectors for geometry
        self.g2 = np.zeros(3)
        self.myNhat = np.zeros(3)
        self.myNI = np.zeros(self.SL_NUM_NODE)
        self.dcrd1 = np.zeros(3)
        self.dcrd2 = np.zeros(3)
        self.dcrd3 = np.zeros(3)
        self.dcrd4 = np.zeros(3)
        
        # Node pointers
        self.theNodes = [None] * self.SL_NUM_NODE
        
        # Stiffness and force vectors
        self.tangentStiffness = np.zeros((self.SL_NUM_DOF, self.SL_NUM_DOF))
        self.internalForces = np.zeros(self.SL_NUM_DOF)
        self.theVector = np.zeros(self.SL_NUM_DOF)
    
    def setDomain(self, domain):
        for i in range(self.SL_NUM_NODE):
            self.theNodes[i] = domain.getNode(self.myExternalNodes[i])
            if self.theNodes[i] is None:
                print(f"Node {self.myExternalNodes[i]} not found in domain")
                return
        
        self.dcrd1 = np.array(self.theNodes[0].getCrds())
        self.dcrd2 = np.array(self.theNodes[1].getCrds())
        self.dcrd3 = np.array(self.theNodes[2].getCrds())
        self.dcrd4 = np.array(self.theNodes[3].getCrds())
    
    def updateBase(self, Xi, Eta):
        oneMinusEta = 1.0 - Eta
        onePlusEta = 1.0 + Eta
        oneMinusXi = 1.0 - Xi
        onePlusXi = 1.0 + Xi
        
        # Calculate vectors g1 and g2 (geometry still in 3D)
        self.g1 = 0.25 * (oneMinusEta * (self.dcrd2 - self.dcrd1) + onePlusEta * (self.dcrd3 - self.dcrd4))
        self.g2 = 0.25 * (onePlusXi * (self.dcrd3 - self.dcrd2) + oneMinusXi * (self.dcrd4 - self.dcrd1))
        
        # Shape functions
        self.myNI = np.array([
            0.25 * oneMinusXi * oneMinusEta,
            0.25 * onePlusXi * oneMinusEta,
            0.25 * onePlusXi * onePlusEta,
            0.25 * oneMinusXi * onePlusEta
        ])
        
        # Normal vector as cross product of g1 and g2
        self.myNhat = np.array([
            self.g1[1] * self.g2[2] - self.g1[2] * self.g2[1],
            self.g1[2] * self.g2[0] - self.g1[0] * self.g2[2],
            self.g1[0] * self.g2[1] - self.g1[1] * self.g2[0]
        ])
    
    def getResistingForce(self):
        self.internalForces = np.zeros(self.SL_NUM_DOF)
        
        # Loop over Gauss points
        for i in range(4):
            self.updateBase(self.GsPts[i, 0], self.GsPts[i, 1])
            
            # Loop over nodes
            for j in range(4):
                # Only apply forces to translational DOFs (0,1,2)
                # Rotational DOFs (3,4,5) remain zero for pressure load
                for k in range(3):
                    idx = j * self.SL_NUM_NDF + k
                    self.internalForces[idx] -= (self.mLoadFactor * self.my_pressure * 
                                              self.myNhat[k] * self.myNI[j])
        
        return self.internalForces
    
    def getEquivalentNodalLoads(self):
        """Calculate and return the equivalent nodal loads"""
        resisting_forces = self.getResistingForce()
        # Equivalent nodal loads are the negative of resisting forces
        equivalent_loads = -resisting_forces
        return equivalent_loads
    
    def printEquivalentNodalLoads(self):
        """Print the equivalent nodal loads in a formatted way"""
        loads = self.getEquivalentNodalLoads()
        print("\nEquivalent Nodal Loads for SurfaceLoad Element", self.tag)
        print("------------------------------------------------")
        print(f"Pressure: {self.my_pressure} (load factor: {self.mLoadFactor})")
        print("Node    FX          FY          FZ          MX          MY          MZ")
        print("----------------------------------------------------------------------------")
        
        for i in range(4):
            node_tag = self.myExternalNodes[i]
            start_idx = i * self.SL_NUM_NDF
            fx = loads[start_idx]
            fy = loads[start_idx+1]
            fz = loads[start_idx+2]
            mx = loads[start_idx+3]
            my = loads[start_idx+4]
            mz = loads[start_idx+5]
            print(f"{node_tag:<6}  {fx:10.4f}  {fy:10.4f}  {fz:10.4f}  {mx:10.4f}  {my:10.4f}  {mz:10.4f}")
        
        print("----------------------------------------------------------------------------\n")

# Example usage
if __name__ == "__main__":
    # Create a simple domain
    class Node:
        def __init__(self, tag, crds):
            self.tag = tag
            self.crds = crds
        
        def getCrds(self):
            return self.crds
    
    class Domain:
        def __init__(self):
            self.nodes = {}
        
        def addNode(self, node):
            self.nodes[node.tag] = node
        
        def getNode(self, tag):
            return self.nodes.get(tag)
    
    # Create a domain and add nodes
    domain = Domain()
    domain.addNode(Node(1, [0, 0, 0]))
    domain.addNode(Node(2, [1, 0, 0]))
    domain.addNode(Node(3, [1, 1, 0]))
    domain.addNode(Node(4, [0, 1, 0]))
    
    # Create a surface load element with pressure = 10.0
    pressure = 10.0
    surface_load = SurfaceLoad6DOF(1, 1, 2, 3, 4, pressure)
    surface_load.setDomain(domain)
    
    # Print equivalent nodal loads
    surface_load.printEquivalentNodalLoads()
    
    # Change load factor and print again
    surface_load.mLoadFactor = 2.0
    print("\nAfter changing load factor to 2.0:")
    surface_load.printEquivalentNodalLoads()

import math
import numpy as np

class TriSurfaceLoad6DOF:
    # Class variables
    oneOverRoot3 = 1.0 / math.sqrt(3.0)
    GsPts = np.array([[1.0/3]])  # Single integration point for triangular element
    
    # Constants
    SL_NUM_NODE = 3  # Triangular element
    SL_NUM_NDF = 6   # 6 DOF per node (3 translations + 3 rotations)
    SL_NUM_DOF = SL_NUM_NODE * SL_NUM_NDF
    
    def __init__(self, tag, nd1, nd2, nd3, pressure, rhoH=0.0):
        self.tag = tag
        self.myExternalNodes = np.array([nd1, nd2, nd3])
        self.my_pressure = pressure
        self.rhoH = rhoH
        self.mLoadFactor = 1.0
        
        # Initialize vectors and matrices
        self.g1 = np.zeros(3)  # Geometric vectors (3D)
        self.g2 = np.zeros(3)
        self.myNhat = np.zeros(3)
        self.myNI = np.zeros(self.SL_NUM_NODE)
        self.dcrd1 = np.zeros(3)
        self.dcrd2 = np.zeros(3)
        self.dcrd3 = np.zeros(3)
        
        # Node pointers
        self.theNodes = [None] * self.SL_NUM_NODE
        
        # Stiffness, mass, and force vectors
        self.tangentStiffness = np.zeros((self.SL_NUM_DOF, self.SL_NUM_DOF))
        self.mass = np.zeros((self.SL_NUM_DOF, self.SL_NUM_DOF))
        self.damp = np.zeros((self.SL_NUM_DOF, self.SL_NUM_DOF))
        self.internalForces = np.zeros(self.SL_NUM_DOF)
    
    def setDomain(self, domain):
        for i in range(self.SL_NUM_NODE):
            self.theNodes[i] = domain.getNode(self.myExternalNodes[i])
            if self.theNodes[i] is None:
                print(f"Node {self.myExternalNodes[i]} not found in domain")
                return
        
        self.dcrd1 = np.array(self.theNodes[0].getCrds())
        self.dcrd2 = np.array(self.theNodes[1].getCrds())
        self.dcrd3 = np.array(self.theNodes[2].getCrds())
    
    def updateBase(self, Xi, Eta):
        # Calculate vectors g1 and g2 (geometry still in 3D)
        self.g1 = self.dcrd2 - self.dcrd1
        self.g2 = self.dcrd3 - self.dcrd1
        
        # CST shape functions (triangular coordinates)
        self.myNI = np.array([
            Xi,
            Eta,
            1 - Xi - Eta
        ])
        
        # Normal vector as cross product of g1 and g2
        self.myNhat = np.array([
            self.g1[1] * self.g2[2] - self.g1[2] * self.g2[1],
            self.g1[2] * self.g2[0] - self.g1[0] * self.g2[2],
            self.g1[0] * self.g2[1] - self.g1[1] * self.g2[0]
        ])
        
        # Normalize (area is norm/2)
        self.myNhat = self.myNhat / 2
    
    def getResistingForce(self):
        self.internalForces = np.zeros(self.SL_NUM_DOF)
        
        # Loop over Gauss points (only 1 for triangle)
        for i in range(1):
            self.updateBase(self.GsPts[i, 0], self.GsPts[i, 0])
            
            # Loop over nodes
            for j in range(3):
                # Only apply forces to translational DOFs (0,1,2)
                # Rotational DOFs (3,4,5) remain zero for pressure load
                for k in range(3):
                    idx = j * self.SL_NUM_NDF + k
                    self.internalForces[idx] -= (self.mLoadFactor * self.my_pressure * 
                                              self.myNhat[k] * self.myNI[j])
        
        return self.internalForces
    
    def getMass(self):
        area = np.linalg.norm(self.myNhat)
        self.mass = np.zeros((self.SL_NUM_DOF, self.SL_NUM_DOF))
        
        if self.rhoH > 0:
            # Distribute mass equally to each node (lumped mass matrix)
            nodal_mass = self.rhoH * area / 3
            for i in range(self.SL_NUM_NODE):
                for j in range(3):  # Only diagonal terms for translational DOFs
                    idx = i * self.SL_NUM_NDF + j
                    self.mass[idx, idx] = nodal_mass
        
        return self.mass
    
    def getResistingForceIncInertia(self):
        self.internalForces = self.getResistingForce()
        
        if self.rhoH > 0:
            # Get nodal accelerations
            accel = np.zeros(self.SL_NUM_DOF)
            for i in range(self.SL_NUM_NODE):
                a = self.theNodes[i].getAccel()
                for j in range(3):  # Only translational DOFs
                    idx = i * self.SL_NUM_NDF + j
                    accel[idx] = a[j]
            
            # Add inertia forces
            self.internalForces -= self.getMass() @ accel
        
        return self.internalForces
    
    def getEquivalentNodalLoads(self):
        """Calculate and return the equivalent nodal loads"""
        resisting_forces = self.getResistingForce()
        # Equivalent nodal loads are the negative of resisting forces
        equivalent_loads = -resisting_forces
        return equivalent_loads
    
    def printEquivalentNodalLoads(self):
        """Print the equivalent nodal loads in a formatted way"""
        loads = self.getEquivalentNodalLoads()
        print("\nEquivalent Nodal Loads for TriSurfaceLoad Element", self.tag)
        print("------------------------------------------------")
        print(f"Pressure: {self.my_pressure} (load factor: {self.mLoadFactor})")
        print("Node    FX          FY          FZ          MX          MY          MZ")
        print("----------------------------------------------------------------------------")
        
        for i in range(3):
            node_tag = self.myExternalNodes[i]
            start_idx = i * self.SL_NUM_NDF
            fx = loads[start_idx]
            fy = loads[start_idx+1]
            fz = loads[start_idx+2]
            mx = loads[start_idx+3]
            my = loads[start_idx+4]
            mz = loads[start_idx+5]
            print(f"{node_tag:<6}  {fx:10.4f}  {fy:10.4f}  {fz:10.4f}  {mx:10.4f}  {my:10.4f}  {mz:10.4f}")
        
        print("----------------------------------------------------------------------------\n")

# Example usage
if __name__ == "__main__":
    # Create a simple domain
    class Node:
        def __init__(self, tag, crds):
            self.tag = tag
            self.crds = crds
            self.accel = np.zeros(3)
        
        def getCrds(self):
            return self.crds
        
        def getAccel(self):
            return self.accel
    
    class Domain:
        def __init__(self):
            self.nodes = {}
        
        def addNode(self, node):
            self.nodes[node.tag] = node
        
        def getNode(self, tag):
            return self.nodes.get(tag)
    
    # Create a domain and add nodes
    domain = Domain()
    domain.addNode(Node(1, [0, 0, 0]))
    domain.addNode(Node(2, [1, 0, 0]))
    domain.addNode(Node(3, [0, 1, 0]))
    
    # Create a triangular surface load element with pressure = 10.0
    pressure = 10.0
    rhoH = 2.0  # Mass per unit area
    tri_surface_load = TriSurfaceLoad6DOF(1, 1, 2, 3, pressure, rhoH)
    tri_surface_load.setDomain(domain)
    
    # Print equivalent nodal loads
    tri_surface_load.printEquivalentNodalLoads()
    
    # Set some accelerations to demonstrate inertia forces
    domain.nodes[1].accel = np.array([1.0, 0.0, 0.0])
    domain.nodes[2].accel = np.array([0.0, 1.0, 0.0])
    domain.nodes[3].accel = np.array([0.0, 0.0, 1.0])
    
    print("\nWith accelerations and mass (rhoH = 2.0):")
    forces_with_inertia = tri_surface_load.getResistingForceIncInertia()
    tri_surface_load.printEquivalentNodalLoads()
