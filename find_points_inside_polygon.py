
import random
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt

# Step 1: Define 3D nodes (z = 0 for 2D area computation)
nodes_3d = [
    (0, 0, 0),
    (4, 0, 0),
    (4, 3, 0),
    (2, 5, 0),
    (0, 3, 0)
]

# Step 2: Project to 2D (XY plane)
nodes_2d = [(x, y) for x, y, z in nodes_3d]

# Step 3: Create Polygon and compute area
polygon = Polygon(nodes_2d)
area = polygon.area
print(f"Area of the closed polygon: {area}")

# Step 4: Generate a random point
# x_rand = random.uniform(-1, 6)
# y_rand = random.uniform(-1, 6)
x_rand = 2
y_rand=0
random_point = Point(x_rand, y_rand)
print(f"Random point: ({x_rand:.2f}, {y_rand:.2f})")

# Step 5: Check location of the point
if polygon.contains(random_point):
    status = "INSIDE"
elif polygon.touches(random_point):
    status = "ON THE EDGE"
else:
    status = "OUTSIDE"
print(f"The random point is {status} the polygon.")

# Step 6: Plot everything
x, y = zip(*nodes_2d + [nodes_2d[0]])  # Close the polygon loop
plt.figure(figsize=(6, 6))
plt.plot(x, y, 'b-o', label='Polygon')
plt.fill(x, y, alpha=0.2, color='blue')

# Plot the random point
color = 'go' if status == "ON THE EDGE" else ('ro' if status == "INSIDE" else 'rx')
plt.plot(x_rand, y_rand, color, markersize=10, label='Random Point')

# Annotate
plt.title(f'Random Point is {status}')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
plt.grid(True)
plt.legend()
plt.show()
