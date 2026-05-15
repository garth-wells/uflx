"""Linear form using symfem."""

import symfem
import sympy

mesh_points = [[x, y] for y in [0, 1] for x in [0, 0.3, 1]]
mesh_triangles = [
    [0, 1, 3],
    [1, 4, 3],
    [1, 2, 4],
    [2, 5, 4],
]

"""
import matplotlib.pyplot as plt

for t in mesh_triangles:
    plt.plot([mesh_points[i][0] for i in t + [t[0]]], [mesh_points[i][1] for i in t + [t[0]]], "-")
plt.plot([i[0] for i in mesh_points], [i[1] for i in mesh_points], "ro")
plt.show()
"""

print(mesh_points)

e = symfem.create_element("triangle", "Lagrange", 1)

vector = [0 for _ in mesh_points]

f = symfem.symbols.x[0]

for t in mesh_triangles:
    triangle = symfem.create_reference("triangle", [mesh_points[i] for i in t])
    functions = e.map_to_cell(triangle.vertices)
    print(f"local vector for {t}")

    for v_dof, v in zip(t, functions):
        entry = (v * f).integral(triangle)
        print(float(entry), end=" ")
        vector[v_dof] += float(entry)
    print()
print(vector)
