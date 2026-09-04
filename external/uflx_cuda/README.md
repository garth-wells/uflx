# uflx-cuda

Experimental CUDA source generation for geometry transforms identified from
UFLx forms.

The initial implementation recognises scalar Poisson stiffness forms on a
three-dimensional tetrahedral domain:

```python
form = inner(grad(u), grad(v)) * dx
source = generate_cuda_geometry_kernel(form)
```

It emits a standalone templated `__global__` kernel that computes the weighted
six-component symmetric tensor

```text
G = weight * abs(det(J)) * inv(J) * inv(J).T
```

at quadrature points. Its runtime interface follows the geometry kernel in
`dolfinx-gpu-solvers`: one block processes one selected cell, coordinate dofs
are staged in dynamic shared memory, and one thread computes each quadrature
point. The caller supplies global coordinates, the geometry dofmap, the full
coordinate-element tabulation `[4, nq, ncdofs]`, quadrature weights, and the
selected entity list.

The output layout is component-major within each cell:
`[cell, component, quadrature_point]`, with components ordered
`(G00, G01, G02, G11, G12, G22)`. Launch with at least `nq` threads per block
and `3 * ncdofs * sizeof(T)` bytes of dynamic shared memory.

Current restrictions are deliberate: a single `dx` integral, scalar
`inner(grad(u), grad(v))`, a tetrahedral cell, and geometric dimension three.
