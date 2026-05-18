# Minimal UFLx language definition

This note describes the smallest useful language model behind legacy UFLx. It intentionally avoids implementation details, compiler algorithms, and broad redesign ideas.

## 1. What UFLx is

UFLx is a symbolic language embedded in Python for writing finite-element variational forms.

A UFLx program is ordinary Python code that constructs symbolic objects. The important result is usually a `Form`, such as a bilinear form, linear form, or functional.

UFLx is not a full PDE language. It does not own the actual mesh, boundary markers, numerical values of coefficients, or assembly procedure. Those belong to the surrounding problem-solving environment, such as DOLFIN or DOLFINx.

In one sentence:

> UFLx is a typed symbolic expression language for finite-element forms.

## 2. The minimal object model

The minimal language needs only a small number of object kinds.

### Finite elements

A finite element describes the local finite-dimensional field used on a cell.

At minimum, an element is defined by:

- the cell type,
- the element family,
- the polynomial degree or equivalent order information,
- the value shape,
- and, for mapped elements, the required pullback or mapping behaviour.

Examples:

```python
P1 = FiniteElement("Lagrange", triangle, 1, ())
vP2 = FiniteElement("Lagrange", triangle, 2, (2,))
```

The element says what kind of discrete field is being represented. It does not store actual coefficient values.

### Domains and meshes

A domain or mesh is the integration domain of a form. For a minimal language definition, a domain only needs to provide:

- topological dimension,
- geometric dimension,
- and a coordinate element, if geometry is represented using a finite element.

The actual mesh is still external to UFLx.

### Function spaces

A function space combines a domain with a finite element.
Conceptually:

```python
V = FunctionSpace(mesh, element)
```

In older examples, `Argument` and `Coefficient` are sometimes created directly from an element. This is not allowed.

### Terminals

Terminals are the leaves of UFLx expression trees.

The essential terminals are:

- `Argument`, representing an arbitrary basis/test/trial function,
- `Coefficient`, representing a user-supplied discrete function,
- literals such as numbers,
- and geometric quantities such as coordinates, normals, Jacobians, cell diameter, facet area, etc.

Example:

```python
v = TestFunction(V, name="v")
u = TrialFunction(V, name="u")
f = Coefficient(V, name="f")
```

> Legacy note: There is no `Constant` in UFLx as a means to represent literal value which is not known at compile time.

UFLx supports named symbolic terminals. This allows for better error messages and better postprocessing of expressions. In addition, it improves the stability of the expression tree.

> Legacy note: In legacy UFL, symbolic objects do not have names. It requires expensive and error prone renumbering of terminals in expressions before lowering and code generation.

### Expressions

An expression is either a terminal or an operator applied to other expressions.

Examples of operators are addition, multiplication, indexing, `grad`, `div`, `dot`, `inner`, `jump`, and `avg`.

Expressions are symbolic. They are not evaluated when they are created.

A minimal expression carries these static properties:

- tensor shape,
- free indices,
- index dimensions,
- associated domain or domains,
- scalar type (real or complex),
- continuity information,
- (possibly) physical dimension.

These properties are enough to reject many meaningless expressions before lowering or code generation.

### Measures

A measure describes where and how an integrand is integrated.

TODO.

The standard measures are:

- `dx` for cell integrals,
- `ds` for exterior-facet integrals,
- `dS` for interior-facet integrals.

Examples:

```python
dx
ds
dS
```

### Integrals

An integral is made by multiplying a scalar expression by a measure.

Examples:

```python
v * f * dx
inner(grad(v), grad(u)) * dx
jump(v) * dS
```

A valid integral has a scalar-valued integrand with no unresolved free indices.

### Forms

A form is a multi-linear functional $a(u, v, w, \ldots): V_1 \times V_2 \times V_3 \times \ldots \to \mathbb{R}$.

Arguments determine the arity of a form. A form with no arguments is a functional. A form with one argument is a linear form. A form with two arguments is a bilinear form.

> Legacy note: In contrast to the legacy UFL, in UFLx a form does not need to be defined as a sum of integrals. For example, the following is a valid UFLx form: `a = (u * dx) * (v * dx)`, or `a = PointEvaluation(u, (0.5, 0.5)) * PointEvaluation(v, (0.5, 0.5))`.

Examples:

```python
a = inner(grad(v), grad(u)) * dx
L = v * f * dx
F = a - L
```


## 3. Minimal validity rules

The language is small, but it is typed.

- Addition and subtraction require compatible shapes and compatible free indices.
- Multiplication follows scalar, tensor, matrix-vector, matrix-matrix, and indexed-expression rules.
- Tensor contractions such as `dot` and `inner` require compatible tensor shapes.
- Spatial derivatives require a spatial domain.
- Facet quantities such as `FacetNormal` are only meaningful in facet integration contexts.
- Interior-facet expressions may use side restrictions:
    ```python
    u("+")
    u("-")
    ```
- Operators such as `jump` and `avg` are shorthand for combinations of such restrictions.
- An integral is valid only when the final integrand is scalar and has no free indices.
- A form is valid when all of its integrals are valid and their domains, arguments, and coefficients are mutually consistent.

## 4. What transformation procedures does the UFLx core provide

There are two key transformation procedures that the UFLx core provides:

1. Automatic/symbolic differentiation.

   Similar to UFL, UFLx provides a symbolic differentiation via a directional (Gateaux) derivative. The derivative of an Expression, Integral or a Form is lazy/eagerly evaluated. That means the derivative is stored as an expression tree operator, and the actual differentiation is performed when the derivative is used in a context where it needs to be evaluated, such as lowering or code generation.

   Example:

   ```python
    a = inner(grad(v), grad(u)) * dx
    da = derivative(a, u, du)
    ```

2. Expression and integral transformation between different configurations.

   > Legacy note: In legacy UFL, there is a very central perspective on the configuration of where objects are defined. Usually, all expressions are defined in the physical configuration, and UFL applies pullbacks to the entire expression tree.

   In UFLx, configuration where expressions are defined are a consequence of the domain. Since the domain is part of the Expression's static attributes, it is propagated through the expression tree based on the language rules.

   TODO.

## 5. What belongs to the language

The minimal language includes:

- symbolic finite-element objects,
- symbolic domains/function spaces,
- symbolic terminals,
- tensor-valued expressions,
- index notation,
- differential operators,
- DG side restrictions,
- integration measures,
- integrals,
- and forms.

This is enough to describe ordinary weak forms such as mass, stiffness, elasticity, Stokes, mixed Poisson, and DG forms.

## 6. What does not belong to the minimal language

The following are not primitive language concepts:

- actual mesh connectivity,
- actual coefficient vectors,
- basis tabulation,
- quadrature-loop generation,
- form compilation,
- assembly,
- solver setup,
- boundary-condition application,
- expression simplification strategy,
- code generation.

These are consumers or transformations of the symbolic language, not the core language itself.

## 7. Extensibility of the attributes

> Legacy note: In legacy UFL, the set of static attributes is fixed and hardcoded. This led to many language features being hacked in (e.g. complex numbers support, dual spaces, etc.).

In UFLx, the set of static attributes is extensible. This allows for cleaner language extensions in the future, such as support for physical dimensions, dual spaces, error estimates, custom quadrature rules, custom integrability/continuity rules, etc.

## 8. Examples

### Poisson equation

```python
coordinate_element = FiniteElement("Lagrange", triangle, 1, shape=(2,))
element = FiniteElement("Lagrange", triangle, 1)

mesh = Mesh(coordinate_element)
V = FunctionSpace(mesh, element)

u = TrialFunction(V, name="u")
v = TestFunction(V, name="v")
f = Coefficient(V, name="f")

a = inner(grad(v), grad(u)) * dx
L = v * f * dx
```

