# AGENTS.md

Guidance for agents (Claude Code and others) working in this repository.

## What this is

UFLx is an experimental, minimal reimplementation of UFL (Unified Form Language): a
symbolic language embedded in Python for writing finite-element variational forms
(e.g. `inner(grad(v), grad(u)) * dx`). It does not own meshes, boundary markers,
coefficient values, or assembly — those belong to a surrounding solver environment
(e.g. DOLFIN/DOLFINx). See [LANGUAGE.md](LANGUAGE.md) for the full language spec and
[DESIGN_CHOICES.md](DESIGN_CHOICES.md) for naming/API conventions — read both before
adding to the object model.

## Commands

Install for development (editable):
```
pip install -e ".[ci,lint,docs]"
```

Lint (must pass with zero diff/errors, matches CI):
```
ruff check .
ruff format --check .          # use `ruff format .` to auto-fix
mypy . --exclude external
mypy external/codegeneration
mypy external/basix_uflx
```

Run the core test suite:
```
pytest test/
pytest test/test_forms.py::test_name    # single test
pytest -n auto test/                    # parallel, matches CI (needs pytest-xdist)
```

The `external/` packages (`external/codegeneration`, `external/basix_uflx`) are
separate installable packages excluded from the main `uflx` package; install and test
them independently:
```
pip install external/codegeneration external/basix_uflx
pytest external/codegeneration/test external/basix_uflx/test
```

Build docs:
```
cd doc && make html
```

## Architecture

Expression trees (the symbolic core) are not plain Python object graphs — they are
backed by an explicit DAG in `uflx/graphs/` (`Graph` subclasses `networkx.DiGraph`;
nodes are `GraphNode`s exposing `successors` and `init_args`). Whole-tree rewrites
(substituting terminals, extracting real/imaginary parts, applying push-forward/
pull-back maps) go through `uflx/graphs/algorithms/replace.py`, which walks the graph
and calls `reconstruct_node` (`uflx/graphs/algorithms/reconstruct.py`) to rebuild any
node whose successors changed, by re-invoking `node.__class__(*args)` with replaced
`init_args`. Anything that needs to transform an expression (`complex.py`, `maps.py`,
`geometry.py`) is built on this `replace` primitive rather than ad hoc tree-walking —
follow that pattern for new transformations.

The object model is a layered chain of abstract base classes, each file named after
its plural concept with `Abstract*` base classes (concrete classes are defined by
consumers, e.g. `test/conftest.py`'s `LagrangeElement`, or the `basix_uflx` extension):

- `entities.py` — `AbstractEntity`: topological mesh entities (point/interval/
  triangle/.../hexahedron), defined recursively via their sub-entities.
- `finite_elements.py` — `AbstractFiniteElement` / `AbstractReferenceMappedFiniteElement`:
  what basis functions look like on a cell; depends on `entities` and `maps`.
- `maps.py` — `AbstractReferenceMap` (e.g. `IdentityReferenceMap`): push-forward/
  pull-back between reference and physical cells, implemented as graph rewrites.
- `domains.py` — `AbstractDomain`: the integration domain (topological/geometric
  dimension + coordinate element); the actual mesh stays external to UFLx.
- `function_spaces.py` — `AbstractFunctionSpace`: standard FE spaces (domain +
  element), constant spaces (shape + scalar type), or non-FE spaces (domain + shape,
  no element). Do not construct `Argument`/`Coefficient` directly from an element —
  they must come from a `FunctionSpace`.
- `functions.py` / `basis_functions.py` — `AbstractFunction` and reference/physical
  basis functions evaluated at points.
- `expressions.py` — `AbstractExpression`: the base of every symbolic node
  (`BinaryOperator`, `UnaryOperator`, terminals). Carries value shape, free indices,
  domain, scalar type as static, extensible attributes (unlike legacy UFL's fixed
  attribute set).
- `operators.py`, `tensors.py`, `points.py`, `geometry.py`, `complex.py` — concrete
  operators (`inner`, `grad`, ...), tensor/vector construction, point sets, geometric
  quantities (spatial coordinates, Jacobians, ...), and complex-number support
  (`re`/`im` via the `ComplexValued` protocol), all built as `AbstractExpression`
  subclasses / graph rewrites over them.
- `integrals.py` — measures (`dx`, `ds`, `dS`) and `Integral` (`expr * measure`); a
  `Form` is a sum/combination of integrals whose integrand must be scalar with no
  free indices.

Symbolic (Gateaux) differentiation and expression/integral transformation between
domain configurations are the two core transformation procedures UFLx provides on top
of this object model (see LANGUAGE.md §4) — both are lazy, expressed as graph
operators rather than eagerly evaluated.

## Python style

Follows FEniCS project conventions (see [dolfinx](https://github.com/FEniCS/dolfinx)'s
`AGENTS.md`), adapted for this pure-Python, MIT-licensed package:

- **Formatting/linting**: `ruff check` and `ruff format --check`, configured in
  `pyproject.toml`. Line length 100, 4-space indent. Rule set includes pydocstyle
  (`D`, Google convention), pycodestyle, pyflakes, isort, pyupgrade,
  flake8-import-conventions, NumPy-specific rules, and a few more (`RUF`, `FLY`,
  `LOG`, `ISC`). Run the formatter/linter locally before calling a change done —
  don't rely on CI to catch formatting.
- **Import order** (ruff's isort, default grouping): future → standard-library →
  third-party (`networkx`, `pytest`) → first-party (`uflx`) → local-folder.
- **Docstrings**: Google style (`Args:`, `Returns:`, etc.), required on modules,
  classes and methods (`D` rules are enforced, `__init__.py` is exempted from unused
  imports only — not docstrings).
- **Type hints**: used throughout and checked with `mypy` (`[tool.mypy]` in
  `pyproject.toml`); not currently enforced by a ruff annotation rule, but new public
  functions/methods should still be annotated to match the rest of the codebase.
- **File header**: most files under `uflx/` (the core object-model layer —
  `entities.py`, `finite_elements.py`, `domains.py`, `function_spaces.py`,
  `functions.py`, `expressions.py`, `operators.py`, `integrals.py`,
  `basis_functions.py`, `__init__.py`) start with:

  ```python
  # Copyright (C) <year> <author(s)>
  #
  # This file is part of UFLx (https://www.fenicsproject.org)
  #
  # SPDX-License-Identifier:    MIT
  ```

  followed by a module docstring. Some newer/infrastructure files (`uflx/graphs/`
  and its `algorithms/` submodule, `maps.py`, `tensors.py`, `points.py`,
  `geometry.py`, `complex.py`, `utils.py`) currently omit it — match the header
  convention for new files in the core object-model layer; add your name/year rather
  than replacing existing authors when making a substantive change to a file that has
  one.

## Conventions (see DESIGN_CHOICES.md)

- Module files are named as the plural of what they contain (e.g. `domains.py`)
  precisely so the module (`uflx.domains`) and a same-named factory function
  (`uflx.domain`) can coexist.
- A class's factory/constructor function has the same name as the class, lowercased
  (mid-word capitals become `_x`, e.g. a hypothetical `FooX` class pairs with
  `foo_x()`).
- `__init__` methods take no default values; if defaults are useful, add a separate
  factory function instead of defaulting the constructor.
