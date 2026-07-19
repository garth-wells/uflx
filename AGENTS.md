# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project overview

UFLx is an experimental, minimal reimplementation of UFL (Unified Form Language): a
symbolic language embedded in Python for writing finite-element variational forms. It
does not own meshes, boundary markers, coefficient values, or assembly — those belong
to a surrounding solver environment (e.g. DOLFIN/DOLFINx). See [LANGUAGE.md](LANGUAGE.md)
for the language model and [DESIGN_CHOICES.md](DESIGN_CHOICES.md) for naming/API
conventions that must be followed when adding code.

Requires Python >= 3.11. Runtime dependencies are `networkx` and `numpy`.

## Layout

- `uflx/` — the package. Notable modules: `expressions.py`, `operators.py`,
  `finite_elements.py`, `function_spaces.py`, `functions.py`, `domains.py`,
  `entities.py`, `integrals.py`, `geometry.py`, `maps.py`, `quadrature.py`,
  `scalars.py`, `points.py`, `complex.py`, and the `graphs/` subpackage (expression
  graph algorithms).
- `uflx/test_utils/` — shared test fixtures/helpers used by both `test/` and
  external test suites.
- `test/` — unit tests (pytest).
- `examples/` — small standalone example forms (e.g. `mass.py`, `stiffness.py`,
  `linear_form.py`).
- `external/` — external code (e.g. `codegeneration`) with its own tests, excluded
  from the main package build (see `pyproject.toml`).
- `doc/` — Sphinx documentation source.

## Build, lint, and test

```bash
pip install -e .[ci,lint,docs]

ruff check .
ruff format --check .
mypy uflx

pytest test/
```

CI (`.github/workflows/pythonapp.yml`) also builds and tests `external/codegeneration`
and builds the Sphinx docs. Some tests are gated on the `UFLX_ENABLE_TESTING=1`
environment variable, as set in CI.

## Conventions

- Follow [DESIGN_CHOICES.md](DESIGN_CHOICES.md): module file names are plural (e.g.
  `domains.py` module, `domain` function); a class's plain-lowercase constructor
  function (e.g. `finite_element` for `FiniteElement`) has no default arguments —
  add a separate convenience function if defaults are needed.
- Docstrings follow the Google convention (`ruff` pydocstyle rule enforced).
- Line length is 100 columns; run `ruff format` before committing.
- Keep this file (and `CLAUDE.md`, which mirrors it) in sync when project structure,
  build commands, or conventions change.
