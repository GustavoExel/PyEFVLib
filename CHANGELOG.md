# Changelog

## 1.1.0 (unreleased)

- Modern setuptools metadata, Python >=3.11 and updated scientific-library ranges.
- Reproducible Pixi environments for current and minimum dependencies.
- NumPy object-alias replacement and SciPy sparse-coordinate compatibility.
- Backward-compatible property access for historical applications and the reservoir simulator.
- Meshio CellBlock support, topological XDMF dimensions and multi-block region tags.
- Deterministic grid stencils, clear missing-CSR-entry errors and ndarray dense conversion.
- Correct evaluation of named mathematical functions in boundary/initial expressions.
- Optional output paths and independent default numerical settings.
- Regression tests and Windows/Linux CI configuration.
- Local conda recipe and documented release/submission requirements.

The obsolete CgnsSaver imports in archived applications were removed or replaced with meshio. CGNS export remains subject to meshio's format and cell-type support. PETSc and the full research-workspace examples have not been validated by the default test suite.
