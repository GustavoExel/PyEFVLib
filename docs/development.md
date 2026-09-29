# Development

This fork requires Python 3.11 or newer. Install [Pixi](https://pixi.prefix.dev/latest/installation/) and run from the repository root:

```sh
pixi install
pixi run check
pixi run -e minimum check
pixi run build
pixi run check-dist
```

Commit `pixi.lock` with dependency changes. It resolves Windows and Linux environments. The default environment uses current dependencies; `minimum` exercises the supported lower bounds. A Windows local run does not validate Linux execution.

NumPy, SciPy, meshio, h5py and xmltodict are runtime dependencies. pandas is not required by the core library. Matplotlib supports examples; PETSc remains optional and is not tested in the default environment.
The Pixi development environment also includes pandas for historical converters and research scripts.

## Using the fork in the reservoir simulator

Create an isolated, updated environment instead of upgrading the old Python 3.9 environment in place. Replace the simulator's old PyEFVLib dependency with the following Pixi dependency during development:

```toml
[pypi-dependencies]
PyEFVLib = { path = "../PyEFVLib", editable = true }
```

For a pyproject.toml manifest, use `[tool.pixi.pypi-dependencies]`. Do not keep another entry for the old PyEFVLib. After release, pin a release or immutable Git commit instead of relying on a sibling directory.

Both `propertyData.get(region_handle, property_name)` and `propertyData.get(property_name, handle=region_handle)` are supported. This preserves historical applications and the simulator's newer API. Updating this library does not automatically remove restrictions in the simulator or torch-flash: resolve and test those environments separately.

Tests cover control-volume conservation, linear-gradient reconstruction, property access, expressions, dense/CSR systems and meshio output. Compilation checks do not prove numerical correctness of every research script in `workspace`; these may require additional inputs and dependencies.
