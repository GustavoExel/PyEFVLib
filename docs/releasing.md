# Release and conda-forge preparation

The proposed next version is 1.1.0. No external release or upload has been performed. Pixi manages environments; conda-forge publishes packages through reviewed feedstocks.

1. Run the checks in `development.md`, review the changes, and test the consuming simulator in its own updated environment. Confirm the version and responsible maintainers; `lucscst` in the recipe must accept maintenance responsibilities.
2. Create an approved release and a publicly downloadable source archive. The upstream package already exists on PyPI: do not upload there without ownership and coordination.
3. `recipe/meta.yaml` is a local-build recipe. Before submission, replace `source: path: ..` with the release URL and the actual checksum:

```yaml
source:
  url: https://github.com/ThermoPhase-FCSRG/PyEFVLib/archive/refs/tags/v{{ version }}.tar.gz
  sha256: ACTUAL_SHA256_OF_THE_DOWNLOADED_RELEASE_ARCHIVE
```

This checksum is a placeholder, not a valid value. Calculate the hash of the downloaded archive with `Get-FileHash -Algorithm SHA256 <archive>` in PowerShell. Never submit a local path or placeholder checksum.

4. Test with `conda build recipe -c conda-forge` in an environment containing conda-build; repeat using the release source. A recipe file alone is not proof of a successful conda build.
5. Check whether a feedstock already exists and coordinate with its maintainers. Otherwise submit the reviewed recipe to `conda-forge/staged-recipes`. Publication follows review and CI, not a Pixi command.

Follow the [official conda-forge contribution process](https://conda-forge.org/docs/maintainer/adding_pkgs/).
