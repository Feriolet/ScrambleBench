# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project does **NOT** adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html), because of the v prefix.

## [Unreleased]
## [v0.1.1]
- None

## [v0.1.0] - 2026-07-17

### Added

- Add Checkpoints to run ScrambleBench parts
- Add easydock PLIF
- User can generate PDF report based on their ScrambleBench parameter runs.
- Unify parameters as one YAML input files
- Validate of user's input as guardrails
- Add genbench3d `Validity3D` breakdown.
- Add examples to run `ScrambleBench`
- Add a subset of `POKMOL` dataset as an example without protein-ligand compelex input.
- logger using python's `logging` module
- Add docstring and typehints as documentation
- Add Changelog.md (lol)

### Changed

- Refactor codes and architecture to lean towards OOP because of the type of project as a workflow
- User can now use YAML as argparse input instead of the previous version where user need to hardcode parameters and arguments
- Integrate diversity metric as part of the pipeline using a wrapper
- Used latest rdkit version `2026.2`. Does not have significant impact as of now.
- Edit README.md
- Execute redocking using `/tmp` folder as default to hide intermediate files.
- Instead of hardcoding virtual hit criteria, user can now define their own criteria using `df.query()` syntax

### Removed

- Determination of the best virtual hit compound
- Code to regenerate Main and Supplementary Figures from manuscript. Please refer to `v0.0.1` version.

## [v0.0.1] - 2026-2-13

### Added

- Barebone code to run ScrambleBench
- Manuscript Figures as part of manuscript submission to journal