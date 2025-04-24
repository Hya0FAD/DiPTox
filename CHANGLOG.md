# Changelog

## [0.11.0] - 2025-03-28
### Added
- Split original `remove_salts` functionality into separate operations (`remove_salts`, `remove_solvents`, `remove_mixtures`)
- Added consistency check for chemical composition after `remove_salts`
- Defined 9 new removable solvents with post-removal consistency verification
- Implemented customizable solvent list management

### Changed
- Redefined 22 types of removable salts
- Renamed salt list management function (`manage_default_salt`)