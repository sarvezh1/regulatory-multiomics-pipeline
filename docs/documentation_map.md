# Documentation map

## Start here

- `README.md` gives the public overview, requirements, quick start, scope, and validation boundaries.
- `CHANGELOG.md` summarizes public versions and notable capabilities.
- `docs/usage.md` contains commands for every implemented raw, processed, and integration route.
- `docs/output_reference.md` defines the main outputs, identifiers, coordinate systems, missing values, effect directions, and component statuses.

## Scientific design and inputs

- `docs/scientific_scope.md`: frozen scientific scope and reference policy.
- `docs/dataset_audit.md`: audited studies, accessions, roles, and limitations.
- `docs/tool_stack.md`: executed tool versions and container policy.
- `docs/validation_plan.md`: validation strategy and test boundaries.

## Implemented branches

- `docs/rna_branch.md`
- `docs/atac_branch.md`
- `docs/wgbs_branch.md`
- `docs/h3k27ac_branch.md`
- `docs/integration_branch.md`

`schemas/` contains the machine-readable contract definitions. Branch documents explain scientific behavior; the schemas remain authoritative for required fields.

## Reproducibility and public release

- `docs/reporting.md`: unified report and provenance architecture.
- `docs/reproducibility.md`: deterministic content, replay, resume, and validated backend.
- `docs/public_data_sources.md`: synthetic, bounded public, processed public, and external-validation data policy.
- `docs/public_file_inventory.md`: future repository inclusion and exclusion policy.
- `docs/external_validation.md`: cCRE concordance and CRISPR benchmark methods and limitations.

## Validation records

- [RNA validation](validation/rna.md)
- [ATAC validation](validation/atac.md)
- [WGBS validation](validation/wgbs.md)
- [H3K27ac validation](validation/h3k27ac.md)
- [Integration validation](validation/integration.md)
- [Reproducibility validation](validation/reproducibility.md)
- [External validation](validation/external.md)

Heavy genuine-data and Docker validation are local-only; hosted CI performs structural, hygiene, contract, and lightweight unit checks.
