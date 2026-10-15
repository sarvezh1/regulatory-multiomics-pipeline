$ErrorActionPreference = 'Stop'

python tests/validate_external_validation_outputs.py --root .
if ($LASTEXITCODE -ne 0) { throw "External-validation output validation failed" }

python -m unittest tests.test_external_validation -v
if ($LASTEXITCODE -ne 0) { throw "External-validation unit validation failed" }
