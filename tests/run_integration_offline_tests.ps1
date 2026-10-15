$ErrorActionPreference = 'Stop'
python tests/create_integration_fixtures.py
python -m unittest tests/test_integration_contracts.py
$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')
$base = "cd $repoLinux && nextflow run main.nf -profile test --assay INTEGRATION --fasta tests/fixtures/integration/reference.fa --gtf tests/fixtures/integration/reference.gtf -resume"
wsl.exe bash -lc "$base --integration_manifest tests/fixtures/integration/manifest.full.tsv --outdir tests/output/integration_synthetic --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED"
if ($LASTEXITCODE -ne 0) { throw 'Full integration workflow failed' }
wsl.exe bash -lc "$base --integration_manifest tests/fixtures/integration/manifest.partial.tsv --outdir tests/output/integration_partial --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED"
if ($LASTEXITCODE -ne 0) { throw 'Partial integration workflow failed' }
wsl.exe bash -lc "$base --integration_manifest tests/fixtures/integration/manifest.mismatch.tsv --outdir tests/output/integration_mismatch --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED"
if ($LASTEXITCODE -ne 0) { throw 'Contrast-mismatch integration workflow failed' }
python tests/assert_integration_outputs.py --full tests/output/integration_synthetic --partial tests/output/integration_partial --mismatch tests/output/integration_mismatch
