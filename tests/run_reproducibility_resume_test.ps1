$ErrorActionPreference = 'Stop'
$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')

if (-not (Test-Path tests/output/reproducibility_replay_a/integration/provenance/integration.analysis_status.tsv)) {
    throw 'Resume test requires completed replay A integration evidence'
}

$resumeCommand = "cd $repoLinux && nextflow run main.nf -profile test -work-dir tests/work/reproducibility_replay_a/integration -with-trace tests/output/reproducibility_resume_unchanged_trace.tsv --assay INTEGRATION --integration_manifest tests/fixtures/integration/manifest.full.tsv --fasta tests/fixtures/integration/reference.fa --gtf tests/fixtures/integration/reference.gtf --outdir tests/output/reproducibility_replay_a --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED -resume 2>&1"
$resumeOutput = & wsl.exe bash -lc $resumeCommand
$resumeExit = $LASTEXITCODE
$resumeOutput | Set-Content -Path tests/output/reproducibility_resume_unchanged.log
if ($resumeExit -ne 0) { throw 'Unchanged resume execution failed' }
if (($resumeOutput -join "`n") -notmatch 'cached:') { throw 'Unchanged resume did not report cached process reuse' }

$originalTrace = Import-Csv -Delimiter "`t" tests/output/reproducibility_resume_unchanged_trace.tsv
$originalHash = ($originalTrace | Where-Object { $_.name -like '*RUN_REGULATORY_INTEGRATION*' } | Select-Object -Last 1).hash
$changedCommand = "cd $repoLinux && nextflow run main.nf -profile test -work-dir tests/work/reproducibility_replay_a/integration -with-trace tests/output/reproducibility_resume_changed_trace.tsv --assay INTEGRATION --integration_manifest tests/fixtures/integration/manifest.full.tsv --fasta tests/fixtures/integration/reference.fa --gtf tests/fixtures/integration/reference.gtf --outdir tests/output/reproducibility_resume_changed --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED --distal_max_distance 50000 -resume 2>&1"
$changedOutput = & wsl.exe bash -lc $changedCommand
$changedExit = $LASTEXITCODE
$changedOutput | Set-Content -Path tests/output/reproducibility_resume_changed.log
if ($changedExit -ne 0) { throw 'Changed-parameter resume execution failed' }
$changedTrace = Import-Csv -Delimiter "`t" tests/output/reproducibility_resume_changed_trace.tsv
$changedHash = ($changedTrace | Where-Object { $_.name -like '*RUN_REGULATORY_INTEGRATION*' } | Select-Object -Last 1).hash
if (-not $originalHash -or -not $changedHash -or $originalHash -eq $changedHash) {
    throw 'Scientifically relevant parameter did not invalidate the integration task hash'
}

Write-Output "Resume verification passed: unchanged cache reused; integration task hash changed from $originalHash to $changedHash."
