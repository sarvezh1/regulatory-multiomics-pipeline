param([switch]$ContinueAfterFirst)
$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param([Parameter(Mandatory=$true)][string]$Command,[Parameter(Mandatory=$true)][string[]]$Arguments)
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Command exited with status $LASTEXITCODE" }
}

function Remove-ReplayTree {
    param([Parameter(Mandatory=$true)][string]$Path)
    $resolvedParent = [System.IO.Path]::GetFullPath((Split-Path -Parent $Path))
    $expectedParent = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'output'))
    $expectedWork = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'work'))
    if ($resolvedParent -ne $expectedParent -and $resolvedParent -ne $expectedWork) {
        throw "Refusing to remove a replay path outside tests/output or tests/work"
    }
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
}

Invoke-Checked python @('tests/create_synthetic_fixtures.py')
Invoke-Checked python @('tests/create_atac_fixtures.py')
Invoke-Checked python @('tests/create_wgbs_fixtures.py')
Invoke-Checked python @('tests/create_chip_fixtures.py')
Invoke-Checked python @('tests/create_integration_fixtures.py')

$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')

function Invoke-Replay {
    param([Parameter(Mandatory=$true)][string]$Label)
    $outWin = Join-Path $PSScriptRoot "output/reproducibility_replay_$Label"
    $workWin = Join-Path $PSScriptRoot "work/reproducibility_replay_$Label"
    $reportWin = Join-Path $PSScriptRoot "output/reproducibility_report_$Label"
    Remove-ReplayTree $outWin
    Remove-ReplayTree $workWin
    Remove-ReplayTree $reportWin

    $out = "tests/output/reproducibility_replay_$Label"
    $work = "tests/work/reproducibility_replay_$Label"
    $runs = @(
        "cd $repoLinux && nextflow run main.nf -profile test -work-dir $work/rna --assay RNA --input tests/fixtures/generated/raw_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir $out --reference_scope synthetic --run_de false",
        "cd $repoLinux && nextflow run main.nf -profile test -work-dir $work/atac --assay ATAC --atac_input_route raw --input tests/fixtures/atac/raw_samplesheet.tsv --fasta tests/fixtures/atac/reference.fa --gtf tests/fixtures/atac/reference.gtf --blacklist tests/fixtures/atac/blacklist.bed --outdir $out --reference_scope synthetic --run_da false",
        "cd $repoLinux && nextflow run main.nf -profile test -work-dir $work/wgbs --assay WGBS --wgbs_input_route raw --input tests/fixtures/wgbs/raw_samplesheet.tsv --fasta tests/fixtures/wgbs/reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir $out --reference_scope synthetic --run_dm false",
        "cd $repoLinux && nextflow run main.nf -profile test -work-dir $work/chip --assay H3K27ac --chip_input_route raw --input tests/fixtures/chip/raw_single.tsv --fasta tests/fixtures/chip/reference.fa --gtf tests/fixtures/chip/reference.gtf --blacklist tests/fixtures/chip/blacklist.bed --outdir $out --reference_scope synthetic --chip_effective_genome_size 2000000 --run_chip_diff false",
        "cd $repoLinux && nextflow run main.nf -profile test -work-dir $work/integration --assay INTEGRATION --integration_manifest tests/fixtures/integration/manifest.full.tsv --fasta tests/fixtures/integration/reference.fa --gtf tests/fixtures/integration/reference.gtf --outdir $out --integration_mode both --integration_match_type EXACT_SAMPLE_MATCHED"
    )
    foreach ($run in $runs) { Invoke-Checked wsl.exe @('bash','-lc',$run) }

    $report = "tests/output/reproducibility_report_$Label"
    Invoke-Checked python @(
        'bin/build_run_report.py', '--results', $out, '--source-root', '.', '--output', $report,
        '--requested-assay', 'RNA,ATAC,WGBS,H3K27ac,INTEGRATION', '--run-mode', 'clean_room_synthetic_replay',
        '--nextflow-version', '26.04.6',
        '--input-record', 'tests/fixtures/generated/raw_samplesheet.tsv',
        '--input-record', 'tests/fixtures/atac/raw_samplesheet.tsv',
        '--input-record', 'tests/fixtures/wgbs/raw_samplesheet.tsv',
        '--input-record', 'tests/fixtures/chip/raw_single.tsv',
        '--input-record', 'tests/fixtures/integration/manifest.full.tsv',
        '--reference', 'RNA_FASTA=tests/fixtures/generated/reference.fa',
        '--reference', 'RNA_GTF=tests/fixtures/generated/reference.gtf',
        '--reference', 'ATAC_FASTA=tests/fixtures/atac/reference.fa',
        '--reference', 'ATAC_GTF=tests/fixtures/atac/reference.gtf',
        '--reference', 'ATAC_BLACKLIST=tests/fixtures/atac/blacklist.bed',
        '--reference', 'WGBS_FASTA=tests/fixtures/wgbs/reference.fa',
        '--reference', 'WGBS_GTF=tests/fixtures/wgbs/reference.gtf',
        '--reference', 'CHIP_FASTA=tests/fixtures/chip/reference.fa',
        '--reference', 'CHIP_GTF=tests/fixtures/chip/reference.gtf',
        '--reference', 'CHIP_BLACKLIST=tests/fixtures/chip/blacklist.bed',
        '--reference', 'INTEGRATION_FASTA=tests/fixtures/integration/reference.fa',
        '--reference', 'INTEGRATION_GTF=tests/fixtures/integration/reference.gtf',
        '--parameter-file', 'tests/fixtures/reproducibility_parameters.json',
        '--tool-containers', 'assets/tool_containers.tsv',
        '--timestamp-utc', '2026-10-13T00:00:00+00:00'
    )
    Invoke-Checked python @('bin/validate_run_report.py', $report)
}

if ($ContinueAfterFirst) {
    if (-not (Test-Path tests/output/reproducibility_replay_a/integration/provenance/integration.analysis_status.tsv)) {
        throw 'Cannot continue: completed replay A evidence is absent'
    }
    Invoke-Checked python @(
        'bin/build_run_report.py', '--results', 'tests/output/reproducibility_replay_a', '--source-root', '.', '--output', 'tests/output/reproducibility_report_a',
        '--requested-assay', 'RNA,ATAC,WGBS,H3K27ac,INTEGRATION', '--run-mode', 'clean_room_synthetic_replay', '--nextflow-version', '26.04.6',
        '--input-record', 'tests/fixtures/generated/raw_samplesheet.tsv', '--input-record', 'tests/fixtures/atac/raw_samplesheet.tsv',
        '--input-record', 'tests/fixtures/wgbs/raw_samplesheet.tsv', '--input-record', 'tests/fixtures/chip/raw_single.tsv',
        '--input-record', 'tests/fixtures/integration/manifest.full.tsv', '--reference', 'RNA_FASTA=tests/fixtures/generated/reference.fa',
        '--reference', 'RNA_GTF=tests/fixtures/generated/reference.gtf', '--reference', 'ATAC_FASTA=tests/fixtures/atac/reference.fa',
        '--reference', 'ATAC_GTF=tests/fixtures/atac/reference.gtf', '--reference', 'ATAC_BLACKLIST=tests/fixtures/atac/blacklist.bed',
        '--reference', 'WGBS_FASTA=tests/fixtures/wgbs/reference.fa', '--reference', 'WGBS_GTF=tests/fixtures/wgbs/reference.gtf',
        '--reference', 'CHIP_FASTA=tests/fixtures/chip/reference.fa', '--reference', 'CHIP_GTF=tests/fixtures/chip/reference.gtf',
        '--reference', 'CHIP_BLACKLIST=tests/fixtures/chip/blacklist.bed', '--reference', 'INTEGRATION_FASTA=tests/fixtures/integration/reference.fa',
        '--reference', 'INTEGRATION_GTF=tests/fixtures/integration/reference.gtf', '--parameter-file', 'tests/fixtures/reproducibility_parameters.json',
        '--tool-containers', 'assets/tool_containers.tsv', '--timestamp-utc', '2026-10-13T00:00:00+00:00'
    )
    Invoke-Checked python @('bin/validate_run_report.py', 'tests/output/reproducibility_report_a')
} else {
    Invoke-Replay 'a'
}
Invoke-Replay 'b'
Invoke-Checked python @('tests/compare_reproducibility_replays.py','tests/output/reproducibility_replay_a','tests/output/reproducibility_replay_b','tests/output/reproducibility_report_a','tests/output/reproducibility_report_b')

Invoke-Checked powershell.exe @('-NoProfile','-ExecutionPolicy','Bypass','-File','tests/run_reproducibility_resume_test.ps1')

Write-Output 'Reproducibility clean-room replay, deterministic comparison, and resume checks passed.'
