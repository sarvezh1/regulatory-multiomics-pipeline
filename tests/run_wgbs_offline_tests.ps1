$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param([Parameter(Mandatory=$true)][string]$Command,[Parameter(Mandatory=$true)][string[]]$Arguments)
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Command exited with status $LASTEXITCODE" }
}

Invoke-Checked python @('tests/create_wgbs_fixtures.py')
Invoke-Checked python @('-m','unittest','tests/test_wgbs_contracts.py','-v')

$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')
$base = "cd $repoLinux && nextflow run main.nf -profile test --assay WGBS --wgbs_input_route"
$runs = @(
    "$base raw --input tests/fixtures/wgbs/raw_samplesheet.tsv --fasta tests/fixtures/wgbs/reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_raw --run_dm false -resume",
    "$base processed --input tests/fixtures/wgbs/processed_samplesheet.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_processed --run_dm true --contrast condition,treatment,control -resume",
    "$base processed --input tests/fixtures/wgbs/one_vs_one.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_1v1 --run_dm true --contrast condition,treatment,control -resume",
    "$base processed --input tests/fixtures/wgbs/processed_samplesheet.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_no_contrast --run_dm true -resume",
    "$base processed --input tests/fixtures/wgbs/rank_deficient.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_rank --run_dm true --design_formula '~ batch + condition' --contrast condition,treatment,control -resume",
    "$base processed --input tests/fixtures/wgbs/percentage_only.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_percentage --run_dm true --contrast condition,treatment,control -resume",
    "$base processed --input tests/fixtures/wgbs/low_coverage.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_low_coverage --run_dm true --contrast condition,treatment,control -resume",
    "$base processed --input tests/fixtures/wgbs/zero_coverage.tsv --fasta tests/fixtures/wgbs/processed_reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_guard_zero_coverage --run_dm true --contrast condition,treatment,control -resume",
    "$base raw --input tests/fixtures/wgbs/zero_align.tsv --fasta tests/fixtures/wgbs/reference.fa --gtf tests/fixtures/wgbs/reference.gtf --outdir tests/output/wgbs_zero_align --run_dm false -resume"
)
foreach ($run in $runs) { Invoke-Checked wsl.exe @('bash','-lc',$run) }
Invoke-Checked python @('tests/assert_wgbs_outputs.py')
