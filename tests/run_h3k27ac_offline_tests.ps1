$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param([Parameter(Mandatory=$true)][string]$Command,[Parameter(Mandatory=$true)][string[]]$Arguments)
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Command exited with status $LASTEXITCODE" }
}

Invoke-Checked python @('tests/create_chip_fixtures.py')
Invoke-Checked python @('-m','unittest','tests/test_chip_contracts.py','-v')

$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')
$rawBase = "cd $repoLinux && nextflow run main.nf -profile test --assay H3K27ac --chip_input_route raw --fasta tests/fixtures/chip/reference.fa --gtf tests/fixtures/chip/reference.gtf --blacklist tests/fixtures/chip/blacklist.bed --reference_scope synthetic --chip_effective_genome_size 2000000 --run_chip_diff false"
$processedBase = "cd $repoLinux && nextflow run main.nf -profile test --assay H3K27ac --chip_input_route processed --fasta tests/fixtures/chip/reference.fa --gtf tests/fixtures/chip/reference.gtf --blacklist tests/fixtures/chip/blacklist.bed --reference_scope synthetic --run_chip_diff true"
$runs = @(
    "$rawBase --input tests/fixtures/chip/raw_single.tsv --outdir tests/output/chip_raw_single -resume",
    "$rawBase --input tests/fixtures/chip/raw_paired.tsv --outdir tests/output/chip_raw_paired -resume",
    "$processedBase --input tests/fixtures/chip/processed.tsv --outdir tests/output/chip_processed --contrast condition,treatment,control -resume",
    "$processedBase --input tests/fixtures/chip/one_vs_one.tsv --outdir tests/output/chip_guard_1v1 --contrast condition,treatment,control -resume",
    "$processedBase --input tests/fixtures/chip/processed.tsv --outdir tests/output/chip_guard_no_contrast -resume",
    "$processedBase --input tests/fixtures/chip/rank_deficient.tsv --outdir tests/output/chip_guard_rank --design_formula '~ batch + condition' --contrast condition,treatment,control -resume"
)
foreach ($run in $runs) { Invoke-Checked wsl.exe @('bash','-lc',$run) }
Invoke-Checked python @('tests/assert_chip_outputs.py')
