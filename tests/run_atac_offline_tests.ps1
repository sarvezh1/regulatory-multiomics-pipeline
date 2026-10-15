$ErrorActionPreference = 'Stop'
function Invoke-Checked { param([string]$Command,[string[]]$Arguments); & $Command @Arguments; if($LASTEXITCODE -ne 0){throw "$Command exited with status $LASTEXITCODE"} }
Invoke-Checked python @('tests/create_atac_fixtures.py')
Invoke-Checked python @('-m','unittest','discover','-s','tests','-p','test_*.py','-v')
$repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')
$base="cd $repoLinux && nextflow run main.nf -profile test --assay ATAC --atac_input_route processed --fasta tests/fixtures/atac/reference.fa --gtf tests/fixtures/atac/reference.gtf --blacklist tests/fixtures/atac/blacklist.bed"
$runs=@(
 "$base --input tests/fixtures/atac/processed_samplesheet.tsv --outdir tests/output/atac_processed --run_da true --contrast condition,treatment,control -resume",
 "$base --input tests/fixtures/atac/one_vs_one.tsv --outdir tests/output/atac_guard_1v1 --run_da true --contrast condition,treatment,control -resume",
 "$base --input tests/fixtures/atac/processed_samplesheet.tsv --outdir tests/output/atac_guard_no_contrast --run_da true -resume",
 "$base --input tests/fixtures/atac/rank_deficient.tsv --outdir tests/output/atac_guard_rank --run_da true --design_formula '~ batch + condition' --contrast condition,treatment,control -resume",
 "cd $repoLinux && nextflow run main.nf -profile test --assay ATAC --atac_input_route raw --input tests/fixtures/atac/raw_samplesheet.tsv --fasta tests/fixtures/atac/reference.fa --gtf tests/fixtures/atac/reference.gtf --blacklist tests/fixtures/atac/blacklist.bed --outdir tests/output/atac_raw --reference_scope synthetic --run_da false -resume"
)
foreach($run in $runs){Invoke-Checked wsl.exe @('bash','-lc',$run)}
Invoke-Checked python @('tests/assert_atac_outputs.py')
