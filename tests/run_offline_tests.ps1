$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Command,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Command exited with status $LASTEXITCODE"
    }
}

function Remove-TestArtifacts {
    $repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
    $paths = @(
        'tests/fixtures/generated',
        'tests/fixtures/atac/raw',
        'tests/fixtures/chip/raw',
        'tests/fixtures/wgbs/raw',
        'tests/output',
        'tests/work',
        'tests/reporting',
        '.nextflow',
        'work',
        'results',
        'output',
        'outputs',
        'raw_inputs.present',
        'processed_inputs.present'
    )
    foreach ($relative in $paths) {
        $path = Join-Path $repo $relative
        if (Test-Path -LiteralPath $path) {
            Remove-Item -LiteralPath $path -Recurse -Force
        }
    }
    Get-ChildItem -LiteralPath $repo -Filter '.nextflow.log*' -Force -File |
        Remove-Item -Force
    Get-ChildItem -LiteralPath $repo -Filter '__pycache__' -Force -Directory -Recurse |
        Sort-Object FullName -Descending |
        ForEach-Object { Remove-Item -LiteralPath $_.FullName -Recurse -Force }
}

try {
    Remove-TestArtifacts
    Invoke-Checked python @('tests/create_synthetic_fixtures.py')
    Invoke-Checked python @('-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_*.py', '-v')

    $repoWindows = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
    $repoLinux = '/mnt/' + $repoWindows.Substring(0,1).ToLower() + $repoWindows.Substring(2).Replace('\','/')
    $runs = @(
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/raw_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/raw -resume",
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/processed --run_de true --contrast condition,treatment,control -resume",
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/one_vs_one.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_1v1 --run_de true --contrast condition,treatment,control -resume",
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_no_contrast --run_de true -resume",
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_rank --run_de true --design_formula '~ batch + condition' --contrast condition,treatment,control -resume",
        "cd $repoLinux && nextflow run main.nf -profile test --input tests/fixtures/generated/technical_replicates.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/technical_replicates --run_de true --contrast condition,treatment,control -resume"
    )
    foreach ($run in $runs) {
        Invoke-Checked wsl.exe @('bash', '-lc', $run)
    }

    Invoke-Checked python @('tests/assert_outputs.py')
}
finally {
    Remove-TestArtifacts
}
