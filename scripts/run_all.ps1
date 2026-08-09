[CmdletBinding()]
param(
    [int]$Seed = 7,
    [int]$Replicates = 10000,
    [switch]$DoNotOpenOutput
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
$VenvPath = Join-Path $RepoRoot ".venv"
$VenvPython = Join-Path $VenvPath "Scripts\python.exe"
$OutputPath = Join-Path $RepoRoot "outputs"

if (-not (Test-Path -LiteralPath $VenvPython)) {
    python -m venv $VenvPath
}

& $VenvPython -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed" }

Push-Location $RepoRoot
try {
    & $VenvPython -m pip install -e ".[dev]"
    if ($LASTEXITCODE -ne 0) { throw "Package installation failed" }

    & $VenvPython -m compileall -q src scripts tests
    if ($LASTEXITCODE -ne 0) { throw "Compilation check failed" }

    & $VenvPython -m pytest
    if ($LASTEXITCODE -ne 0) { throw "Test suite failed" }

    & $VenvPython -m correlated_abiogenesis.appendix_b `
        --seed $Seed `
        --replicates $Replicates `
        --output-dir $OutputPath
    if ($LASTEXITCODE -ne 0) { throw "Appendix B reproduction failed" }
}
finally {
    Pop-Location
}

Write-Host ""
Write-Host "Repository validation complete"
Write-Host "Tests: PASS"
Write-Host "Master seed: $Seed"
Write-Host "Estimator replicates: $Replicates"
Write-Host "Outputs: $OutputPath"

if (-not $DoNotOpenOutput) {
    Invoke-Item -LiteralPath $OutputPath
}
