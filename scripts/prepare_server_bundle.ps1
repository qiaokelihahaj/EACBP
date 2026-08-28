param(
    [string]$BundleName = "server_bundle"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$DistRoot = Join-Path $ProjectRoot "dist"
$BundleRoot = Join-Path $DistRoot $BundleName
$StagingRoot = Join-Path $DistRoot (".{0}.staging-{1}" -f $BundleName, [Guid]::NewGuid().ToString("N"))

$files = @(
    "README.md",
    "PROJECT.md",
    "ORIGINAL_REQUEST.md",
    "TEST_INFRA.md",
    "TEST_READY.md",
    "SERVER_MIGRATION.md",
    "pyproject.toml",
    "environment.yml",
    "scripts",
    "eacbp",
    "tests"
)

New-Item -ItemType Directory -Path $DistRoot -Force | Out-Null
New-Item -ItemType Directory -Path $StagingRoot | Out-Null

try {
    foreach ($relative in $files) {
        $source = Join-Path $ProjectRoot $relative
        if (-not (Test-Path -LiteralPath $source)) {
            continue
        }
        $destination = Join-Path $StagingRoot $relative
        if ((Get-Item -LiteralPath $source).PSIsContainer) {
            Get-ChildItem -LiteralPath $source -Recurse -File | Where-Object {
                $_.FullName -notmatch "\\__pycache__\\" -and
                $_.FullName -notmatch "\\.pytest_cache\\" -and
                $_.FullName -notmatch "\\.ipynb_checkpoints\\" -and
                $_.Name -notmatch "\.pyc$" -and
                $_.Length -gt 0
            } | ForEach-Object {
                $relativeFile = $_.FullName.Substring($source.Length).TrimStart("\")
                $targetFile = Join-Path $destination $relativeFile
                New-Item -ItemType Directory -Path (Split-Path -Parent $targetFile) -Force | Out-Null
                Copy-Item -LiteralPath $_.FullName -Destination $targetFile
            }
        } else {
            Copy-Item -LiteralPath $source -Destination $destination
        }
    }

    $manifestFiles = Get-ChildItem -LiteralPath $StagingRoot -File -Recurse | ForEach-Object {
        [ordered]@{
            path = $_.FullName.Substring($StagingRoot.Length + 1).Replace("\", "/")
            bytes = $_.Length
            sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName).Hash.ToLowerInvariant()
        }
    }
    $manifest = [ordered]@{
        schema_version = 1
        artifact_type = "eacbp-server-bundle"
        created_utc = [DateTime]::UtcNow.ToString("o")
        source_directory = Split-Path -Leaf $ProjectRoot
        files = $manifestFiles
    }
    $manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $StagingRoot "bundle_manifest.json") -Encoding utf8

    if (Test-Path -LiteralPath $BundleRoot) {
        Remove-Item -LiteralPath $BundleRoot -Recurse -Force
    }
    Move-Item -LiteralPath $StagingRoot -Destination $BundleRoot
    Write-Output "Created reproducible bundle: $BundleRoot"
} catch {
    if (Test-Path -LiteralPath $StagingRoot) {
        Remove-Item -LiteralPath $StagingRoot -Recurse -Force
    }
    throw
}
