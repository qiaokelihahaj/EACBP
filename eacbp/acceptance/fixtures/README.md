# Research acceptance fixtures

This directory contains only a fill-in manifest template. No real h5ad dataset was supplied for P1, and no fixture here should be read as biologically validated. Synthetic regression data is generated inside the dedicated acceptance tests and is marked `dataset_kind="synthetic"` in their result assertions.

Copy `research_manifest.template.json`, then replace the source, license, file digest, design fields, and exact reference environment values. The CLI rejects a real-dataset manifest without the pinned dependency list and a SHA256-verified environment lockfile.
