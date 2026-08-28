"""
Artifact URI parsing and formatting utilities.
Format: <scheme>://<study_id>/<artifact_name>/<version>
Examples:
  adata://AD_mouse_001/raw/v1
  adata://AD_mouse_001/microglia_subset/v4
  table://AD_mouse_001/pseudobulk_deg/v1
  fig://AD_mouse_001/paga_trajectory/v1
"""

import re
from typing import Tuple


URI_REGEX = re.compile(r"^([a-zA-Z0-9_\-]+)://([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-\./]+)/([a-zA-Z0-9_\-]+)$")


class ArtifactURI:
    def __init__(self, scheme: str, study_id: str, name: str, version: str):
        self.scheme = scheme.lower()
        self.study_id = study_id
        self.name = name.strip("/")
        self.version = version

    @classmethod
    def parse(cls, uri_str: str) -> "ArtifactURI":
        match = URI_REGEX.match(uri_str.strip())
        if not match:
            raise ValueError(f"Invalid artifact URI format: '{uri_str}'. Expected: scheme://study_id/artifact_name/version")
        scheme, study_id, name, version = match.groups()
        return cls(scheme=scheme, study_id=study_id, name=name, version=version)

    def to_string(self) -> str:
        return f"{self.scheme}://{self.study_id}/{self.name}/{self.version}"

    def __str__(self) -> str:
        return self.to_string()

    def __repr__(self) -> str:
        return f"ArtifactURI({self.to_string()})"

    def next_version(self) -> "ArtifactURI":
        """Increment integer version (e.g. v1 -> v2) or branch if alphabetic (e.g. v4a -> v4b)."""
        if self.version.startswith("v") and self.version[1:].isdigit():
            curr_num = int(self.version[1:])
            new_ver = f"v{curr_num + 1}"
        else:
            new_ver = f"{self.version}_next"
        return ArtifactURI(self.scheme, self.study_id, self.name, new_ver)

    def branch(self, branch_suffix: str) -> "ArtifactURI":
        """Create a branched version, e.g., v3 -> v4_harmony or v4a."""
        return ArtifactURI(self.scheme, self.study_id, self.name, f"{self.version}_{branch_suffix}")
