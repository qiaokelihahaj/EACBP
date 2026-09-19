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


URI_REGEX = re.compile(
    r"^([a-zA-Z0-9_\-]+)://([a-zA-Z0-9_\-]+)/([a-zA-Z0-9_\-./]+)/([a-zA-Z0-9_\-]+)$"
)
_COMPONENT_REGEX = re.compile(r"^[A-Za-z0-9_-]+$")
_NAME_COMPONENT_REGEX = re.compile(r"^[A-Za-z0-9_.-]+$")


class ArtifactURI:
    def __init__(self, scheme: str, study_id: str, name: str, version: str):
        self.scheme = scheme.lower()
        self.study_id = study_id
        self.name = name
        self.version = version
        self._validate_components()

    def _validate_components(self) -> None:
        """Validate URI components before they can become filesystem paths.

        Artifact names may contain slash-separated namespaces, but path traversal
        components, backslashes, empty segments, and platform-specific absolute
        path forms are rejected.  Keeping the validation here means callers that
        only parse a URI receive the same safety guarantees as the storage layer.
        """
        if not _COMPONENT_REGEX.fullmatch(self.scheme):
            raise ValueError(f"Invalid artifact URI scheme: {self.scheme!r}")
        if not _COMPONENT_REGEX.fullmatch(self.study_id):
            raise ValueError(f"Invalid artifact study_id: {self.study_id!r}")
        if "\\" in self.name or self.name.startswith("/") or self.name.endswith("/"):
            raise ValueError(f"Invalid artifact name path: {self.name!r}")
        name_parts = self.name.split("/")
        if not name_parts or any(
            part in {"", ".", ".."} or not _NAME_COMPONENT_REGEX.fullmatch(part)
            for part in name_parts
        ):
            raise ValueError(
                f"Invalid artifact name path: {self.name!r}; path traversal is not allowed"
            )
        if not _COMPONENT_REGEX.fullmatch(self.version):
            raise ValueError(f"Invalid artifact version: {self.version!r}")

    @classmethod
    def parse(cls, uri_str: str) -> "ArtifactURI":
        if not isinstance(uri_str, str):
            raise TypeError("Artifact URI must be a string")
        cleaned = uri_str.strip()
        match = URI_REGEX.fullmatch(cleaned)
        if not match:
            raise ValueError(
                f"Invalid artifact URI format: '{uri_str}'. "
                "Expected: scheme://study_id/artifact_name/version"
            )
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
