"""Resume and checkpoint helpers for the scientific execution loop."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.orchestrator.checkpoint import fingerprint
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


@dataclass
class ResumeLookup:
    """A journal/registry lookup result."""

    result: Optional[TaskResult]
    saved: Optional[Dict[str, Any]]
    reused: bool = False


class ResumeManager:
    """Build signatures, validate saved outputs, and create journal entries."""

    def __init__(self, artifact_registry: ArtifactRegistry):
        self.artifact_registry = artifact_registry

    @staticmethod
    def signature(
        task: TaskContract,
        input_hashes: Mapping[str, str],
        manifest: StudyManifest,
        environment: Mapping[str, Any],
    ) -> str:
        return fingerprint({
            "task": task.model_dump(),
            "inputs": dict(input_hashes),
            "manifest": manifest.model_dump(),
            "environment": dict(environment),
        })

    def input_hashes(self, task: TaskContract) -> Dict[str, str]:
        return {
            uri: self.artifact_registry.get(uri)[0].sha256_hash
            for uri in task.input_artifacts
        }

    def lookup(self, journal: Any, task: TaskContract, signature: str) -> ResumeLookup:
        """Return a validated saved success, if one can be safely reused.

        Failed journal entries are superseded by a committed registry receipt
        when available.  This preserves the crash-recovery behavior of the
        original loop while rejecting signatures from a different run.
        """

        saved = deepcopy(journal.entries.get(task.task_id))
        receipt = self.artifact_registry.get_task_commit(signature)
        if (
            not saved
            or (
                saved.get("signature") == signature
                and saved.get("result", {}).get("status")
                in (TaskStatus.EXECUTION_FAILURE.value, TaskStatus.METHOD_FAILURE.value)
            )
        ):
            saved = receipt or saved
        if saved and saved.get("signature") != signature:
            raise ValueError("Resume configuration/input differs from the saved task; use a new run.")
        if not saved or saved.get("result", {}).get("status") != TaskStatus.SUCCESS.value:
            return ResumeLookup(result=None, saved=saved, reused=False)

        for uri, digest in saved.get("output_hashes", {}).items():
            meta, _ = self.artifact_registry.get(uri)
            if meta.sha256_hash != digest:
                raise ValueError(f"Resume output metadata changed: {uri}")
        return ResumeLookup(
            result=TaskResult.model_validate(saved["result"]),
            saved=saved,
            reused=True,
        )

    def output_hashes(self, result: TaskResult) -> Dict[str, str]:
        return {
            uri: self.artifact_registry.get_metadata(uri).sha256_hash
            for uri in result.output_artifacts
            if self.artifact_registry.exists(uri)
        }

    def journal_entry(
        self,
        *,
        signature: str,
        result: TaskResult,
        phase: str,
        environment: Mapping[str, Any],
    ) -> Dict[str, Any]:
        return {
            "signature": signature,
            "result": result.model_dump(),
            "phase": phase,
            "environment": deepcopy(dict(environment)),
            "output_hashes": self.output_hashes(result),
        }
