"""Local curated knowledge retrieval exposed through the capability contract."""
import pandas as pd
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.knowledge.engine import KnowledgeEngine


class KnowledgeRetrievalCapability(BaseCapability):
    """
    Capability wrapping Multi-Source Knowledge Engine for Discovery and Prior-Guided knowledge retrieval.
    """

    def __init__(self, implementation_id: str = "knowledge_engine_discovery_v1"):
        super().__init__(
            capability_name="knowledge_retrieval",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["scRNA", "spatial", "genomics"],
            accepts_types=[ArtifactType.TABLE, ArtifactType.ANNDATA, ArtifactType.GENE_LIST],
            output_types=[ArtifactType.TABLE, ArtifactType.JSON],
            suitable_for=["literature_retrieval", "pathway_enrichment", "prior_guided_hypothesis_testing"],
        )
        self.engine = KnowledgeEngine()

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0] if contract.input_artifacts else ""
        deg_genes = []

        if in_uri and registry.exists(in_uri):
            meta, payload = registry.get(in_uri)
            if meta.type == ArtifactType.TABLE:
                df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
                if "gene" in df.columns:
                    if "fdr_q_value" in df:
                        df = df[pd.to_numeric(df["fdr_q_value"], errors="coerce").between(0, .05, inclusive="left")]
                    else:
                        df = df.iloc[:0]
                    deg_genes = df["gene"].dropna().head(20).tolist()
            elif meta.type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
                data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
                if "gene_name" in data.var.columns:
                    deg_genes = data.var["gene_name"].dropna().head(20).tolist()

        is_prior = bool(contract.parameters.get("prior_guided", False)) or ("prior" in (contract.method or ""))
        hypotheses = contract.parameters.get("hypotheses", [])
        hypothesis = hypotheses[0] if isinstance(hypotheses, list) and hypotheses else str(hypotheses)

        # Manifest representation for knowledge engine
        study_manifest = StudyManifest(
            study_id=contract.parameters.get("study_id", "knowledge_study"),
            biological_design=BiologicalDesign(
                species=contract.parameters.get("species", "unknown"),
                tissue=contract.parameters.get("tissue", "unknown"),
                disease=contract.parameters.get("disease"),
                target_cell_types=contract.parameters.get("cell_types", []),
            ),
        )

        if is_prior:
            report = self.engine.execute_prior_guided(
                manifest=study_manifest,
                hypothesis=hypothesis,
                target_genes=contract.parameters.get("target_genes") or deg_genes[:5],
            )
        else:
            report = self.engine.execute_discovery(
                manifest=study_manifest,
                deg_genes=deg_genes,
                top_n_genes=20,
            )

        # Build table payload
        table_rows = []
        for go in report.go_enrichments:
            table_rows.append({
                "category": "GO_Biological_Process",
                "id": go.go_id,
                "name": go.term,
                "p_value": go.p_value,
                "fdr_q_value": go.fdr_q_value,
                "fold_enrichment": go.fold_enrichment,
                "genes": ", ".join(go.genes) if hasattr(go, "genes") else "",
            })
        for pw in report.pathway_enrichments:
            table_rows.append({
                "category": "Reactome_Pathway",
                "id": pw.pathway_id,
                "name": pw.pathway_name,
                "p_value": pw.p_value,
                "fdr_q_value": pw.fdr_q_value,
                "fold_enrichment": pw.fold_enrichment,
                "genes": ", ".join(pw.genes) if hasattr(pw, "genes") else "",
            })
        for lit in report.literature_evidence:
            table_rows.append({
                "category": "Local_Curated_Literature",
                "id": lit.pmid or lit.doi or "unidentified",
                "name": lit.title,
                "p_value": None,
                "fdr_q_value": None,
                "fold_enrichment": None,
                "genes": ", ".join(lit.matched_keywords) if hasattr(lit, "matched_keywords") else "",
            })

        evidence_df = pd.DataFrame(table_rows) if table_rows else pd.DataFrame([{"category": "No_evidence", "name": "No evidence"}])

        sid = contract.parameters.get("study_id", "study")
        # Target branches provide explicit expected outputs; retain the
        # historical study-level URIs for single-target and broad studies.
        table_uri = (
            contract.expected_outputs[0]
            if contract.expected_outputs
            else f"table://{sid}/knowledge_evidence/v1"
        )
        json_uri = (
            contract.expected_outputs[1]
            if len(contract.expected_outputs) > 1
            else f"json://{sid}/knowledge_report/v1"
        )

        # Register versioned artifacts
        registry.register(
            uri_str=table_uri,
            payload=evidence_df,
            artifact_type=ArtifactType.TABLE,
            study_id=sid,
            created_by_task=contract.task_id,
            operation="knowledge_enrichment_table",
            parent_uris=[in_uri] if in_uri else [],
            summary_metrics={
                "n_enrichments": len(evidence_df),
                "prior_guided": report.prior_guided,
                "target_cell_type": contract.parameters.get("target_cell_type"),
                "target_branch": contract.parameters.get("target_branch"),
            },
        )

        registry.register(
            uri_str=json_uri,
            payload=report.model_dump(mode="json"),
            artifact_type=ArtifactType.JSON,
            study_id=sid,
            created_by_task=contract.task_id,
            operation="knowledge_report_json",
            parent_uris=[in_uri] if in_uri else [],
            summary_metrics={
                "mode": report.mode,
                "prior_guided": report.prior_guided,
                "source_mode": "local_curated_unverified",
                "target_cell_type": contract.parameters.get("target_cell_type"),
                "target_branch": contract.parameters.get("target_branch"),
            },
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=contract.method or self.implementation_id,
            output_artifacts=[table_uri, json_uri],
            metrics={
                "report": report.model_dump(mode="json"),
                "mode": report.mode,
                "prior_guided": report.prior_guided,
                "evidence_nodes": [e.model_dump(mode="json") for e in report.evidence_nodes],
                "target_genes": report.target_genes,
                "summary": report.summary,
            },
        )
