"""
Knowledge Engine module orchestrating multi-source retrieval in Discovery Mode and Prior-Guided Mode.
Enforces epistemic tagging, generates structured EvidenceNodes, and synthesizes traceable knowledge reports.
"""

from typing import List, Dict, Any, Optional
import re
from datetime import datetime, timezone
from pydantic import BaseModel, Field

from eacbp.schemas.study import StudyManifest
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    ClaimNode,
    ClaimType,
    LanguageTier,
    ConfidenceScore,
)
from eacbp.knowledge.literature import LiteratureRetriever, LiteratureEvidence
from eacbp.knowledge.biological_db import BiologicalDBRetriever, GOEnrichment, PathwayEnrichment, GeneInfo


class KnowledgeReport(BaseModel):
    """Structured report produced by the Knowledge Engine."""
    mode: str = Field(..., description="'discovery' vs 'prior_guided'")
    prior_guided: bool = Field(False, description="Whether execution was guided by prior hypotheses")
    hypothesis_tested: Optional[str] = Field(None, description="Hypothesis statement tested in prior-guided mode")
    target_genes: List[str] = Field(default_factory=list, description="Target genes analyzed")
    literature_evidence: List[LiteratureEvidence] = Field(default_factory=list)
    go_enrichments: List[GOEnrichment] = Field(default_factory=list)
    pathway_enrichments: List[PathwayEnrichment] = Field(default_factory=list)
    gene_annotations: Dict[str, GeneInfo] = Field(default_factory=dict)
    evidence_nodes: List[EvidenceNode] = Field(default_factory=list)
    epistemic_tags: List[str] = Field(default_factory=list)
    summary: str = Field(..., description="High-level narrative summary of knowledge retrieval findings")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class KnowledgeEngine:
    """
    Central coordinator for the Knowledge Plane.
    Orchestrates multi-source literature, Gene Ontology ORA, Reactome pathways, and NCBI Gene annotations
    under Discovery Mode (unbiased) and Prior-Guided Mode (mandatory epistemic tagging).
    """

    def __init__(
        self,
        literature_retriever: Optional[LiteratureRetriever] = None,
        biological_db_retriever: Optional[BiologicalDBRetriever] = None,
    ):
        self.literature_retriever = literature_retriever or LiteratureRetriever()
        self.biological_db_retriever = biological_db_retriever or BiologicalDBRetriever()

    def execute_discovery(
        self,
        manifest: StudyManifest,
        deg_genes: List[str],
        top_n_genes: int = 20,
        max_literature: int = 5,
        fdr_cutoff: float = 0.05,
    ) -> KnowledgeReport:
        """
        Executes Discovery Mode knowledge retrieval.
        Data-driven analysis across top DEGs without biasing feature selection by prior hypotheses.

        Args:
            manifest: StudyManifest defining species, disease, and cell types.
            deg_genes: List of empirically identified differentially expressed genes.
            top_n_genes: Number of top DEGs to evaluate.
            max_literature: Max literature articles to retrieve.
            fdr_cutoff: FDR threshold for GO and pathway enrichments.

        Returns:
            KnowledgeReport with mode='discovery', prior_guided=False.
        """
        target_genes = [g.strip() for g in deg_genes[:top_n_genes] if g.strip()]
        species = manifest.biological_design.species
        disease = manifest.biological_design.disease
        cell_types = manifest.biological_design.target_cell_types

        # 1. Biological Database ORA (GO & Reactome)
        go_results = self.biological_db_retriever.query_go(
            gene_list=target_genes,
            fdr_cutoff=fdr_cutoff,
        )
        pathway_results = self.biological_db_retriever.query_pathway(
            gene_list=target_genes,
            fdr_cutoff=fdr_cutoff,
        )

        # 2. Gene Entity Lookups
        gene_annotations = {}
        for g in target_genes[:10]:
            info = self.biological_db_retriever.query_gene(g, species=species)
            if info:
                gene_annotations[g] = info

        # 3. Data-driven literature retrieval
        search_terms = target_genes[:5]
        if disease:
            search_terms.append(disease)
        if cell_types:
            search_terms.extend(cell_types)

        lit_results = self.literature_retriever.query_literature(
            query_terms=search_terms,
            max_results=max_literature,
        )

        # 4. Extract structured EvidenceNodes
        evidence_nodes: List[EvidenceNode] = []
        
        # Pathway and GO evidence nodes
        go_nodes = self.biological_db_retriever.to_evidence_nodes(
            enrichments=go_results[:3],
            task_id="task_knowledge_discovery_go",
            prior_guided=False,
        )
        evidence_nodes.extend(go_nodes)

        pw_nodes = self.biological_db_retriever.to_evidence_nodes(
            enrichments=pathway_results[:3],
            task_id="task_knowledge_discovery_reactome",
            prior_guided=False,
        )
        evidence_nodes.extend(pw_nodes)

        # Literature evidence nodes
        for lit in lit_results[:3]:
            evidence_nodes.append(
                self.literature_retriever.to_evidence_node(
                    evidence=lit,
                    task_id="task_knowledge_discovery_literature",
                    prior_guided=False,
                )
            )

        # 5. Narrative summary
        top_go = [f"{g.term} (FDR={g.fdr_q_value:.2e})" for g in go_results[:3]]
        top_pw = [f"{p.pathway_name} (FDR={p.fdr_q_value:.2e})" for p in pathway_results[:3]]
        top_lit = [f"{l.title} ({l.journal}, {l.year})" for l in lit_results[:2]]

        summary_lines = [
            f"Discovery Mode Knowledge Analysis for study '{manifest.study_id}' (unbiased data-driven).",
            f"Evaluated {len(target_genes)} top empirical DEGs: {', '.join(target_genes[:8])}.",
            f"Identified {len(go_results)} enriched GO terms (Top: {'; '.join(top_go) if top_go else 'None'}).",
            f"Identified {len(pathway_results)} enriched Reactome pathways (Top: {'; '.join(top_pw) if top_pw else 'None'}).",
            f"Retrieved {len(lit_results)} supporting publications (Top: {'; '.join(top_lit) if top_lit else 'None'}).",
        ]
        summary = "\n".join(summary_lines)

        return KnowledgeReport(
            mode="discovery",
            prior_guided=False,
            hypothesis_tested=None,
            target_genes=target_genes,
            literature_evidence=lit_results,
            go_enrichments=go_results,
            pathway_enrichments=pathway_results,
            gene_annotations=gene_annotations,
            evidence_nodes=evidence_nodes,
            epistemic_tags=["mode:unbiased_discovery"],
            summary=summary,
        )

    def execute_prior_guided(
        self,
        manifest: StudyManifest,
        hypothesis: str,
        target_genes: Optional[List[str]] = None,
        deg_genes: Optional[List[str]] = None,
        max_literature: int = 5,
        fdr_cutoff: float = 0.05,
    ) -> KnowledgeReport:
        """
        Executes Prior-Guided Mode knowledge retrieval.
        Evaluates a specific user-defined biological hypothesis (e.g. DAM TREM2-APOE activation axis)
        with mandatory epistemic tagging across evidence nodes and claim metadata.

        Args:
            manifest: StudyManifest.
            hypothesis: User-defined hypothesis string.
            target_genes: List of hypothesized genes (e.g. ['Trem2', 'Apoe', 'Clec7a']).
            deg_genes: Optional empirical DEGs to cross-validate against hypothesis.
            max_literature: Max literature articles to retrieve.
            fdr_cutoff: FDR threshold for GO and pathway enrichments.

        Returns:
            KnowledgeReport with mode='prior_guided', prior_guided=True, and epistemic tagging.
        """
        # 1. Resolve target genes from input or hypothesis text
        if not target_genes:
            extracted = []
            tokens = re.split(r'[\s,\-_\(\)/]+', hypothesis)
            for t in tokens:
                clean_t = t.strip().upper()
                if clean_t and (clean_t in self.biological_db_retriever._genes or self.biological_db_retriever.query_gene(clean_t)):
                    extracted.append(clean_t)
            target_genes = extracted if extracted else ["TREM2", "APOE", "CLEC7A", "TYROBP"]

        species = manifest.biological_design.species
        disease = manifest.biological_design.disease or "Alzheimer's Disease"
        cell_types = manifest.biological_design.target_cell_types

        # 2. Targeted Literature Retrieval
        search_terms = [hypothesis] + target_genes
        if disease:
            search_terms.append(disease)
        if cell_types:
            search_terms.extend(cell_types)

        lit_results = self.literature_retriever.query_literature(
            query_terms=search_terms,
            max_results=max_literature,
        )

        # 3. Targeted Biological DB ORA on Prior Gene Set
        eval_genes = list(target_genes)
        if deg_genes:
            # Add intersecting empirical DEGs
            deg_upper = [g.upper() for g in deg_genes]
            intersecting = [g for g in eval_genes if g.upper() in deg_upper]
            if intersecting:
                eval_genes = list(set(eval_genes + intersecting))

        go_results = self.biological_db_retriever.query_go(
            gene_list=eval_genes,
            fdr_cutoff=fdr_cutoff,
        )
        pathway_results = self.biological_db_retriever.query_pathway(
            gene_list=eval_genes,
            fdr_cutoff=fdr_cutoff,
        )

        # 4. Gene Entity Lookups
        gene_annotations = {}
        for g in target_genes:
            info = self.biological_db_retriever.query_gene(g, species=species)
            if info:
                gene_annotations[g] = info

        # 5. Extract structured EvidenceNodes with MANDATORY EPISTEMIC TAGS
        epistemic_tag = "[PRIOR-GUIDED HYPOTHESIS TESTING]"
        evidence_nodes: List[EvidenceNode] = []

        go_nodes = self.biological_db_retriever.to_evidence_nodes(
            enrichments=go_results[:3],
            task_id="task_knowledge_prior_go",
            prior_guided=True,
            hypothesis=hypothesis,
        )
        evidence_nodes.extend(go_nodes)

        pw_nodes = self.biological_db_retriever.to_evidence_nodes(
            enrichments=pathway_results[:3],
            task_id="task_knowledge_prior_reactome",
            prior_guided=True,
            hypothesis=hypothesis,
        )
        evidence_nodes.extend(pw_nodes)

        for lit in lit_results[:3]:
            evidence_nodes.append(
                self.literature_retriever.to_evidence_node(
                    evidence=lit,
                    task_id="task_knowledge_prior_literature",
                    prior_guided=True,
                    hypothesis=hypothesis,
                )
            )

        # 6. Narrative summary with mandatory epistemic disclosure
        top_go = [f"{g.term} (FDR={g.fdr_q_value:.2e})" for g in go_results[:3]]
        top_pw = [f"{p.pathway_name} (FDR={p.fdr_q_value:.2e})" for p in pathway_results[:3]]
        top_lit = [f"{l.title} ({l.journal}, {l.year})" for l in lit_results[:2]]

        summary_lines = [
            f"{epistemic_tag}: Evaluated targeted prior hypothesis '{hypothesis}'.",
            f"Target Hypothesized Gene Axis: {', '.join(target_genes)}.",
            f"Prior Pathway Support: {len(pathway_results)} pathways (Top: {'; '.join(top_pw) if top_pw else 'None'}).",
            f"Prior Gene Ontology Support: {len(go_results)} GO terms (Top: {'; '.join(top_go) if top_go else 'None'}).",
            f"Prior Literature Citations: {len(lit_results)} articles (Top: {'; '.join(top_lit) if top_lit else 'None'}).",
            f"Epistemic Status: This analysis reflects confirmatory hypothesis evaluation rather than unbiased discovery.",
        ]
        summary = "\n".join(summary_lines)

        return KnowledgeReport(
            mode="prior_guided",
            prior_guided=True,
            hypothesis_tested=hypothesis,
            target_genes=target_genes,
            literature_evidence=lit_results,
            go_enrichments=go_results,
            pathway_enrichments=pathway_results,
            gene_annotations=gene_annotations,
            evidence_nodes=evidence_nodes,
            epistemic_tags=[epistemic_tag, f"{epistemic_tag}: {hypothesis.upper()}"],
            summary=summary,
        )

    def synthesize_claim(
        self,
        report: KnowledgeReport,
        claim_id: str = "C_knowledge_001",
    ) -> ClaimNode:
        """
        Synthesizes a traceable ClaimNode backed by the KnowledgeReport evidence nodes.
        Maintains strict 4-tier language calibration and epistemic provenance.
        """
        eids = [node.evidence_id for node in report.evidence_nodes]
        
        # Calculate confidence from evidence nodes
        scores = [n.score for n in report.evidence_nodes]
        mech_score = float(sum(scores) / len(scores)) if scores else 0.5
        overall_conf = round(0.45 * mech_score + 0.35 * 0.8, 3)

        if report.prior_guided:
            tag = "[PRIOR-GUIDED HYPOTHESIS TESTING]"
            statement = f"{tag}: Prior-guided evaluation supports the biological plausibility of the hypothesized axis '{report.hypothesis_tested}'."
            lang_tier = LanguageTier.LEVEL_4_HYPOTHESIS
            prov_summary = f"{tag} Claim {claim_id} evaluates user prior hypothesis '{report.hypothesis_tested}' supported by {len(report.evidence_nodes)} knowledge evidence items."
        else:
            statement = f"Data-driven pathway enrichment and literature context support activation of {report.target_genes[:3]} in {report.mode} analysis."
            lang_tier = LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION
            prov_summary = f"Claim {claim_id} supported by {len(report.evidence_nodes)} discovery knowledge evidence items."

        return ClaimNode(
            claim_id=claim_id,
            statement=statement,
            language_tier=lang_tier,
            claim_type=ClaimType.MECHANISTIC_HYPOTHESIS if report.prior_guided else ClaimType.REGULATORY,
            causal_status="observational",
            support_evidence_ids=eids,
            confidence=ConfidenceScore(
                association=0.7,
                mechanistic=round(mech_score, 3),
                causal=0.0,
                overall=overall_conf,
            ),
            provenance_summary=prov_summary,
            created_at=datetime.now(timezone.utc),
        )
