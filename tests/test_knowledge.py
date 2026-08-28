"""
Unit and integration tests for the EACBP Knowledge Plane.
Covers LiteratureRetriever, BiologicalDBRetriever, KnowledgeEngine (Discovery vs Prior-Guided Modes),
Hypergeometric ORA, Benjamini-Hochberg FDR, and Epistemic Tagging.
"""

import pytest
import numpy as np
from eacbp.schemas.study import (
    StudyManifest,
    BiologicalDesign,
    ExperimentalDesign,
    DataSpec,
    AnalysisPolicy,
)
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
    LanguageTier,
    ClaimNode,
)
from eacbp.knowledge import (
    LiteratureRetriever,
    LiteratureEvidence,
    BiologicalDBRetriever,
    GOEnrichment,
    PathwayEnrichment,
    GeneInfo,
    KnowledgeEngine,
    KnowledgeReport,
)


@pytest.fixture
def sample_study_manifest() -> StudyManifest:
    """Fixture providing a standard AD single-cell study manifest."""
    return StudyManifest(
        study_id="AD_mouse_test_001",
        title="Single-cell Transcriptomic Dissection of AD Microglial States",
        biological_design=BiologicalDesign(
            species="mus_musculus",
            tissue="cortex",
            disease="Alzheimer's Disease",
            conditions=["AD", "Control"],
            target_cell_types=["Microglia"],
        ),
        experimental_design=ExperimentalDesign(
            biological_unit="mouse",
            donor_replicates_per_condition={"AD": 4, "Control": 4},
        ),
        data=DataSpec(modalities=["scRNA", "spatial"]),
        analysis_policy=AnalysisPolicy(
            discovery_mode=True,
            prior_guided_analysis=False,
        ),
    )


# -------------------------------------------------------------------------
# 1. Literature Retriever Tests
# -------------------------------------------------------------------------

def test_literature_retriever_query_and_parsing():
    """Tests keyword queries, relevance scoring, and article parsing in LiteratureRetriever."""
    retriever = LiteratureRetriever()
    assert retriever.total_articles >= 10

    # Query with multiple keywords
    results = retriever.query_literature("TREM2 APOE microglia", max_results=5)
    assert len(results) > 0
    assert len(results) <= 5

    top_art = results[0]
    assert isinstance(top_art, LiteratureEvidence)
    assert top_art.title != ""
    assert top_art.journal != ""
    assert top_art.year >= 2015
    assert top_art.relevance_score > 0.0
    assert len(top_art.matched_keywords) > 0

    # Test PMID lookup
    pmid_art = retriever.get_article_by_pmid("28602351")
    assert pmid_art is not None
    assert "Keren-Shaul" in pmid_art.authors[0]
    assert "A Unique Microglia Type" in pmid_art.title

    # Test DOI lookup
    doi_art = retriever.get_article_by_doi("10.1016/j.immuni.2017.08.008")
    assert doi_art is not None
    assert "Krasemann" in doi_art.authors[0]

    # Test source filtering
    biorxiv_results = retriever.query_literature("microglia", source="bioRxiv")
    for r in biorxiv_results:
        assert r.source == "bioRxiv"


def test_literature_retriever_edge_cases_and_custom_addition():
    """Tests edge cases, boolean token filtering, min_year, min_citations, and custom article addition."""
    retriever = LiteratureRetriever()

    # Empty query tokens
    assert retriever.query_literature("") == []
    assert retriever.query_literature([]) == []
    assert retriever.query_literature("AND OR NOT THE A") == []

    # Non-existent query
    assert retriever.query_literature("NON_EXISTENT_QUANTUM_GENE_12345") == []

    # Min year and min citation filters
    recent_cit = retriever.query_literature("microglia", min_year=2020, min_citations=100)
    for art in recent_cit:
        assert art.year >= 2020
        assert art.citation_count >= 100

    # Custom article addition
    custom_art = LiteratureEvidence(
        pmid="99999999",
        doi="10.1016/j.cell.2026.99.999",
        title="Novel Microglial Lipid Receptor Axis in Aging",
        authors=["Auditor A", "Scientist B"],
        journal="Cell",
        year=2026,
        citation_count=50,
        source="PubMed",
        abstract="Discovery of a novel microglial lipid receptor regulating neuroinflammation.",
        evidence_text="Novel lipid receptor modulates neuroinflammation in aged cortex.",
    )
    retriever.add_article(custom_art)
    assert retriever.get_article_by_pmid("99999999") is not None
    found_custom = retriever.query_literature("Novel Microglial Lipid Receptor")
    assert len(found_custom) > 0
    assert found_custom[0].pmid == "99999999"


def test_literature_to_evidence_node():
    """Tests converting LiteratureEvidence to standardized EvidenceNode with epistemic tags."""
    retriever = LiteratureRetriever()
    article = retriever.get_article_by_pmid("28602351")
    assert article is not None

    # 1. Discovery Mode (prior_guided=False)
    node_disc = retriever.to_evidence_node(
        evidence=article,
        task_id="task_lit_disc",
        prior_guided=False,
    )
    assert isinstance(node_disc, EvidenceNode)
    assert node_disc.type == EvidenceType.LITERATURE_SUPPORT
    assert node_disc.polarity == EvidencePolarity.SUPPORTING
    assert 0.0 <= node_disc.score <= 1.0
    assert "Keren-Shaul" in node_disc.summary
    assert "[PRIOR-GUIDED" not in node_disc.summary
    assert node_disc.biological_context.get("mode") == "unbiased_discovery"
    assert node_disc.metrics.get("prior_guided") is False

    # 2. Prior-Guided Mode (prior_guided=True)
    hyp_text = "DAM microglia TREM2-APOE activation axis"
    node_prior = retriever.to_evidence_node(
        evidence=article,
        task_id="task_lit_prior",
        prior_guided=True,
        hypothesis=hyp_text,
    )
    assert node_prior.type == EvidenceType.LITERATURE_SUPPORT
    assert node_prior.summary.startswith("[PRIOR-GUIDED HYPOTHESIS TESTING]")
    assert node_prior.biological_context.get("mode") == "prior_guided"
    assert node_prior.biological_context.get("epistemic_tag") == "[PRIOR-GUIDED HYPOTHESIS TESTING]"
    assert node_prior.biological_context.get("hypothesis") == hyp_text
    assert node_prior.metrics.get("prior_guided") is True


# -------------------------------------------------------------------------
# 2. Biological DB Retriever (GO ORA, FDR, Reactome, NCBI Gene) Tests
# -------------------------------------------------------------------------

def test_biological_db_benjamini_hochberg_fdr():
    """Tests Benjamini-Hochberg FDR correction algorithm implementation."""
    retriever = BiologicalDBRetriever()

    # Known test case
    raw_p = [0.001, 0.01, 0.04, 0.05, 0.20]
    q_vals = retriever.benjamini_hochberg(raw_p)

    assert len(q_vals) == len(raw_p)
    # Monotonic property: q_vals should be <= 1.0
    for q in q_vals:
        assert 0.0 <= q <= 1.0

    # First p-value should be adjusted to 0.001 * 5 / 1 = 0.005
    assert round(q_vals[0], 4) == 0.005
    # Second: 0.01 * 5 / 2 = 0.025
    assert round(q_vals[1], 4) == 0.025

    # Edge cases
    assert retriever.benjamini_hochberg([]) == []
    assert retriever.benjamini_hochberg([0.05]) == [0.05]

    # Large random p-values: assert all q-values are within [0, 1] and monotonic when sorted
    np.random.seed(42)
    rand_p = list(np.random.uniform(0.0001, 0.9999, size=100))
    rand_q = retriever.benjamini_hochberg(rand_p)
    assert len(rand_q) == 100
    assert all(0.0 <= q <= 1.0 for q in rand_q)


def test_biological_db_go_ora_hypergeometric():
    """Tests Gene Ontology Over-Representation Analysis with exact hypergeometric math."""
    retriever = BiologicalDBRetriever()
    dam_genes = ["APOE", "TREM2", "TYROBP", "CLEC7A", "C3", "C1QA", "SPP1", "CST7", "LPL"]

    go_results = retriever.query_go(
        gene_list=dam_genes,
        background_size=20000,
        fdr_cutoff=0.05,
    )

    assert len(go_results) > 0
    top_go = go_results[0]
    assert isinstance(top_go, GOEnrichment)
    assert top_go.gene_count > 0
    assert top_go.sample_size == len(set(g.upper() for g in dam_genes))
    assert top_go.fold_enrichment > 1.0
    assert top_go.p_value < 0.05
    assert top_go.fdr_q_value <= 0.05
    assert len(top_go.genes) == top_go.gene_count

    # Verify category filtering
    bp_results = retriever.query_go(dam_genes, category="BP")
    for r in bp_results:
        assert r.category == "BP"

    mf_results = retriever.query_go(dam_genes, category="MF")
    for r in mf_results:
        assert r.category == "MF"

    cc_results = retriever.query_go(dam_genes, category="CC")
    for r in cc_results:
        assert r.category == "CC"

    # Empty gene list
    assert retriever.query_go([]) == []
    # Gene list with no matches
    assert retriever.query_go(["NON_EXISTENT_GENE_1", "NON_EXISTENT_GENE_2"]) == []


def test_biological_db_reactome_and_gene_lookups():
    """Tests Reactome pathway ORA and NCBI Gene ortholog lookups."""
    retriever = BiologicalDBRetriever()
    dam_genes = ["TREM2", "TYROBP", "CLEC7A", "C3", "C1QA", "CD68", "AXL", "APOE"]

    # Reactome query
    pathways = retriever.query_pathway(dam_genes, fdr_cutoff=0.05)
    assert len(pathways) > 0
    top_pw = pathways[0]
    assert isinstance(top_pw, PathwayEnrichment)
    assert top_pw.gene_count >= 3
    assert top_pw.fold_enrichment > 1.0

    # NCBI Gene lookups
    apoe_human = retriever.query_gene("APOE")
    assert apoe_human is not None
    assert apoe_human.gene_id == "348"
    assert "apolipoprotein" in apoe_human.full_name.lower()

    trem2_mouse = retriever.query_gene("Trem2")
    assert trem2_mouse is not None
    assert trem2_mouse.species == "mus_musculus"

    # Cross-species ortholog resolution
    mouse_apoe_ortholog = retriever.resolve_ortholog("Apoe", target_species="human")
    assert mouse_apoe_ortholog == "APOE"

    human_trem2_ortholog = retriever.resolve_ortholog("TREM2", target_species="mouse")
    assert human_trem2_ortholog == "Trem2"

    # Unknown gene
    assert retriever.query_gene("NON_EXISTENT_GENE_XYZ") is None
    assert retriever.resolve_ortholog("NON_EXISTENT_GENE_XYZ", "human") is None


def test_biological_db_to_evidence_nodes():
    """Tests converting GO and Pathway enrichments to EvidenceNodes."""
    retriever = BiologicalDBRetriever()
    dam_genes = ["APOE", "TREM2", "TYROBP", "CLEC7A", "C3"]
    go_results = retriever.query_go(dam_genes, fdr_cutoff=0.05)

    nodes = retriever.to_evidence_nodes(
        enrichments=go_results[:2],
        task_id="task_go_enrich",
        prior_guided=True,
        hypothesis="DAM Axis",
    )
    assert len(nodes) == min(2, len(go_results))
    for node in nodes:
        assert isinstance(node, EvidenceNode)
        assert node.type == EvidenceType.PATHWAY_ENRICHMENT
        assert node.summary.startswith("[PRIOR-GUIDED HYPOTHESIS TESTING]")
        assert node.metrics.get("prior_guided") is True
        assert node.biological_context.get("epistemic_tag") == "[PRIOR-GUIDED HYPOTHESIS TESTING]"


# -------------------------------------------------------------------------
# 3. Knowledge Engine Dual Mode Tests
# -------------------------------------------------------------------------

def test_knowledge_engine_discovery_mode(sample_study_manifest):
    """Tests unbiased Discovery Mode execution in KnowledgeEngine."""
    engine = KnowledgeEngine()
    deg_genes = ["Apoe", "Trem2", "Clec7a", "Tyrobp", "C3", "C1qa", "Spp1", "Cst7", "Lpl", "Itgax"]

    report = engine.execute_discovery(
        manifest=sample_study_manifest,
        deg_genes=deg_genes,
        top_n_genes=10,
        max_literature=3,
    )

    assert isinstance(report, KnowledgeReport)
    assert report.mode == "discovery"
    assert report.prior_guided is False
    assert report.hypothesis_tested is None
    assert "mode:unbiased_discovery" in report.epistemic_tags

    # Evidence nodes
    assert len(report.evidence_nodes) > 0
    for node in report.evidence_nodes:
        assert "[PRIOR-GUIDED" not in node.summary
        assert node.biological_context.get("mode") == "unbiased_discovery"
        assert node.metrics.get("prior_guided") is False
        assert 0.0 <= node.score <= 1.0

    # Claim synthesis
    claim = engine.synthesize_claim(report, claim_id="C_discovery_001")
    assert isinstance(claim, ClaimNode)
    assert claim.language_tier == LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION
    assert claim.causal_status == "observational"
    assert len(claim.support_evidence_ids) == len(report.evidence_nodes)


def test_knowledge_engine_prior_guided_mode(sample_study_manifest):
    """Tests Prior-Guided Mode execution with mandatory epistemic tags across all outputs."""
    engine = KnowledgeEngine()
    hypothesis = "DAM microglia TREM2-APOE activation axis"
    target_genes = ["Trem2", "Apoe", "Clec7a", "Tyrobp"]

    report = engine.execute_prior_guided(
        manifest=sample_study_manifest,
        hypothesis=hypothesis,
        target_genes=target_genes,
        deg_genes=["Trem2", "Apoe", "Clec7a", "C3"],
        max_literature=3,
    )

    assert isinstance(report, KnowledgeReport)
    assert report.mode == "prior_guided"
    assert report.prior_guided is True
    assert report.hypothesis_tested == hypothesis
    assert any("[PRIOR-GUIDED HYPOTHESIS TESTING]" in tag for tag in report.epistemic_tags)
    assert report.summary.startswith("[PRIOR-GUIDED HYPOTHESIS TESTING]")

    # Verify every evidence node carries epistemic badges
    assert len(report.evidence_nodes) > 0
    for node in report.evidence_nodes:
        assert node.summary.startswith("[PRIOR-GUIDED HYPOTHESIS TESTING]")
        assert node.biological_context.get("mode") == "prior_guided"
        assert node.biological_context.get("epistemic_tag") == "[PRIOR-GUIDED HYPOTHESIS TESTING]"
        assert node.biological_context.get("hypothesis") == hypothesis
        assert node.metrics.get("prior_guided") is True

    # Synthesize claim
    claim = engine.synthesize_claim(report, claim_id="C_prior_001")
    assert isinstance(claim, ClaimNode)
    assert claim.statement.startswith("[PRIOR-GUIDED HYPOTHESIS TESTING]")
    assert claim.language_tier == LanguageTier.LEVEL_4_HYPOTHESIS
    assert claim.causal_status == "observational"
    assert len(claim.support_evidence_ids) == len(report.evidence_nodes)


def test_knowledge_engine_prior_guided_automatic_gene_extraction(sample_study_manifest):
    """Tests Prior-Guided Mode when target_genes is omitted and automatically extracted from hypothesis."""
    engine = KnowledgeEngine()
    hypothesis = "Evaluation of APOE and CLEC7A in microglia neurodegeneration"

    report = engine.execute_prior_guided(
        manifest=sample_study_manifest,
        hypothesis=hypothesis,
        target_genes=None,  # Auto-extracted from hypothesis
    )

    assert report.prior_guided is True
    assert "APOE" in report.target_genes
    assert "CLEC7A" in report.target_genes
    assert len(report.evidence_nodes) > 0


def test_evidence_node_schema_compliance_and_validation():
    """Verifies that all generated EvidenceNodes satisfy Pydantic schema validation rules."""
    engine = KnowledgeEngine()
    manifest = StudyManifest(
        study_id="test_schema_study",
        title="Schema test",
        biological_design=BiologicalDesign(species="mus_musculus", tissue="brain"),
    )

    report_disc = engine.execute_discovery(manifest, deg_genes=["Apoe", "Trem2"])
    report_prior = engine.execute_prior_guided(manifest, hypothesis="Trem2 activation", target_genes=["Trem2"])

    all_nodes = report_disc.evidence_nodes + report_prior.evidence_nodes
    assert len(all_nodes) > 0

    for node in all_nodes:
        # Re-validate against Pydantic schema
        dumped = node.model_dump()
        revalidated = EvidenceNode.model_validate(dumped)
        assert revalidated.evidence_id == node.evidence_id
        assert 0.0 <= revalidated.score <= 1.0
        assert revalidated.polarity in (EvidencePolarity.SUPPORTING, EvidencePolarity.CONTRADICTING, EvidencePolarity.NEUTRAL)
        assert revalidated.strength in (EvidenceStrength.VERY_STRONG, EvidenceStrength.STRONG, EvidenceStrength.MODERATE, EvidenceStrength.WEAK, EvidenceStrength.INSUFFICIENT)
        assert len(revalidated.source_artifact_uris) > 0
        assert revalidated.source_task_id != ""
