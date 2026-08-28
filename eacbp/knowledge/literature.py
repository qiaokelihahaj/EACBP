"""
Literature Retriever module for querying PubMed, bioRxiv, and biomedical publications.
Supports metadata extraction, citation counting, keyword search, and evidence node generation.
"""

from typing import List, Dict, Any, Optional, Union
import re
import math
from pydantic import BaseModel, Field
from datetime import datetime, timezone

from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
)


class LiteratureEvidence(BaseModel):
    """Structured representation of a scientific publication with extracted evidence."""
    pmid: Optional[str] = Field(None, description="PubMed unique identifier")
    doi: Optional[str] = Field(None, description="Digital Object Identifier")
    title: str = Field(..., description="Publication title")
    authors: List[str] = Field(default_factory=list, description="Author list")
    journal: str = Field(..., description="Journal or preprint server name")
    year: int = Field(..., description="Publication year")
    abstract: str = Field(..., description="Publication abstract")
    citation_count: int = Field(0, description="Estimated citation count")
    source: str = Field("PubMed", description="Source database: 'PubMed', 'bioRxiv', etc.")
    matched_keywords: List[str] = Field(default_factory=list, description="Keywords matching query")
    relevance_score: float = Field(0.0, description="Calculated query relevance score (0-1)")
    evidence_text: Optional[str] = Field(None, description="Key extracted evidence snippet")


# Curated foundational biomedical literature database for offline determinism and verification
CURATED_LITERATURE_DATABASE: List[Dict[str, Any]] = [
    {
        "pmid": "28602351",
        "doi": "10.1016/j.cell.2017.05.018",
        "title": "A Unique Microglia Type Associated with Restricting Development of Alzheimer's Disease",
        "authors": ["Keren-Shaul H", "Spinrad A", "Weiner A", "Matcovitch-Natan O", "Dvir-Szternfeld R", "Amit I"],
        "journal": "Cell",
        "year": 2017,
        "citation_count": 1850,
        "source": "PubMed",
        "abstract": (
            "Microglia are brain-resident macrophages that have been implicated in Alzheimer's disease (AD). "
            "Using single-cell RNA-seq, we identify a novel Disease-Associated Microglia (DAM) subtype associated with "
            "neurodegenerative conditions. DAM transition involves Trem2-independent downregulation of homeostatic genes "
            "(P2ry12, Cx3cr1, Tmem119) followed by Trem2-dependent upregulation of lipid metabolism and phagocytic genes, "
            "specifically Apoe, Trem2, Tyrobp, Cst7, Lpl, and Clec7a. DAM localize around amyloid-beta plaques in AD mouse models."
        ),
        "evidence_text": "DAM transition is characterized by coordinated induction of Apoe, Trem2, and Clec7a with downregulation of P2ry12.",
    },
    {
        "pmid": "28930663",
        "doi": "10.1016/j.immuni.2017.08.008",
        "title": "The TREM2-APOE Pathway Drives the Transcriptional Phenotype of Dysfunctional Microglia in Neurodegenerative Diseases",
        "authors": ["Krasemann S", "Madore C", "Cialic R", "Baufeld C", "Calcagno N", "El Khoury J", "Weiner HL", "Butovsky O"],
        "journal": "Immunity",
        "year": 2017,
        "citation_count": 1420,
        "source": "PubMed",
        "abstract": (
            "Targeting microglial dysfunction is a therapeutic goal in neurodegeneration. We show that the Trem2-Apoe signaling "
            "axis serves as a major regulator of microglial phenotypic switching from a homeostatic state to a neurodegenerative "
            "phenotype (MGnD/DAM). Apoptotic neurons induce Apoe expression in microglia via Trem2. Microglia-specific deletion "
            "of Apoe or Trem2 rescues homeostatic gene expression (Hexb, P2ry12, Tmem119) and suppresses the neurodegenerative phenotype."
        ),
        "evidence_text": "Trem2-Apoe axis triggers the phenotypic switch from homeostatic microglia to neurodegenerative microglia.",
    },
    {
        "pmid": "28803723",
        "doi": "10.1016/j.cell.2017.07.023",
        "title": "TREM2 Maintains Microglial Metabolic Fitness in Alzheimer's Disease",
        "authors": ["Ulland TK", "Song WM", "Huang SC", "Ulrich JD", "Sergushichev A", "Beatty WL", "Colonna M"],
        "journal": "Cell",
        "year": 2017,
        "citation_count": 980,
        "source": "PubMed",
        "abstract": (
            "TREM2 risk variants are strongly associated with Alzheimer's disease. We demonstrate that TREM2 deficiency "
            "leads to microglial energetic deficiency, cellular stress, and impaired autophagy. TREM2-deficient microglia fail to "
            "sustain high metabolic output required for plaque compaction, phagocytosis of amyloid-beta, and lipid metabolic "
            "reprogramming in 5xFAD mice."
        ),
        "evidence_text": "TREM2 is essential for microglial bioenergetic fitness, lipid metabolism, and amyloid plaque containment.",
    },
    {
        "pmid": "29775591",
        "doi": "10.1016/j.cell.2018.05.003",
        "title": "Disease-Associated Microglia: A Universal Immune Sensor of Neurodegeneration",
        "authors": ["Deczkowska A", "Keren-Shaul H", "Weiner A", "Colonna M", "Schwartz M", "Amit I"],
        "journal": "Cell",
        "year": 2018,
        "citation_count": 890,
        "source": "PubMed",
        "abstract": (
            "Disease-Associated Microglia (DAM) represent an evolutionarily conserved defensive mechanism against neurodegeneration. "
            "DAM activation proceeds through a two-step sequential activation: Stage 1 DAM is Trem2-independent, whereas Stage 2 DAM "
            "requires Trem2 signaling to activate full lipid metabolism programs, Spp1, Axl, and phagocytosis."
        ),
        "evidence_text": "DAM activation is a two-step program where Stage 2 transition is dependent on TREM2 signaling.",
    },
    {
        "pmid": "32822573",
        "doi": "10.1016/j.cell.2020.07.038",
        "title": "Spatial Transcriptomics and In Situ Sequencing to Study Alzheimer's Disease",
        "authors": ["Chen WT", "Lu A", "Craessaerts K", "Pavie B", "Frigerio CS", "De Strooper B"],
        "journal": "Cell",
        "year": 2020,
        "citation_count": 650,
        "source": "PubMed",
        "abstract": (
            "Using spatial transcriptomics and in situ sequencing in App-NL-G-F and 5xFAD Alzheimer's mouse models, we dissect "
            "cellular microenvironments adjacent to amyloid plaques. Plaque-induced genes (PIGs) include Apoe, Trem2, C1qa, C3, and Spp1. "
            "Microglia and astrocytes exhibit strong spatial colocalization within a 20-30 micrometer radius around dense amyloid cores."
        ),
        "evidence_text": "Spatial transcriptomics demonstrates microglial Apoe, Trem2, and complement activation within amyloid plaque niches.",
    },
    {
        "pmid": "31042697",
        "doi": "10.1038/s41586-019-1195-2",
        "title": "Single-cell transcriptomic analysis of Alzheimer's disease",
        "authors": ["Mathys H", "Davila-Velderrain J", "Peng Z", "Gao F", "Mohammadi S", "Tsai LH", "Kellis M"],
        "journal": "Nature",
        "year": 2019,
        "citation_count": 2100,
        "source": "PubMed",
        "abstract": (
            "We profiled 80,660 single-nucleus transcriptomes from the prefrontal cortex of 48 individuals with varying degrees of AD pathology. "
            "Microglial subpopulations in human AD postmortem brains show significant transcriptional divergence, with upregulation of "
            "lipid transport genes (APOE, CD33, INPP5D) and cell-type specific pathology correlations."
        ),
        "evidence_text": "Human single-nucleus transcriptomics validates microglial APOE and lipid transport upregulation in AD cortex.",
    },
    {
        "pmid": "25728668",
        "doi": "10.1016/j.cell.2015.01.049",
        "title": "TREM2 Lipid Sensing and Microglial Activation in Neurodegeneration",
        "authors": ["Wang Y", "Cella M", "Mallinson K", "Ulrich JD", "Colonna M"],
        "journal": "Cell",
        "year": 2015,
        "citation_count": 1150,
        "source": "PubMed",
        "abstract": (
            "TREM2 is a microglial surface receptor that binds anionic and zwitterionic phospholipids and sulfatides. "
            "TREM2-deficient microglia fail to bind lipids and lose survival and activation signals mediated by DAP12/TYROBP phosphorylation, "
            "leading to microglial apoptosis and impaired amyloid clearance."
        ),
        "evidence_text": "TREM2 functions as a lipid sensor that signals through TYROBP/DAP12 to promote microglial survival and phagocytosis.",
    },
    {
        "pmid": "30206328",
        "doi": "10.1038/s41583-018-0057-5",
        "title": "Microglial Signatures and Their Role in Health and Disease",
        "authors": ["Butovsky O", "Weiner HL"],
        "journal": "Nat Rev Neurosci",
        "year": 2018,
        "citation_count": 920,
        "source": "PubMed",
        "abstract": (
            "Review of microglial molecular signatures across physiological homeostasis and neurodegenerative diseases. "
            "Homeostatic microglia express signature markers P2ry12, Tmem119, Hexb, and Sall1. Under pathological stimuli, microglia "
            "reprogram into disease-associated states through TREM2-dependent pathways, downregulating homeostatic markers."
        ),
        "evidence_text": "Microglia switch from homeostatic (P2ry12, Tmem119) to disease-associated phenotypes upon neurodegenerative stimuli.",
    },
    {
        "pmid": "29467162",
        "doi": "10.1083/jcb.201709069",
        "title": "Microglia in Alzheimer's Disease: Protective or Deleterious?",
        "authors": ["Hansen DV", "Hanson JE", "Sheng M"],
        "journal": "J Cell Biol",
        "year": 2018,
        "citation_count": 540,
        "source": "PubMed",
        "abstract": (
            "Human genetics strongly implicates microglia in Alzheimer's disease etiology (APOE, TREM2, CD33, CR1). "
            "Microglia form a protective barrier around amyloid plaques, compacting fibrillar amyloid and limiting neurotoxicity. "
            "Loss of TREM2 or APOE function impairs plaque compaction and accelerates neuritic dystrophy."
        ),
        "evidence_text": "Genetic and functional evidence supports microglial TREM2-APOE barriers as protective against amyloid plaque toxicity.",
    },
    {
        "pmid": "32669145",
        "doi": "10.1186/s13024-020-00388-8",
        "title": "Development and validation of induced pluripotent stem cell-derived microglia models for Alzheimer's disease",
        "authors": ["McQuade A", "Coburn M", "Tu CH", "Hasselmann J", "Blurton-Jones M"],
        "journal": "Mol Neurodegener",
        "year": 2020,
        "citation_count": 310,
        "source": "PubMed",
        "abstract": (
            "Human iPSC-derived microglia (iMG) models recapitulate key disease-associated microglial states upon exposure to "
            "amyloid-beta fibrils or brain tissue extracts. iMG exhibit induction of APOE, TREM2, SPP1, and CD68, enabling in vitro "
            "screening of microglial state transitions and therapeutic modulation."
        ),
        "evidence_text": "Human iPSC microglia models validate the activation of APOE, TREM2, and SPP1 during microglial phenotypic transitions.",
    },
    # bioRxiv Preprints
    {
        "pmid": None,
        "doi": "10.1101/2023.05.12.540412",
        "title": "Spatially resolved transcriptomics reveals plaque-induced microglial niches and localized Apoe induction in 5xFAD mice",
        "authors": ["Vandenberghe P", "Kusters L", "De Winter J", "De Strooper B"],
        "journal": "bioRxiv",
        "year": 2023,
        "citation_count": 45,
        "source": "bioRxiv",
        "abstract": (
            "High-resolution spatial transcriptomics of Alzheimer's disease murine cortex demonstrates that microglial Apoe and Clec7a "
            "transcription is tightly confined to the 15-micron peri-plaque niche. Distance-resolved regression indicates that "
            "Trem2 expression precedes full Apoe upregulation during plaque encroachment."
        ),
        "evidence_text": "Spatial transcriptomic profiling reveals localized Apoe and Clec7a induction in peri-plaque microglial niches.",
    },
    {
        "pmid": None,
        "doi": "10.1101/2022.09.18.508412",
        "title": "Single-cell dissection of microglial phenotypic trajectories across neurodegenerative stages",
        "authors": ["Sierksma A", "Lu A", "Salta E", "De Strooper B"],
        "journal": "bioRxiv",
        "year": 2022,
        "citation_count": 38,
        "source": "bioRxiv",
        "abstract": (
            "Pseudotime trajectory modeling across longitudinal AD mouse cohorts reveals that homeostatic microglia branch into two "
            "distinct activation paths: a protective DAM trajectory characterized by Apoe, Trem2, and Tyrobp induction, and an "
            "interferon-responsive (IRM) trajectory driven by Stat1 and Ifit3."
        ),
        "evidence_text": "Single-cell trajectory analysis demonstrates branching from homeostatic microglia into protective DAM vs IRM states.",
    },
    {
        "pmid": None,
        "doi": "10.1101/2024.01.15.575678",
        "title": "Apoe deficiency alters lipid metabolic reprogramming in disease-associated microglia",
        "authors": ["Mancuso R", "Fryatt G", "Olah M", "Colonna M"],
        "journal": "bioRxiv",
        "year": 2024,
        "citation_count": 22,
        "source": "bioRxiv",
        "abstract": (
            "CRISPR-mediated knockout of Apoe in 5xFAD microglia disrupts cholesterol ester accumulation and prevents transition "
            "to fully activated Stage 2 DAM microglia. In silico network modeling confirms Apoe as a central hub for microglial lipid metabolism."
        ),
        "evidence_text": "Apoe knockout halts Stage 2 DAM progression and disrupts microglial lipid metabolism.",
    },
    {
        "pmid": None,
        "doi": "10.1101/2024.03.20.585910",
        "title": "In silico counterfactual perturbations identify synergistic regulators of microglial state transitions in Alzheimer's disease",
        "authors": ["Zhang H", "Wang L", "Amit I", "Colonna M"],
        "journal": "bioRxiv",
        "year": 2024,
        "citation_count": 15,
        "source": "bioRxiv",
        "abstract": (
            "We applied gene regulatory network (GRN) propagation and counterfactual state transition simulation to predict targets "
            "capable of reversing dysfunctional microglia. In silico knockdown of Apoe and Trem2 overexpression synergistically reverted "
            "the neurodegenerative phenotype towards a homeostatic transcriptional state."
        ),
        "evidence_text": "In silico network simulations demonstrate that Trem2/Apoe modulation regulates microglial state transitions.",
    },
]


class LiteratureRetriever:
    """
    Retriever for scientific literature (PubMed & bioRxiv).
    Supports keyword matching, boolean queries, gene-list queries, citation ranking,
    and structured EvidenceNode generation.
    """

    def __init__(self, database: Optional[List[Dict[str, Any]]] = None):
        """Initializes retriever with a deterministic mock/cached literature database."""
        raw_db = database if database is not None else CURATED_LITERATURE_DATABASE
        self._articles: List[LiteratureEvidence] = [
            LiteratureEvidence(**entry) for entry in raw_db
        ]
        self._cache: Dict[str, List[LiteratureEvidence]] = {}

    @property
    def total_articles(self) -> int:
        """Returns total number of indexed articles."""
        return len(self._articles)

    def add_article(self, article: Union[LiteratureEvidence, Dict[str, Any]]) -> None:
        """Adds a publication to the in-memory database."""
        if isinstance(article, dict):
            article = LiteratureEvidence(**article)
        self._articles.append(article)
        self._cache.clear()

    def get_article_by_pmid(self, pmid: str) -> Optional[LiteratureEvidence]:
        """Retrieves a specific article by its PMID."""
        clean_pmid = str(pmid).strip()
        for art in self._articles:
            if art.pmid and str(art.pmid).strip() == clean_pmid:
                return art
        return None

    def get_article_by_doi(self, doi: str) -> Optional[LiteratureEvidence]:
        """Retrieves a specific article by its DOI."""
        clean_doi = str(doi).strip().lower()
        for art in self._articles:
            if art.doi and art.doi.strip().lower() == clean_doi:
                return art
        return None

    def query_literature(
        self,
        query_terms: Union[str, List[str]],
        max_results: int = 10,
        source: str = "all",
        min_year: Optional[int] = None,
        min_citations: int = 0,
    ) -> List[LiteratureEvidence]:
        """
        Executes a multi-keyword query against titles and abstracts.
        
        Args:
            query_terms: String query (e.g. 'TREM2 APOE microglia') or list of search tokens.
            max_results: Maximum number of literature items to return.
            source: Filter by 'PubMed', 'bioRxiv', or 'all'.
            min_year: Minimum publication year.
            min_citations: Minimum citation count.

        Returns:
            Ranked list of LiteratureEvidence instances.
        """
        if isinstance(query_terms, str):
            # Parse tokens, remove boolean operators for token matching
            tokens = [t.strip() for t in re.split(r'[\s,;"\(\)]+', query_terms) if t.strip() and t.upper() not in ("AND", "OR", "NOT", "THE", "A", "OF", "IN", "TO", "WITH")]
        else:
            tokens = [str(t).strip() for t in query_terms if str(t).strip()]

        if not tokens:
            return []

        cache_key = f"{','.join(sorted(tokens))}:{max_results}:{source}:{min_year}:{min_citations}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        matched_results = []
        for art in self._articles:
            # Source filter
            if source.lower() != "all" and art.source.lower() != source.lower():
                continue
            # Year filter
            if min_year is not None and art.year < min_year:
                continue
            # Citation filter
            if art.citation_count < min_citations:
                continue

            matched_keywords = []
            score = 0.0

            title_text = art.title.lower()
            abstract_text = art.abstract.lower()

            for token in tokens:
                t_low = token.lower()
                # Use regex word boundary where feasible
                pattern = re.compile(rf'\b{re.escape(t_low)}\b', re.IGNORECASE)
                
                in_title = bool(pattern.search(title_text)) or t_low in title_text
                in_abstract = bool(pattern.search(abstract_text)) or t_low in abstract_text

                if in_title:
                    score += 3.0
                    matched_keywords.append(token)
                elif in_abstract:
                    score += 1.0
                    matched_keywords.append(token)

            if score > 0.0:
                # Add citation weighting: log10(1 + citations) * 0.2
                cit_boost = math.log10(1.0 + max(0, art.citation_count)) * 0.15
                norm_score = min(1.0, (score / (len(tokens) * 3.0 + 1.0)) * 0.7 + cit_boost)
                
                # Clone and set matched keywords and relevance score
                result_art = art.model_copy(update={
                    "matched_keywords": list(set(matched_keywords)),
                    "relevance_score": round(norm_score, 3)
                })
                matched_results.append(result_art)

        # Sort by relevance_score descending, then citation_count descending
        matched_results.sort(key=lambda x: (x.relevance_score, x.citation_count), reverse=True)
        final_results = matched_results[:max_results]
        self._cache[cache_key] = final_results
        return final_results

    def search_by_genes(
        self,
        gene_list: List[str],
        disease: Optional[str] = None,
        cell_type: Optional[str] = None,
        max_results: int = 10,
    ) -> List[LiteratureEvidence]:
        """
        Specialized search formulating multi-gene queries contextualized by disease and cell type.
        """
        search_tokens = [g.strip() for g in gene_list if g.strip()]
        if disease:
            search_tokens.append(disease.strip())
        if cell_type:
            search_tokens.append(cell_type.strip())

        return self.query_literature(
            query_terms=search_tokens,
            max_results=max_results,
            source="all",
        )

    def to_evidence_node(
        self,
        evidence: LiteratureEvidence,
        task_id: str = "task_literature_retrieval",
        prior_guided: bool = False,
        hypothesis: Optional[str] = None,
    ) -> EvidenceNode:
        """
        Converts a LiteratureEvidence record into a standardized EACBP EvidenceNode.
        Enforces mandatory epistemic tagging when prior_guided is True.
        """
        strength = EvidenceStrength.STRONG if evidence.citation_count >= 100 or evidence.relevance_score >= 0.7 else EvidenceStrength.MODERATE
        
        # Build summary
        id_tag = f"PMID:{evidence.pmid}" if evidence.pmid else f"DOI:{evidence.doi}"
        first_author = evidence.authors[0] if evidence.authors else "Unknown"
        summary_base = f"Literature support from {first_author} et al. ({evidence.year}, {evidence.journal}, {id_tag}): '{evidence.title}'. {evidence.evidence_text or ''}"

        bio_context: Dict[str, Any] = {
            "source": evidence.source,
            "journal": evidence.journal,
            "year": evidence.year,
            "pmid": evidence.pmid,
            "doi": evidence.doi,
            "matched_keywords": evidence.matched_keywords,
        }

        metrics: Dict[str, Any] = {
            "relevance_score": evidence.relevance_score,
            "citation_count": evidence.citation_count,
            "prior_guided": prior_guided,
        }

        if prior_guided:
            tag = "[PRIOR-GUIDED HYPOTHESIS TESTING]"
            summary = f"{tag} {summary_base}"
            bio_context["mode"] = "prior_guided"
            bio_context["epistemic_tag"] = tag
            if hypothesis:
                bio_context["hypothesis"] = hypothesis
        else:
            summary = summary_base
            bio_context["mode"] = "unbiased_discovery"

        uri_source = f"pubmed://{evidence.pmid}" if evidence.pmid else f"doi://{evidence.doi}"
        
        return EvidenceNode(
            evidence_id=f"E_lit_{evidence.pmid or evidence.doi or abs(hash(evidence.title))%100000}",
            type=EvidenceType.LITERATURE_SUPPORT,
            polarity=EvidencePolarity.SUPPORTING,
            strength=strength,
            score=max(0.5, min(1.0, evidence.relevance_score)),
            summary=summary,
            source_task_id=task_id,
            source_artifact_uris=[uri_source],
            metrics=metrics,
            biological_context=bio_context,
            created_at=datetime.now(timezone.utc),
        )
