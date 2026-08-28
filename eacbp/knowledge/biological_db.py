"""
Biological Database Retriever module.
Provides Gene Ontology (GO) over-representation analysis (ORA) with hypergeometric tests and
Benjamini-Hochberg FDR correction, Reactome pathway mapping, and NCBI Gene / ortholog resolution.
"""

from typing import List, Dict, Any, Optional, Union, Tuple
from pydantic import BaseModel, Field
from datetime import datetime, timezone
import scipy.stats as stats

from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
)


class GOEnrichment(BaseModel):
    """Result of Gene Ontology Over-Representation Analysis."""
    go_id: str = Field(..., description="GO term identifier (e.g. GO:0006954)")
    term: str = Field(..., description="GO term name")
    category: str = Field(..., description="GO category: BP, MF, or CC")
    gene_count: int = Field(..., description="Number of query genes in this term (k)")
    term_size: int = Field(..., description="Total annotated genes in this term (K)")
    background_size: int = Field(20000, description="Total background universe genes (M)")
    sample_size: int = Field(..., description="Total query genes in analysis (n)")
    fold_enrichment: float = Field(..., description="(k/n) / (K/M)")
    p_value: float = Field(..., description="Hypergeometric p-value (survival function)")
    fdr_q_value: float = Field(..., description="Benjamini-Hochberg adjusted p-value (FDR)")
    genes: List[str] = Field(default_factory=list, description="List of matched query genes")


class PathwayEnrichment(BaseModel):
    """Result of Pathway Over-Representation Analysis."""
    pathway_id: str = Field(..., description="Pathway identifier (e.g. R-HSA-2173782)")
    pathway_name: str = Field(..., description="Pathway name")
    source: str = Field("Reactome", description="Pathway database source")
    gene_count: int = Field(..., description="Number of query genes in pathway (k)")
    pathway_size: int = Field(..., description="Total genes in pathway (K)")
    background_size: int = Field(20000, description="Total background universe genes (M)")
    sample_size: int = Field(..., description="Total query genes in analysis (n)")
    fold_enrichment: float = Field(..., description="(k/n) / (K/M)")
    p_value: float = Field(..., description="Hypergeometric p-value")
    fdr_q_value: float = Field(..., description="Benjamini-Hochberg adjusted p-value (FDR)")
    genes: List[str] = Field(default_factory=list, description="List of matched query genes")


class GeneInfo(BaseModel):
    """NCBI Gene entity and cross-species ortholog annotation."""
    gene_symbol: str = Field(..., description="Official gene symbol")
    gene_id: str = Field(..., description="NCBI Gene ID")
    species: str = Field(..., description="Species, e.g. 'homo_sapiens', 'mus_musculus'")
    full_name: str = Field(..., description="Full descriptive gene name")
    aliases: List[str] = Field(default_factory=list, description="Known gene aliases and synonyms")
    summary: str = Field(..., description="Functional biological summary")
    orthologs: Dict[str, str] = Field(default_factory=dict, description="Cross-species ortholog mapping {'human': '...', 'mouse': '...'}")
    go_terms: List[str] = Field(default_factory=list, description="Associated key GO terms")
    pathways: List[str] = Field(default_factory=list, description="Associated key pathways")


# Curated Gene Ontology Term Database
CURATED_GO_TERMS: List[Dict[str, Any]] = [
    # Biological Process (BP)
    {
        "go_id": "GO:0006954",
        "term": "inflammatory response",
        "category": "BP",
        "genes": ["APOE", "TREM2", "TYROBP", "CLEC7A", "C3", "C1QA", "C1QB", "C1QC", "IL1B", "TNF", "CX3CR1", "TLR4", "AIF1", "CD68", "SPP1", "CST7"],
        "term_size": 650,
    },
    {
        "go_id": "GO:0006911",
        "term": "phagocytosis, engulfment",
        "category": "BP",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "C3", "C1QA", "C1QB", "C1QC", "CD68", "AXL", "MERTK", "ITGAX", "APOE", "LPL"],
        "term_size": 320,
    },
    {
        "go_id": "GO:0006629",
        "term": "lipid metabolic process",
        "category": "BP",
        "genes": ["APOE", "TREM2", "LPL", "LIPA", "ABCA1", "ABCG1", "APOC1", "FABP5", "CH25H", "SOAT1", "LDLR", "LRP1"],
        "term_size": 980,
    },
    {
        "go_id": "GO:0050727",
        "term": "regulation of inflammatory response",
        "category": "BP",
        "genes": ["TREM2", "TYROBP", "APOE", "IL10", "TGFB1", "SOCS3", "CX3CL1", "CX3CR1", "C3", "CD68"],
        "term_size": 420,
    },
    {
        "go_id": "GO:0001775",
        "term": "cell activation",
        "category": "BP",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "SPP1", "ITGAX", "CD68", "CD14", "AIF1", "CSF1R"],
        "term_size": 850,
    },
    {
        "go_id": "GO:0030198",
        "term": "extracellular matrix organization",
        "category": "BP",
        "genes": ["SPP1", "MMP9", "MMP12", "COL1A1", "COL1A2", "FN1", "TIMP1", "LRP1"],
        "term_size": 480,
    },
    {
        "go_id": "GO:0006956",
        "term": "complement activation",
        "category": "BP",
        "genes": ["C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFH", "CR1", "C3AR1"],
        "term_size": 190,
    },
    {
        "go_id": "GO:0042552",
        "term": "myelination",
        "category": "BP",
        "genes": ["MBP", "MOG", "PLP1", "MAG", "CNP", "CLDN11"],
        "term_size": 210,
    },
    {
        "go_id": "GO:0007268",
        "term": "chemical synaptic transmission",
        "category": "BP",
        "genes": ["SYN1", "SYP", "GRIN1", "GRIN2B", "GRIA1", "SNAP25", "VAMP2", "GABRA1", "GABBR1"],
        "term_size": 780,
    },
    {
        "go_id": "GO:0048699",
        "term": "generation of neurons",
        "category": "BP",
        "genes": ["NEUROD1", "DCX", "SOX2", "PAX6", "ASCL1", "MAP2", "RBFOX3"],
        "term_size": 1100,
    },
    {
        "go_id": "GO:0001819",
        "term": "positive regulation of cytokine production",
        "category": "BP",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "TLR4", "IL1B", "TNF", "SPP1", "AIF1"],
        "term_size": 530,
    },
    {
        "go_id": "GO:0043066",
        "term": "negative regulation of apoptotic process",
        "category": "BP",
        "genes": ["TREM2", "TYROBP", "BCL2", "AKT1", "APOE", "CSF1R", "AXL"],
        "term_size": 620,
    },
    # Molecular Function (MF)
    {
        "go_id": "GO:0005543",
        "term": "phospholipid binding",
        "category": "MF",
        "genes": ["TREM2", "APOE", "LPL", "SYT1", "ANXA1", "ANXA5"],
        "term_size": 240,
    },
    {
        "go_id": "GO:0005044",
        "term": "scavenger receptor activity",
        "category": "MF",
        "genes": ["TREM2", "CD68", "MSR1", "MARCO", "CD36", "SCARB1", "LRP1"],
        "term_size": 110,
    },
    {
        "go_id": "GO:0005125",
        "term": "cytokine activity",
        "category": "MF",
        "genes": ["IL1B", "TNF", "IL6", "TGFB1", "CXCL12", "CCL2", "SPP1"],
        "term_size": 290,
    },
    {
        "go_id": "GO:0005515",
        "term": "protein binding",
        "category": "MF",
        "genes": ["APOE", "TREM2", "TYROBP", "CLEC7A", "C3", "SPP1", "CX3CR1", "CD68", "LPL", "CST7"],
        "term_size": 4500,
    },
    # Cellular Component (CC)
    {
        "go_id": "GO:0005886",
        "term": "plasma membrane",
        "category": "CC",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "CX3CR1", "P2RY12", "TMEM119", "ITGAX", "CD68", "CSF1R", "AXL"],
        "term_size": 3800,
    },
    {
        "go_id": "GO:0005764",
        "term": "lysosome",
        "category": "CC",
        "genes": ["CD68", "HEXB", "CTSB", "CTSD", "LAMP1", "LAMP2", "GBA", "NPC1", "CST7", "LIPA"],
        "term_size": 510,
    },
    {
        "go_id": "GO:0005576",
        "term": "extracellular region",
        "category": "CC",
        "genes": ["APOE", "C1QA", "C1QB", "C1QC", "C3", "SPP1", "CST7", "LPL", "TGFB1", "TNF"],
        "term_size": 2800,
    },
    {
        "go_id": "GO:0016575",
        "term": "histone H4-K16 acetylation",
        "category": "BP",
        "genes": ["KAT8", "KANSL1", "MSL1", "MSL2", "MSL3", "H4C1", "KAT5", "EP300"],
        "term_size": 120,
    },
    {
        "go_id": "GO:0006325",
        "term": "chromatin organization and remodeling",
        "category": "BP",
        "genes": ["KAT8", "KANSL1", "MSL1", "SMARCA4", "ARID1A", "CHD4", "BAX", "TRP53", "H2AX"],
        "term_size": 480,
    },
    {
        "go_id": "GO:0000077",
        "term": "DNA damage checkpoint and response",
        "category": "BP",
        "genes": ["KAT8", "TP53", "TRP53", "H2AX", "CHEK1", "CHEK2", "ATM", "ATR", "CDKN1A", "BAX"],
        "term_size": 350,
    },
    {
        "go_id": "GO:0007049",
        "term": "cell cycle mitotic progression",
        "category": "BP",
        "genes": ["CDK1", "TOP2A", "CCNB1", "CDKN1A", "MKI67", "AURKA", "KAT8"],
        "term_size": 620,
    },
    {
        "go_id": "GO:0006915",
        "term": "apoptotic process regulation",
        "category": "BP",
        "genes": ["BAX", "BCL2", "CASP3", "CASP9", "TRP53", "CDKN1A", "CYCS"],
        "term_size": 750,
    },
    {
        "go_id": "GO:0004402",
        "term": "histone acetyltransferase activity",
        "category": "MF",
        "genes": ["KAT8", "KAT5", "EP300", "CREBBP", "KAT2A", "KAT6A"],
        "term_size": 85,
    },
]


# Curated Reactome Pathways Database
CURATED_REACTOME_PATHWAYS: List[Dict[str, Any]] = [
    {
        "pathway_id": "R-HSA-2173782",
        "pathway_name": "Microglia Pathogen Phagocytosis",
        "source": "Reactome",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "C1QA", "C1QB", "C1QC", "C3", "CD68", "AXL", "MERTK", "APOE"],
        "pathway_size": 140,
    },
    {
        "pathway_id": "R-HSA-168256",
        "pathway_name": "Innate Immune System",
        "source": "Reactome",
        "genes": ["TREM2", "TYROBP", "CLEC7A", "C1QA", "C1QB", "C1QC", "C3", "TLR4", "IL1B", "TNF", "CX3CR1", "CD68", "SPP1"],
        "pathway_size": 1200,
    },
    {
        "pathway_id": "R-HSA-174824",
        "pathway_name": "Plasma Lipoprotein Assembly, Remodeling, and Clearance",
        "source": "Reactome",
        "genes": ["APOE", "LPL", "LIPA", "ABCA1", "ABCG1", "APOC1", "LDLR", "LRP1", "CH25H"],
        "pathway_size": 160,
    },
    {
        "pathway_id": "R-HSA-166658",
        "pathway_name": "Complement Cascade",
        "source": "Reactome",
        "genes": ["C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFH", "C3AR1", "C5AR1", "CR1"],
        "pathway_size": 110,
    },
    {
        "pathway_id": "R-HSA-380108",
        "pathway_name": "Chemokine Receptors Bind Chemokines",
        "source": "Reactome",
        "genes": ["CX3CR1", "CX3CL1", "CCL2", "CCR2", "CXCL12", "CXCR4"],
        "pathway_size": 95,
    },
    {
        "pathway_id": "R-HSA-5688426",
        "pathway_name": "Deubiquitination",
        "source": "Reactome",
        "genes": ["USP7", "USP8", "USP14", "UCHL1", "BAP1", "STAMBP"],
        "pathway_size": 230,
    },
    {
        "pathway_id": "R-HSA-112316",
        "pathway_name": "Neuronal System",
        "source": "Reactome",
        "genes": ["SYN1", "SYP", "GRIN1", "GRIN2B", "GRIA1", "SNAP25", "VAMP2", "GABRA1", "GABBR1"],
        "pathway_size": 540,
    },
    {
        "pathway_id": "R-HSA-6798695",
        "pathway_name": "Neutrophil Degranulation",
        "source": "Reactome",
        "genes": ["S100A8", "S100A9", "MMP9", "MPO", "ITGAM", "CD14", "C3"],
        "pathway_size": 480,
    },
    {
        "pathway_id": "R-HSA-2132295",
        "pathway_name": "MHC Class II Antigen Presentation",
        "source": "Reactome",
        "genes": ["HLA-DRA", "HLA-DRB1", "HLA-DQA1", "CD74", "CTSS", "IFI30"],
        "pathway_size": 125,
    },
    {
        "pathway_id": "R-HSA-9006934",
        "pathway_name": "Signaling by Receptor Tyrosine Kinases",
        "source": "Reactome",
        "genes": ["CSF1R", "AXL", "MERTK", "EGFR", "FGFR1", "PDGFRB", "FLT1"],
        "pathway_size": 450,
    },
    {
        "pathway_id": "R-HSA-3247509",
        "pathway_name": "Chromatin modifying enzymes and Histone Acetylation",
        "source": "Reactome",
        "genes": ["KAT8", "KANSL1", "MSL1", "MSL2", "MSL3", "EP300", "CREBBP", "KAT5", "H4C1"],
        "pathway_size": 180,
    },
    {
        "pathway_id": "R-HSA-69278",
        "pathway_name": "Cell Cycle Checkpoints and DNA Replication",
        "source": "Reactome",
        "genes": ["CDK1", "TOP2A", "CDKN1A", "CCNB1", "TP53", "TRP53", "KAT8", "H2AX"],
        "pathway_size": 290,
    },
    {
        "pathway_id": "R-HSA-109581",
        "pathway_name": "Apoptosis and Programmed Cell Death",
        "source": "Reactome",
        "genes": ["BAX", "BCL2", "CASP3", "CASP9", "TRP53", "CYCS", "CDKN1A"],
        "pathway_size": 210,
    },
]


# Curated NCBI Gene & Ortholog Database
CURATED_NCBI_GENES: Dict[str, Dict[str, Any]] = {
    "KAT8": {
        "gene_symbol": "KAT8",
        "gene_id": "84148",
        "species": "homo_sapiens",
        "full_name": "lysine acetyltransferase 8",
        "aliases": ["MOF", "MYST1", "hMOF"],
        "summary": "Catalytic subunit of MSL and NSL complexes; specifically mediates histone H4 lysine 16 acetylation (H4K16ac) to promote chromatin accessibility and transcriptional activation.",
        "orthologs": {"human": "KAT8", "mouse": "Kat8", "rat": "Kat8"},
        "go_terms": ["GO:0016575", "GO:0006325", "GO:0000077", "GO:0004402"],
        "pathways": ["R-HSA-3247509", "R-HSA-69278"],
    },
    "KAT8_MOUSE": {
        "gene_symbol": "Kat8",
        "gene_id": "67957",
        "species": "mus_musculus",
        "full_name": "lysine acetyltransferase 8",
        "aliases": ["Mof", "Myst1", "AU043380"],
        "summary": "Mouse ortholog of KAT8; essential for H4K16ac epigenetic mark, stem cell self-renewal, and neural/germline lineage differentiation.",
        "orthologs": {"human": "KAT8", "mouse": "Kat8", "rat": "Kat8"},
        "go_terms": ["GO:0016575", "GO:0006325", "GO:0000077", "GO:0004402"],
        "pathways": ["R-HSA-3247509", "R-HSA-69278"],
    },
    "KANSL1_MOUSE": {
        "gene_symbol": "Kansl1",
        "gene_id": "226535",
        "species": "mus_musculus",
        "full_name": "KAT8 regulatory NSL complex subunit 1",
        "aliases": ["Kdvs", "mKIAA1267"],
        "summary": "Subunit of the NSL histone acetyltransferase complex that directs Kat8 to active promoters and coordinates global transcriptional regulation.",
        "orthologs": {"human": "KANSL1", "mouse": "Kansl1"},
        "go_terms": ["GO:0016575", "GO:0006325"],
        "pathways": ["R-HSA-3247509"],
    },
    "CDKN1A_MOUSE": {
        "gene_symbol": "Cdkn1a",
        "gene_id": "12575",
        "species": "mus_musculus",
        "full_name": "cyclin-dependent kinase inhibitor 1A (P21)",
        "aliases": ["p21", "CIP1", "WAF1"],
        "summary": "Potent cyclin-dependent kinase inhibitor that mediates p53-dependent G1/S cell cycle arrest and cellular stress response.",
        "orthologs": {"human": "CDKN1A", "mouse": "Cdkn1a"},
        "go_terms": ["GO:0007049", "GO:0000077"],
        "pathways": ["R-HSA-69278"],
    },
    "BAX_MOUSE": {
        "gene_symbol": "Bax",
        "gene_id": "12028",
        "species": "mus_musculus",
        "full_name": "BCL2-associated X protein",
        "aliases": ["bax"],
        "summary": "Pro-apoptotic BCL-2 family member that promotes mitochondrial outer membrane permeabilization and apoptosis in response to DNA damage or chromatin stress.",
        "orthologs": {"human": "BAX", "mouse": "Bax"},
        "go_terms": ["GO:0006915", "GO:0000077"],
        "pathways": ["R-HSA-109581"],
    },
    "CDK1_MOUSE": {
        "gene_symbol": "Cdk1",
        "gene_id": "12534",
        "species": "mus_musculus",
        "full_name": "cyclin-dependent kinase 1",
        "aliases": ["Cdc2", "Cdc2a"],
        "summary": "Key catalytic kinase driving the G2/M phase transition of the cell cycle; downregulated upon loss of epigenetic transcriptional maintenance.",
        "orthologs": {"human": "CDK1", "mouse": "Cdk1"},
        "go_terms": ["GO:0007049"],
        "pathways": ["R-HSA-69278"],
    },
    "APOE": {
        "gene_symbol": "APOE",
        "gene_id": "348",
        "species": "homo_sapiens",
        "full_name": "apolipoprotein E",
        "aliases": ["AD2", "LPG", "APO-E"],
        "summary": "Major apolipoprotein of the chylomicron; essential for the normal catabolism of triglyceride-rich lipoprotein constituents and implicated as a critical risk factor in Alzheimer's disease.",
        "orthologs": {"human": "APOE", "mouse": "Apoe", "rat": "Apoe"},
        "go_terms": ["GO:0006629", "GO:0006954", "GO:0005576"],
        "pathways": ["R-HSA-174824", "R-HSA-2173782"],
    },
    "APOE_MOUSE": {
        "gene_symbol": "Apoe",
        "gene_id": "11816",
        "species": "mus_musculus",
        "full_name": "apolipoprotein E",
        "aliases": ["Apo-E", "AI255918"],
        "summary": "Mouse ortholog of APOE; highly induced in disease-associated microglia and lipid-accumulating reactive macrophages.",
        "orthologs": {"human": "APOE", "mouse": "Apoe", "rat": "Apoe"},
        "go_terms": ["GO:0006629", "GO:0006954", "GO:0005576"],
        "pathways": ["R-HSA-174824", "R-HSA-2173782"],
    },
    "TREM2": {
        "gene_symbol": "TREM2",
        "gene_id": "54209",
        "species": "homo_sapiens",
        "full_name": "triggering receptor expressed on myeloid cells 2",
        "aliases": ["PLOSL2", "TREM-2", "HDLS"],
        "summary": "Innate immune receptor on myeloid cells that binds phospholipids, sulfatides, and amyloid-beta; signals through TYROBP/DAP12 to sustain microglial survival and phagocytosis.",
        "orthologs": {"human": "TREM2", "mouse": "Trem2", "rat": "Trem2"},
        "go_terms": ["GO:0006911", "GO:0006954", "GO:0005543", "GO:0005886"],
        "pathways": ["R-HSA-2173782", "R-HSA-168256"],
    },
    "TREM2_MOUSE": {
        "gene_symbol": "Trem2",
        "gene_id": "83433",
        "species": "mus_musculus",
        "full_name": "triggering receptor expressed on myeloid cells 2",
        "aliases": ["Trem2a", "Trem2b", "Trem2c"],
        "summary": "Mouse ortholog of TREM2; master regulator of Stage 2 DAM microglial transition and plaque containment.",
        "orthologs": {"human": "TREM2", "mouse": "Trem2", "rat": "Trem2"},
        "go_terms": ["GO:0006911", "GO:0006954", "GO:0005543", "GO:0005886"],
        "pathways": ["R-HSA-2173782", "R-HSA-168256"],
    },
    "TYROBP": {
        "gene_symbol": "TYROBP",
        "gene_id": "7305",
        "species": "homo_sapiens",
        "full_name": "TYRO protein tyrosine kinase binding protein",
        "aliases": ["DAP12", "KARAP", "PLOSL"],
        "summary": "Transmembrane signaling adaptor containing an ITAM motif that associates with TREM2 and CLEC7A to transduce activation and survival signals in microglia.",
        "orthologs": {"human": "TYROBP", "mouse": "Tyrobp", "rat": "Tyrobp"},
        "go_terms": ["GO:0006954", "GO:0006911", "GO:0005886"],
        "pathways": ["R-HSA-2173782", "R-HSA-168256"],
    },
    "CLEC7A": {
        "gene_symbol": "CLEC7A",
        "gene_id": "64581",
        "species": "homo_sapiens",
        "full_name": "C-type lectin domain containing 7A",
        "aliases": ["BGR", "CLECSF12", "DECTIN1"],
        "summary": "Pattern recognition receptor and major marker of disease-associated microglia; recognizes beta-glucans and damage-associated molecular patterns.",
        "orthologs": {"human": "CLEC7A", "mouse": "Clec7a", "rat": "Clec7a"},
        "go_terms": ["GO:0006954", "GO:0006911", "GO:0005886"],
        "pathways": ["R-HSA-2173782", "R-HSA-168256"],
    },
    "SPP1": {
        "gene_symbol": "SPP1",
        "gene_id": "6696",
        "species": "homo_sapiens",
        "full_name": "secreted phosphoprotein 1",
        "aliases": ["BNSP", "OPN", "BSPI", "ETA-1"],
        "summary": "Osteopontin; cytokine and extracellular matrix protein strongly upregulated in activated microglia and peri-plaque niches.",
        "orthologs": {"human": "SPP1", "mouse": "Spp1", "rat": "Spp1"},
        "go_terms": ["GO:0030198", "GO:0001775", "GO:0005576"],
        "pathways": ["R-HSA-168256"],
    },
    "C3": {
        "gene_symbol": "C3",
        "gene_id": "718",
        "species": "homo_sapiens",
        "full_name": "complement C3",
        "aliases": ["ARMD9", "CPAMD1"],
        "summary": "Central component of the complement system; plays key roles in opsonization, microglial phagocytic pruning, and neuroinflammation.",
        "orthologs": {"human": "C3", "mouse": "C3", "rat": "C3"},
        "go_terms": ["GO:0006956", "GO:0006911", "GO:0005576"],
        "pathways": ["R-HSA-166658", "R-HSA-2173782"],
    },
    "C1QA": {
        "gene_symbol": "C1QA",
        "gene_id": "712",
        "species": "homo_sapiens",
        "full_name": "complement C1q A chain",
        "aliases": ["C1Q-A"],
        "summary": "Initiator of the classical complement pathway, secreted by microglia to tag synapses and apoptotic debris for phagocytosis.",
        "orthologs": {"human": "C1QA", "mouse": "C1qa", "rat": "C1qa"},
        "go_terms": ["GO:0006956", "GO:0006911", "GO:0005576"],
        "pathways": ["R-HSA-166658", "R-HSA-2173782"],
    },
    "CX3CR1": {
        "gene_symbol": "CX3CR1",
        "gene_id": "1524",
        "species": "homo_sapiens",
        "full_name": "C-X3-C motif chemokine receptor 1",
        "aliases": ["CCR1L", "CMKBRV1", "GPR13", "V28"],
        "summary": "Fractalkine receptor highly expressed on homeostatic microglia; essential for neuron-microglia communication and downregulated in DAM.",
        "orthologs": {"human": "CX3CR1", "mouse": "Cx3cr1", "rat": "Cx3cr1"},
        "go_terms": ["GO:0050727", "GO:0005886"],
        "pathways": ["R-HSA-380108", "R-HSA-168256"],
    },
    "P2RY12": {
        "gene_symbol": "P2RY12",
        "gene_id": "64805",
        "species": "homo_sapiens",
        "full_name": "purinergic receptor P2Y12",
        "aliases": ["ADPR-L", "BDPLT8", "HORK3", "P2T(AC)", "P2Y(12)", "P2Y(AC)", "P2Y12", "SP1999"],
        "summary": "Canonical homeostatic microglial purinergic receptor sensing extracellular nucleotides (ATP/ADP); suppressed during microglial disease activation.",
        "orthologs": {"human": "P2RY12", "mouse": "P2ry12", "rat": "P2ry12"},
        "go_terms": ["GO:0005886"],
        "pathways": ["R-HSA-168256"],
    },
    "TMEM119": {
        "gene_symbol": "TMEM119",
        "gene_id": "338773",
        "species": "homo_sapiens",
        "full_name": "transmembrane protein 119",
        "aliases": ["OBIF"],
        "summary": "Specific cell-surface marker for homeostatic brain-resident microglia.",
        "orthologs": {"human": "TMEM119", "mouse": "Tmem119", "rat": "Tmem119"},
        "go_terms": ["GO:0005886"],
        "pathways": ["R-HSA-168256"],
    },
    "CD68": {
        "gene_symbol": "CD68",
        "gene_id": "968",
        "species": "homo_sapiens",
        "full_name": "CD68 molecule",
        "aliases": ["GP110", "LAMP4", "SCARD1"],
        "summary": "Transmembrane glycoprotein heavily expressed in endosomes/lysosomes of macrophages and activated phagocytic microglia.",
        "orthologs": {"human": "CD68", "mouse": "Cd68", "rat": "Cd68"},
        "go_terms": ["GO:0005764", "GO:0006911", "GO:0005044"],
        "pathways": ["R-HSA-2173782"],
    },
}


class BiologicalDBRetriever:
    """
    Biological Database Retriever for Gene Ontology, Reactome Pathways, and NCBI Gene / Orthologs.
    Executes Over-Representation Analysis (ORA) using scipy.stats.hypergeom and Benjamini-Hochberg FDR.
    """

    def __init__(
        self,
        go_database: Optional[List[Dict[str, Any]]] = None,
        pathway_database: Optional[List[Dict[str, Any]]] = None,
        gene_database: Optional[Dict[str, Dict[str, Any]]] = None,
    ):
        self._go_terms = go_database if go_database is not None else CURATED_GO_TERMS
        self._pathways = pathway_database if pathway_database is not None else CURATED_REACTOME_PATHWAYS
        self._genes = gene_database if gene_database is not None else CURATED_NCBI_GENES

    @staticmethod
    def benjamini_hochberg(p_values: List[float]) -> List[float]:
        """
        Calculates Benjamini-Hochberg False Discovery Rate (FDR) adjusted p-values (q-values).
        Enforces monotonic step-up correction capped at 1.0.
        """
        m = len(p_values)
        if m == 0:
            return []
        if m == 1:
            return [min(1.0, float(p_values[0]))]

        # Pair with original index and sort ascending
        indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
        q_values = [1.0] * m

        # Compute raw adjusted: p_(i) * m / i (1-based i)
        raw_adj = []
        for rank, (orig_idx, p_val) in enumerate(indexed_p, start=1):
            adj = (p_val * m) / rank
            raw_adj.append((orig_idx, adj))

        # Enforce step-up monotonicity: min_{j >= i} raw_adj[j]
        current_min = 1.0
        for orig_idx, adj in reversed(raw_adj):
            current_min = min(current_min, adj)
            q_values[orig_idx] = max(0.0, min(1.0, current_min))

        return q_values

    def query_go(
        self,
        gene_list: List[str],
        category: Optional[str] = None,
        background_size: int = 20000,
        fdr_cutoff: float = 0.05,
    ) -> List[GOEnrichment]:
        """
        Runs hypergeometric Over-Representation Analysis (ORA) against Gene Ontology terms.
        
        Args:
            gene_list: Input list of query genes (e.g. top upregulated DEGs).
            category: Optional filter for 'BP', 'MF', or 'CC'.
            background_size: Total number of genes in reference genome (M).
            fdr_cutoff: Maximum FDR q-value threshold.

        Returns:
            List of statistically significant GOEnrichment results sorted by FDR q-value ascending.
        """
        # Canonicalize query genes to uppercase
        clean_genes = set(g.strip().upper() for g in gene_list if g.strip())
        n = len(clean_genes)
        if n == 0:
            return []

        raw_results = []
        p_vals = []

        for item in self._go_terms:
            if category and item["category"].upper() != category.upper():
                continue

            term_genes = set(g.upper() for g in item["genes"])
            matched = list(clean_genes.intersection(term_genes))
            k = len(matched)
            K = item["term_size"]
            M = background_size

            # Hypergeometric survival function: P(X >= k) = sf(k-1, M, K, n)
            if k > 0:
                p_val = float(stats.hypergeom.sf(k - 1, M, K, n))
                # Fold enrichment: (k / n) / (K / M)
                fold_enrich = (k / n) / (K / M) if (n > 0 and K > 0) else 0.0
            else:
                p_val = 1.0
                fold_enrich = 0.0

            raw_results.append({
                "go_id": item["go_id"],
                "term": item["term"],
                "category": item["category"],
                "gene_count": k,
                "term_size": K,
                "background_size": M,
                "sample_size": n,
                "fold_enrichment": round(fold_enrich, 3),
                "p_value": p_val,
                "genes": sorted(matched),
            })
            p_vals.append(p_val)

        # Multiple testing correction
        q_vals = self.benjamini_hochberg(p_vals)

        enriched_list = []
        for res, q_val in zip(raw_results, q_vals):
            if res["gene_count"] > 0 and (fdr_cutoff is None or q_val <= fdr_cutoff):
                enriched_list.append(GOEnrichment(
                    go_id=res["go_id"],
                    term=res["term"],
                    category=res["category"],
                    gene_count=res["gene_count"],
                    term_size=res["term_size"],
                    background_size=res["background_size"],
                    sample_size=res["sample_size"],
                    fold_enrichment=res["fold_enrichment"],
                    p_value=res["p_value"],
                    fdr_q_value=q_val,
                    genes=res["genes"],
                ))

        # Sort by FDR ascending, then p-value ascending
        enriched_list.sort(key=lambda x: (x.fdr_q_value, x.p_value))
        return enriched_list

    def query_pathway(
        self,
        gene_list: List[str],
        source: str = "Reactome",
        background_size: int = 20000,
        fdr_cutoff: float = 0.05,
    ) -> List[PathwayEnrichment]:
        """
        Runs hypergeometric Over-Representation Analysis (ORA) against Reactome pathways.
        """
        clean_genes = set(g.strip().upper() for g in gene_list if g.strip())
        n = len(clean_genes)
        if n == 0:
            return []

        raw_results = []
        p_vals = []

        for pw in self._pathways:
            if source.lower() != "all" and pw.get("source", "Reactome").lower() != source.lower():
                continue

            pw_genes = set(g.upper() for g in pw["genes"])
            matched = list(clean_genes.intersection(pw_genes))
            k = len(matched)
            K = pw["pathway_size"]
            M = background_size

            if k > 0:
                p_val = float(stats.hypergeom.sf(k - 1, M, K, n))
                fold_enrich = (k / n) / (K / M) if (n > 0 and K > 0) else 0.0
            else:
                p_val = 1.0
                fold_enrich = 0.0

            raw_results.append({
                "pathway_id": pw["pathway_id"],
                "pathway_name": pw["pathway_name"],
                "source": pw.get("source", "Reactome"),
                "gene_count": k,
                "pathway_size": K,
                "background_size": M,
                "sample_size": n,
                "fold_enrichment": round(fold_enrich, 3),
                "p_value": p_val,
                "genes": sorted(matched),
            })
            p_vals.append(p_val)

        q_vals = self.benjamini_hochberg(p_vals)

        enriched_list = []
        for res, q_val in zip(raw_results, q_vals):
            if res["gene_count"] > 0 and (fdr_cutoff is None or q_val <= fdr_cutoff):
                enriched_list.append(PathwayEnrichment(
                    pathway_id=res["pathway_id"],
                    pathway_name=res["pathway_name"],
                    source=res["source"],
                    gene_count=res["gene_count"],
                    pathway_size=res["pathway_size"],
                    background_size=res["background_size"],
                    sample_size=res["sample_size"],
                    fold_enrichment=res["fold_enrichment"],
                    p_value=res["p_value"],
                    fdr_q_value=q_val,
                    genes=res["genes"],
                ))

        enriched_list.sort(key=lambda x: (x.fdr_q_value, x.p_value))
        return enriched_list

    def query_gene(self, gene_symbol: str, species: str = "auto") -> Optional[GeneInfo]:
        """
        Retrieves gene metadata, functional summary, and cross-species ortholog mappings from NCBI Gene.
        """
        sym_clean = gene_symbol.strip()
        sym_upper = sym_clean.upper()
        
        # 1. Exact case match first (e.g. 'Trem2' -> mouse, 'TREM2' -> human)
        for key, entry in self._genes.items():
            if entry["gene_symbol"] == sym_clean:
                if species != "auto" and entry.get("species", "").lower() != species.lower():
                    continue
                return GeneInfo(**entry)

        # 2. Case-insensitive symbol match, respecting species preference or standard genetic naming conventions
        # Standard convention: all-caps for human (TREM2), title-case for mouse (Trem2)
        candidates = []
        for key, entry in self._genes.items():
            if entry["gene_symbol"].upper() == sym_upper:
                if species != "auto":
                    if entry.get("species", "").lower() == species.lower():
                        return GeneInfo(**entry)
                else:
                    candidates.append(entry)

        if candidates:
            # If TitleCase input (e.g. Trem2), prefer mouse
            if sym_clean.istitle():
                for c in candidates:
                    if c.get("species") == "mus_musculus":
                        return GeneInfo(**c)
            # If UPPERCASE input (e.g. TREM2), prefer human
            elif sym_clean.isupper():
                for c in candidates:
                    if c.get("species") == "homo_sapiens":
                        return GeneInfo(**c)
            return GeneInfo(**candidates[0])

        # 3. Check aliases
        for key, entry in self._genes.items():
            if sym_upper in [a.upper() for a in entry.get("aliases", [])]:
                if species != "auto" and entry.get("species", "").lower() != species.lower():
                    continue
                return GeneInfo(**entry)

        # 4. Check ortholog mappings
        for key, entry in self._genes.items():
            orthos = entry.get("orthologs", {})
            for sp, orth_sym in orthos.items():
                if orth_sym.upper() == sym_upper:
                    if species != "auto" and entry.get("species", "").lower() != species.lower():
                        continue
                    return GeneInfo(**entry)

        return None

    def resolve_ortholog(self, gene_symbol: str, target_species: str = "human") -> Optional[str]:
        """
        Resolves cross-species ortholog gene symbol (e.g. mouse Apoe -> human APOE).
        """
        info = self.query_gene(gene_symbol)
        if info and target_species.lower() in info.orthologs:
            return info.orthologs[target_species.lower()]
        return None

    def to_evidence_nodes(
        self,
        enrichments: List[Union[GOEnrichment, PathwayEnrichment]],
        task_id: str = "task_pathway_enrichment",
        prior_guided: bool = False,
        hypothesis: Optional[str] = None,
    ) -> List[EvidenceNode]:
        """
        Converts GO or Pathway enrichment results into standardized EvidenceNode instances.
        Enforces epistemic tagging in prior-guided mode.
        """
        nodes = []
        for enr in enrichments:
            is_go = isinstance(enr, GOEnrichment)
            item_id = enr.go_id if is_go else enr.pathway_id
            name = enr.term if is_go else enr.pathway_name
            
            # Calibrate strength based on FDR
            if enr.fdr_q_value < 1e-4:
                strength = EvidenceStrength.VERY_STRONG
            elif enr.fdr_q_value < 0.01:
                strength = EvidenceStrength.STRONG
            else:
                strength = EvidenceStrength.MODERATE

            score = max(0.5, min(1.0, 1.0 - enr.fdr_q_value))
            genes_str = ", ".join(enr.genes[:5]) + ("..." if len(enr.genes) > 5 else "")
            
            summary_base = (
                f"{'GO Term' if is_go else 'Pathway'} '{name}' ({item_id}) significantly enriched "
                f"(FDR={enr.fdr_q_value:.2e}, Fold={enr.fold_enrichment:.2f}x, {enr.gene_count} genes: {genes_str})."
            )

            metrics = {
                "p_value": enr.p_value,
                "fdr_q_value": enr.fdr_q_value,
                "fold_enrichment": enr.fold_enrichment,
                "gene_count": enr.gene_count,
                "matched_genes": enr.genes,
                "prior_guided": prior_guided,
            }

            bio_context: Dict[str, Any] = {
                "item_id": item_id,
                "name": name,
                "type": "GeneOntology" if is_go else "Reactome",
                "matched_genes": enr.genes,
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

            uri = f"go://{item_id}" if is_go else f"reactome://{item_id}"
            
            nodes.append(EvidenceNode(
                evidence_id=f"E_path_{item_id.replace(':', '_').replace('-', '_')}",
                type=EvidenceType.PATHWAY_ENRICHMENT,
                polarity=EvidencePolarity.SUPPORTING,
                strength=strength,
                score=score,
                summary=summary,
                source_task_id=task_id,
                source_artifact_uris=[uri],
                metrics=metrics,
                biological_context=bio_context,
                created_at=datetime.now(timezone.utc),
            ))

        return nodes
