"""
SCData: Single-Cell Data container supporting AnnData interoperability and standalone matrices.
"""

from typing import Dict, Any, Optional, List
import numpy as np
import pandas as pd


class SCData:
    """Lightweight single-cell data representation compatible with Scanpy AnnData."""

    def __init__(
        self,
        X: np.ndarray,
        obs: pd.DataFrame,
        var: pd.DataFrame,
        obsm: Optional[Dict[str, np.ndarray]] = None,
        uns: Optional[Dict[str, Any]] = None,
    ):
        self.X = np.asarray(X, dtype=np.float32)
        self.obs = obs.copy()
        self.var = var.copy()
        self.obsm = obsm or {}
        self.uns = uns or {}

        if len(self.obs) != self.X.shape[0]:
            raise ValueError(f"obs length ({len(self.obs)}) does not match X rows ({self.X.shape[0]})")
        if len(self.var) != self.X.shape[1]:
            raise ValueError(f"var length ({len(self.var)}) does not match X cols ({self.X.shape[1]})")

    @property
    def n_obs(self) -> int:
        return self.X.shape[0]

    @property
    def n_vars(self) -> int:
        return self.X.shape[1]

    @property
    def shape(self):
        return self.X.shape

    def copy(self) -> "SCData":
        return SCData(
            X=self.X.copy(),
            obs=self.obs.copy(),
            var=self.var.copy(),
            obsm={k: v.copy() for k, v in self.obsm.items()},
            uns={k: v for k, v in self.uns.items()},
        )

    def subset_obs(self, mask: np.ndarray) -> "SCData":
        return SCData(
            X=self.X[mask],
            obs=self.obs.iloc[mask].copy().reset_index(drop=True),
            var=self.var.copy(),
            obsm={k: v[mask] for k, v in self.obsm.items()},
            uns=self.uns.copy(),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "X": self.X,
            "obs": self.obs,
            "var": self.var,
            "obsm": self.obsm,
            "uns": self.uns,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SCData":
        obs = data["obs"] if isinstance(data["obs"], pd.DataFrame) else pd.DataFrame(data["obs"])
        var = data["var"] if isinstance(data["var"], pd.DataFrame) else pd.DataFrame(data["var"])
        return cls(
            X=data["X"],
            obs=obs,
            var=var,
            obsm=data.get("obsm", {}),
            uns=data.get("uns", {}),
        )

    def to_anndata(self):
        """Converts SCData to scanpy/anndata AnnData instance."""
        import anndata as ad
        adata = ad.AnnData(
            X=self.X,
            obs=self.obs,
            var=self.var,
            obsm=self.obsm,
            uns=self.uns,
        )
        return adata

    @classmethod
    def from_anndata(
        cls,
        adata: Any,
        max_cells: Optional[int] = None,
        random_seed: int = 42,
    ) -> "SCData":
        """Builds SCData directly from an anndata.AnnData object."""
        if max_cells is not None and adata.n_obs > max_cells:
            np.random.seed(random_seed)
            sub_idx = np.random.choice(adata.n_obs, size=max_cells, replace=False)
            sub_idx.sort()
            adata = adata[sub_idx].copy()

        # Handle sparse or dense matrix
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)

        obs = adata.obs.copy()
        var = adata.var.copy()

        # Ensure essential columns exist in obs
        if "cell_id" not in obs.columns:
            obs["cell_id"] = obs.index.astype(str)
        if "n_counts" not in obs.columns:
            obs["n_counts"] = np.sum(X, axis=1)
        if "n_genes" not in obs.columns:
            obs["n_genes"] = np.sum(X > 0, axis=1)

        # Ensure essential columns exist in var
        if "gene_name" not in var.columns:
            var["gene_name"] = var.index.astype(str)
        if "n_cells" not in var.columns:
            var["n_cells"] = np.sum(X > 0, axis=0)

        obsm_dict = {}
        for k in adata.obsm.keys():
            obsm_dict[k] = np.asarray(adata.obsm[k])

        uns_dict = dict(adata.uns) if hasattr(adata, "uns") else {}

        return cls(
            X=X,
            obs=obs,
            var=var,
            obsm=obsm_dict,
            uns=uns_dict,
        )

    @classmethod
    def from_h5ad(
        cls,
        file_path: str,
        max_cells: Optional[int] = None,
        random_seed: int = 42,
    ) -> "SCData":
        """Loads single-cell data directly from a .h5ad file."""
        import scanpy as sc
        adata = sc.read_h5ad(file_path)
        return cls.from_anndata(adata, max_cells=max_cells, random_seed=random_seed)

    @classmethod
    def create_synthetic_ad_study(
        cls,
        n_cells: int = 1200,
        n_genes: int = 500,
        n_ad_mice: int = 6,
        n_ctrl_mice: int = 6,
        random_seed: int = 42,
    ) -> "SCData":
        """Generates a realistic synthetic AD vs Control single-cell dataset with known microglia subpopulation."""
        np.random.seed(random_seed)

        # Biological units (mice)
        ad_mice = [f"AD_mouse_{i+1:02d}" for i in range(n_ad_mice)]
        ctrl_mice = [f"Ctrl_mouse_{i+1:02d}" for i in range(n_ctrl_mice)]
        all_mice = ad_mice + ctrl_mice

        # Assign cells to mice & batches
        assigned_mice = np.random.choice(all_mice, size=n_cells)
        conditions = ["AD" if "AD" in m else "control" for m in assigned_mice]
        batches = [f"batch_{(i % 2) + 1}" for i, m in enumerate(assigned_mice)]

        # Cell types: Microglia (~35%), Astrocytes (~25%), Neurons (~25%), Oligodendrocytes (~15%)
        cell_type_probs = [0.35, 0.25, 0.25, 0.15]
        cell_types = np.random.choice(["Microglia", "Astrocytes", "Neurons", "Oligodendrocytes"], size=n_cells, p=cell_type_probs)

        # Genes
        gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
        # Add key marker genes
        marker_map = {
            0: "Apoe",
            1: "Trem2",
            2: "Clec7a",
            3: "Cx3cr1",
            4: "P2ry12",
            5: "Gfap",
            6: "Rbfox3",
            7: "Mog",
            8: "Tmem119",
            9: "Itgax",
        }
        for idx, g in marker_map.items():
            gene_names[idx] = g

        # Expression base matrix (Poisson / Negative Binomial counts)
        base_rate = np.random.gamma(2.0, 1.0, size=(n_cells, n_genes))
        
        # Inject cell-type specific markers
        for i in range(n_cells):
            ct = cell_types[i]
            cond = conditions[i]
            if ct == "Microglia":
                base_rate[i, 3] *= 4.0  # Cx3cr1
                base_rate[i, 4] *= 3.5  # P2ry12
                base_rate[i, 8] *= 3.5  # Tmem119
                # If AD condition, induce disease-associated microglia (DAM) signature in a subset
                if cond == "AD" and np.random.rand() > 0.4:
                    base_rate[i, 0] *= 5.0  # Apoe high
                    base_rate[i, 1] *= 4.0  # Trem2 high
                    base_rate[i, 2] *= 4.5  # Clec7a high
                    base_rate[i, 4] *= 0.2  # P2ry12 downregulated
            elif ct == "Astrocytes":
                base_rate[i, 5] *= 6.0  # Gfap
            elif ct == "Neurons":
                base_rate[i, 6] *= 6.0  # Rbfox3
            elif ct == "Oligodendrocytes":
                base_rate[i, 7] *= 6.0  # Mog

        # Add batch effect to batch_2
        for i in range(n_cells):
            if batches[i] == "batch_2":
                base_rate[i, :50] *= 1.4

        # Sample raw integer counts
        X = np.random.poisson(base_rate).astype(np.float32)

        # QC stats: mito %
        mito_pct = np.random.beta(2, 25, size=n_cells) * 100.0

        obs = pd.DataFrame({
            "cell_id": [f"cell_{i:05d}" for i in range(n_cells)],
            "mouse_id": assigned_mice,
            "condition": conditions,
            "batch": batches,
            "cell_type_ground_truth": cell_types,
            "percent_mito": mito_pct,
            "n_counts": X.sum(axis=1),
            "n_genes": (X > 0).sum(axis=1),
        })

        var = pd.DataFrame({
            "gene_name": gene_names,
            "n_cells": (X > 0).sum(axis=0),
        })

        return cls(X=X, obs=obs, var=var)

    @classmethod
    def create_synthetic_kat8_study(
        cls,
        n_cells: int = 1200,
        n_genes: int = 500,
        n_cko_mice: int = 4,
        n_con_mice: int = 4,
        random_seed: int = 42,
    ) -> "SCData":
        """
        Generates a realistic synthetic Kat8 cKO vs Control single-cell study.
        Models Kat8 (Mof) knockout in neural/germline/immune lineage with H4K16ac transcriptional disruption.
        """
        np.random.seed(random_seed)

        cko_mice = [f"P12_cKO_mouse_{i+1:02d}" for i in range(n_cko_mice)]
        con_mice = [f"P12_con_mouse_{i+1:02d}" for i in range(n_con_mice)]
        all_mice = cko_mice + con_mice

        assigned_mice = np.random.choice(all_mice, size=n_cells)
        conditions = ["cKO" if "cKO" in m else "con" for m in assigned_mice]
        batches = [f"batch_{(i % 2) + 1}" for i, m in enumerate(assigned_mice)]

        # Cell types in P12 tissue (e.g. Progenitors, Differentiating Neurons, Mature Glia, Astrocytes)
        cell_types_pool = ["Progenitors", "Immature_Neurons", "Mature_Neurons", "Astrocytes"]
        cell_type_probs = [0.30, 0.30, 0.25, 0.15]
        cell_types = np.random.choice(cell_types_pool, size=n_cells, p=cell_type_probs)

        # Genes
        gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
        marker_map = {
            0: "Kat8",      # Target gene (Mof)
            1: "Kansl1",    # NSL complex subunit
            2: "Msl1",      # MSL complex subunit
            3: "H4c1",      # Histone H4
            4: "Cdk1",      # Cell cycle kinase (down in cKO)
            5: "Top2a",     # Cell proliferation marker
            6: "Bax",       # Pro-apoptotic factor (up in cKO)
            7: "Cdkn1a",    # p21 cell cycle arrest (up in cKO)
            8: "Sox2",      # Neural stem/progenitor marker
            9: "Dcx",       # Immature neuronal marker
            10: "Rbfox3",   # NeuN mature neuron
            11: "Gfap",     # Astrocyte marker
            12: "Trp53",    # p53 stress pathway
            13: "H2ax",     # DNA damage marker
        }
        for idx, g in marker_map.items():
            gene_names[idx] = g

        base_rate = np.random.gamma(2.5, 1.2, size=(n_cells, n_genes))

        # Cell-type specific baselines
        for i in range(n_cells):
            ct = cell_types[i]
            cond = conditions[i]

            if ct == "Progenitors":
                base_rate[i, 8] *= 5.0  # Sox2
                base_rate[i, 4] *= 3.5  # Cdk1
                base_rate[i, 5] *= 3.5  # Top2a
            elif ct == "Immature_Neurons":
                base_rate[i, 9] *= 6.0  # Dcx
            elif ct == "Mature_Neurons":
                base_rate[i, 10] *= 5.5 # Rbfox3
            elif ct == "Astrocytes":
                base_rate[i, 11] *= 6.0 # Gfap

            # Kat8 cKO effects:
            if cond == "cKO":
                base_rate[i, 0] *= 0.05 # Kat8 strongly knocked out
                base_rate[i, 6] *= 3.2  # Bax upregulated (apoptosis stress)
                base_rate[i, 7] *= 3.8  # Cdkn1a upregulated (cell cycle arrest)
                base_rate[i, 12] *= 2.5 # Trp53 upregulated
                base_rate[i, 13] *= 2.8 # H2ax upregulated (chromatin / DNA stress)
                # Cell proliferation and histone expression downregulated
                base_rate[i, 4] *= 0.3  # Cdk1 down
                base_rate[i, 5] *= 0.3  # Top2a down
                base_rate[i, 3] *= 0.4  # H4c1 transcription altered
            else:
                base_rate[i, 0] *= 3.0  # Kat8 normally expressed in control

        X = np.random.poisson(base_rate).astype(np.float32)
        mito_pct = np.random.beta(2, 28, size=n_cells) * 100.0

        obs = pd.DataFrame({
            "cell_id": [f"cell_{i:05d}" for i in range(n_cells)],
            "mouse_id": assigned_mice,
            "condition": conditions,
            "batch": batches,
            "cell_type_ground_truth": cell_types,
            "percent_mito": mito_pct,
            "n_counts": X.sum(axis=1),
            "n_genes": (X > 0).sum(axis=1),
        })

        var = pd.DataFrame({
            "gene_name": gene_names,
            "n_cells": (X > 0).sum(axis=0),
        })

        return cls(X=X, obs=obs, var=var)
