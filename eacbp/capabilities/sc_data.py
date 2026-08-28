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
