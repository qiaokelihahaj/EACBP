"""
SCData: Single-Cell Data container supporting AnnData interoperability and standalone matrices.
"""

from copy import deepcopy
from typing import Dict, Any, Optional, List, Mapping
import numpy as np
import pandas as pd


class SCData:
    """Lightweight single-cell data representation compatible with Scanpy AnnData."""

    def __init__(
        self,
        X: Any,
        obs: pd.DataFrame,
        var: pd.DataFrame,
        obsm: Optional[Dict[str, np.ndarray]] = None,
        uns: Optional[Dict[str, Any]] = None,
        layers: Optional[Dict[str, Any]] = None,
        raw: Optional[Any] = None,
        obsp: Optional[Dict[str, Any]] = None,
        varm: Optional[Dict[str, Any]] = None,
        varp: Optional[Dict[str, Any]] = None,
    ):
        if hasattr(X, "tocsr"):
            self.X = X.tocsr(copy=True)
        else:
            self.X = np.asarray(X, dtype=np.float32).copy()
        self.obs = obs.copy(deep=True) if isinstance(obs, pd.DataFrame) else pd.DataFrame(obs)
        self.var = var.copy(deep=True) if isinstance(var, pd.DataFrame) else pd.DataFrame(var)
        self.obsm = {k: self._copy_value(v) for k, v in (obsm or {}).items()}
        self.uns = deepcopy(uns) if uns is not None else {}
        self.layers = {k: self._copy_value(v) for k, v in (layers or {}).items()}
        self.raw = self._normalize_raw(raw)
        self.obsp = {k: self._copy_value(v) for k, v in (obsp or {}).items()}
        self.varm = {k: self._copy_value(v) for k, v in (varm or {}).items()}
        self.varp = {k: self._copy_value(v) for k, v in (varp or {}).items()}

        if len(self.obs) != self.X.shape[0]:
            raise ValueError(f"obs length ({len(self.obs)}) does not match X rows ({self.X.shape[0]})")
        if len(self.var) != self.X.shape[1]:
            raise ValueError(f"var length ({len(self.var)}) does not match X cols ({self.X.shape[1]})")
        self._validate_slots()

    @staticmethod
    def _copy_value(value: Any) -> Any:
        """Copy dense, sparse, pandas, and nested slot values safely."""
        if hasattr(value, "copy"):
            try:
                return value.copy(deep=True)
            except TypeError:
                try:
                    return value.copy()
                except TypeError:
                    pass
        return deepcopy(value)

    @classmethod
    def _normalize_raw(cls, raw: Any) -> Optional[Dict[str, Any]]:
        if raw is None:
            return None
        if isinstance(raw, Mapping):
            normalized = deepcopy(dict(raw))
        elif hasattr(raw, "X") and hasattr(raw, "var"):
            normalized = {
                "X": cls._copy_value(raw.X),
                "var": raw.var.copy(deep=True),
                "varm": {
                    key: cls._copy_value(value)
                    for key, value in getattr(raw, "varm", {}).items()
                },
            }
            if hasattr(raw, "uns"):
                normalized["uns"] = deepcopy(dict(raw.uns))
        else:
            raise TypeError(
                "raw must be an AnnData.Raw/AnnData-like object or a dictionary"
            )
        if "X" not in normalized or "var" not in normalized:
            raise ValueError("raw must contain X and var")
        if not isinstance(normalized["var"], pd.DataFrame):
            normalized["var"] = pd.DataFrame(normalized["var"])
        normalized["X"] = cls._copy_value(normalized["X"])
        normalized["varm"] = {
            key: cls._copy_value(value)
            for key, value in normalized.get("varm", {}).items()
        }
        if "uns" in normalized:
            normalized["uns"] = deepcopy(normalized["uns"])
        return normalized

    @staticmethod
    def _slot_shape(value: Any) -> tuple:
        shape = getattr(value, "shape", None)
        if shape is None:
            shape = np.asarray(value).shape
        return tuple(shape)

    def _validate_slots(self) -> None:
        for key, value in self.layers.items():
            if self._slot_shape(value) != self.shape:
                raise ValueError(f"layers[{key!r}] shape does not match X: {self._slot_shape(value)}")
        for key, value in self.obsm.items():
            if len(self._slot_shape(value)) == 0 or self._slot_shape(value)[0] != self.n_obs:
                raise ValueError(f"obsm[{key!r}] first dimension does not match n_obs")
        for key, value in self.obsp.items():
            if self._slot_shape(value) != (self.n_obs, self.n_obs):
                raise ValueError(f"obsp[{key!r}] must be n_obs by n_obs")
        for key, value in self.varm.items():
            if len(self._slot_shape(value)) == 0 or self._slot_shape(value)[0] != self.n_vars:
                raise ValueError(f"varm[{key!r}] first dimension does not match n_vars")
        for key, value in self.varp.items():
            if self._slot_shape(value) != (self.n_vars, self.n_vars):
                raise ValueError(f"varp[{key!r}] must be n_vars by n_vars")
        if self.raw is not None:
            raw_shape = self._slot_shape(self.raw["X"])
            if len(raw_shape) != 2 or raw_shape[0] != self.n_obs:
                raise ValueError("raw.X first dimension does not match n_obs")
            if len(self.raw["var"]) != raw_shape[1]:
                raise ValueError("raw.var length does not match raw.X columns")

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
            obs=self.obs.copy(deep=True),
            var=self.var.copy(deep=True),
            obsm={k: self._copy_value(v) for k, v in self.obsm.items()},
            uns=deepcopy(self.uns),
            layers={k: self._copy_value(v) for k, v in self.layers.items()},
            raw=deepcopy(self.raw),
            obsp={k: self._copy_value(v) for k, v in self.obsp.items()},
            varm={k: self._copy_value(v) for k, v in self.varm.items()},
            varp={k: self._copy_value(v) for k, v in self.varp.items()},
        )

    def subset_obs(self, mask: np.ndarray) -> "SCData":
        mask_array = np.asarray(mask)
        if mask_array.ndim != 1:
            raise ValueError("obs mask must be one-dimensional")
        if mask_array.dtype == bool:
            if len(mask_array) != self.n_obs:
                raise ValueError("boolean obs mask length does not match n_obs")
            positions = np.flatnonzero(mask_array)
        else:
            positions = mask_array.astype(int, copy=False)
            if np.any(positions < 0) or np.any(positions >= self.n_obs):
                raise IndexError("obs subset positions are outside the valid range")

        def subset_rows(value: Any) -> Any:
            if hasattr(value, "iloc") and isinstance(value, pd.DataFrame):
                return value.iloc[positions].copy(deep=True)
            try:
                return self._copy_value(value[positions])
            except (IndexError, TypeError):
                return self._copy_value(np.asarray(value)[positions])

        def subset_pairwise(value: Any) -> Any:
            if hasattr(value, "iloc") and isinstance(value, pd.DataFrame):
                return value.iloc[positions, positions].copy(deep=True)
            if hasattr(value, "tocsr"):
                return value[positions][:, positions].copy()
            array = np.asarray(value)
            return array[np.ix_(positions, positions)].copy()

        subset_raw = None
        if self.raw is not None:
            subset_raw = deepcopy(self.raw)
            subset_raw["X"] = subset_rows(subset_raw["X"])
            if "obs" in subset_raw and isinstance(subset_raw["obs"], pd.DataFrame):
                subset_raw["obs"] = subset_raw["obs"].iloc[positions].copy(deep=True)
        return SCData(
            X=subset_rows(self.X),
            obs=self.obs.iloc[positions].copy(deep=True),
            var=self.var.copy(deep=True),
            obsm={k: subset_rows(v) for k, v in self.obsm.items()},
            uns=deepcopy(self.uns),
            layers={k: subset_rows(v) for k, v in self.layers.items()},
            raw=subset_raw,
            obsp={k: subset_pairwise(v) for k, v in self.obsp.items()},
            varm={k: self._copy_value(v) for k, v in self.varm.items()},
            varp={k: self._copy_value(v) for k, v in self.varp.items()},
        )

    def subset_var(self, mask: np.ndarray) -> "SCData":
        """Subset variables while keeping layers and pairwise slots aligned."""
        mask_array = np.asarray(mask)
        if mask_array.ndim != 1:
            raise ValueError("var mask must be one-dimensional")
        if mask_array.dtype == bool:
            if len(mask_array) != self.n_vars:
                raise ValueError("boolean var mask length does not match n_vars")
            positions = np.flatnonzero(mask_array)
        else:
            positions = mask_array.astype(int, copy=False)
            if np.any(positions < 0) or np.any(positions >= self.n_vars):
                raise IndexError("var subset positions are outside the valid range")

        def subset_cols(value: Any) -> Any:
            if hasattr(value, "iloc") and isinstance(value, pd.DataFrame):
                return value.iloc[:, positions].copy(deep=True)
            if hasattr(value, "tocsr"):
                return value[:, positions].copy()
            return np.asarray(value)[:, positions].copy()

        def subset_rows(value: Any) -> Any:
            if hasattr(value, "iloc") and isinstance(value, pd.DataFrame):
                return value.iloc[positions].copy(deep=True)
            return self._copy_value(value[positions])

        def subset_pairwise(value: Any) -> Any:
            if hasattr(value, "tocsr"):
                return value[positions][:, positions].copy()
            array = np.asarray(value)
            return array[np.ix_(positions, positions)].copy()

        # AnnData.raw is an immutable snapshot in the original variable space;
        # its columns may be a superset or use a different ordering, so a
        # current-X variable subset must leave raw completely intact.
        subset_raw = deepcopy(self.raw)
        return SCData(
            X=self.X[:, positions],
            obs=self.obs.copy(deep=True),
            var=self.var.iloc[positions].copy(deep=True),
            obsm={k: self._copy_value(v) for k, v in self.obsm.items()},
            uns=deepcopy(self.uns),
            layers={k: subset_cols(v) for k, v in self.layers.items()},
            raw=subset_raw,
            obsp={k: self._copy_value(v) for k, v in self.obsp.items()},
            varm={k: subset_rows(v) for k, v in self.varm.items()},
            varp={k: subset_pairwise(v) for k, v in self.varp.items()},
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "X": self._copy_value(self.X),
            "obs": self.obs.copy(deep=True),
            "var": self.var.copy(deep=True),
            "obsm": {k: self._copy_value(v) for k, v in self.obsm.items()},
            "uns": deepcopy(self.uns),
            "layers": {k: self._copy_value(v) for k, v in self.layers.items()},
            "raw": deepcopy(self.raw),
            "obsp": {k: self._copy_value(v) for k, v in self.obsp.items()},
            "varm": {k: self._copy_value(v) for k, v in self.varm.items()},
            "varp": {k: self._copy_value(v) for k, v in self.varp.items()},
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SCData":
        if isinstance(data, cls):
            return data.copy()
        obs = data["obs"] if isinstance(data["obs"], pd.DataFrame) else pd.DataFrame(data["obs"])
        var = data["var"] if isinstance(data["var"], pd.DataFrame) else pd.DataFrame(data["var"])
        return cls(
            X=data["X"],
            obs=obs,
            var=var,
            obsm=data.get("obsm", {}),
            uns=data.get("uns", {}),
            layers=data.get("layers", {}),
            raw=data.get("raw"),
            obsp=data.get("obsp", {}),
            varm=data.get("varm", {}),
            varp=data.get("varp", {}),
        )

    def to_anndata(self):
        """Converts SCData to scanpy/anndata AnnData instance."""
        import anndata as ad
        adata = ad.AnnData(
            X=self._copy_value(self.X),
            obs=self.obs.copy(deep=True),
            var=self.var.copy(deep=True),
            obsm={k: self._copy_value(v) for k, v in self.obsm.items()},
            layers={k: self._copy_value(v) for k, v in self.layers.items()},
            uns=deepcopy(self.uns),
        )
        adata.obsp = {k: self._copy_value(v) for k, v in self.obsp.items()}
        adata.varm = {k: self._copy_value(v) for k, v in self.varm.items()}
        if hasattr(adata, "varp"):
            adata.varp = {k: self._copy_value(v) for k, v in self.varp.items()}
        if self.raw is not None:
            raw_kwargs = {
                "X": self._copy_value(self.raw["X"]),
                "var": self.raw["var"].copy(deep=True),
                "varm": {
                    k: self._copy_value(v)
                    for k, v in self.raw.get("varm", {}).items()
                },
            }
            if "uns" in self.raw:
                raw_kwargs["uns"] = deepcopy(self.raw["uns"])
            adata.raw = ad.AnnData(**raw_kwargs)
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
            rng = np.random.default_rng(random_seed)
            sub_idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
            sub_idx.sort()
            adata = adata[sub_idx].copy()

        # Handle sparse or dense matrix
        X = adata.X
        if hasattr(X, "tocsr"):
            X = X.tocsr(copy=True)
            n_counts = np.asarray(X.sum(axis=1)).flatten()
            n_genes = np.asarray((X > 0).sum(axis=1)).flatten()
            n_cells_per_gene = np.asarray((X > 0).sum(axis=0)).flatten()
        else:
            X = np.asarray(X, dtype=np.float32).copy()
            n_counts = np.sum(X, axis=1)
            n_genes = np.sum(X > 0, axis=1)
            n_cells_per_gene = np.sum(X > 0, axis=0)

        obs = adata.obs.copy()
        var = adata.var.copy()

        # Ensure essential columns exist in obs
        if "cell_id" not in obs.columns:
            obs["cell_id"] = obs.index.astype(str)
        if "n_counts" not in obs.columns:
            obs["n_counts"] = n_counts
        if "n_genes" not in obs.columns:
            obs["n_genes"] = n_genes

        # Ensure essential columns exist in var
        if "gene_name" not in var.columns:
            var["gene_name"] = var.index.astype(str)
        if "n_cells" not in var.columns:
            var["n_cells"] = n_cells_per_gene

        obsm_dict = {}
        for k in adata.obsm.keys():
            obsm_dict[k] = cls._copy_value(adata.obsm[k])

        uns_dict = deepcopy(dict(adata.uns)) if hasattr(adata, "uns") else {}
        layers_dict = {
            k: cls._copy_value(adata.layers[k]) for k in adata.layers.keys()
            if k is not None  # AnnData 0.13 exposes X as the None layer.
        }
        obsp_dict = {
            k: cls._copy_value(adata.obsp[k]) for k in adata.obsp.keys()
        }
        varm_dict = {
            k: cls._copy_value(adata.varm[k]) for k in adata.varm.keys()
        }
        varp_dict = {
            k: cls._copy_value(adata.varp[k]) for k in adata.varp.keys()
        } if hasattr(adata, "varp") else {}
        raw_dict = cls._normalize_raw(adata.raw) if getattr(adata, "raw", None) is not None else None

        return cls(
            X=X,
            obs=obs,
            var=var,
            obsm=obsm_dict,
            uns=uns_dict,
            layers=layers_dict,
            raw=raw_dict,
            obsp=obsp_dict,
            varm=varm_dict,
            varp=varp_dict,
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

        return cls(
            X=X,
            obs=obs,
            var=var,
            uns={"is_simulated": True, "data_origin": "synthetic"},
        )

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

        return cls(
            X=X,
            obs=obs,
            var=var,
            uns={"is_simulated": True, "data_origin": "synthetic"},
        )
