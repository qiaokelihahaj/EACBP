"""Resolve biological identities without promoting specimens or libraries."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import pandas as pd


DONOR_ALIASES = (
    "donor", "donor_id", "subject_id", "subject", "patient_id", "patient",
    "mouse_id", "mouse",
)


@dataclass(frozen=True)
class DonorRole:
    column: str | None
    source: str
    limitation: str

    def as_metrics(self) -> dict[str, Any]:
        return {
            "donor_col": self.column,
            "donor_role_source": self.source,
            "donor_role_limitation": self.limitation,
            "donor_independence_established": False,
        }


def complete_identity(values: pd.Series) -> bool:
    normalized = values.astype("string").str.strip()
    return bool(not values.isna().any() and normalized.ne("").all()
                and not normalized.str.casefold().isin({"nan", "none", "null", "<na>", "nat"}).any())


def resolve_donor_role(obs: pd.DataFrame, params: Mapping[str, Any]) -> DonorRole:
    """An explicit mapping is a recorded user assertion, not proof of independence.

    With an observed biological identity, an explicit specimen/library mapping
    cannot split one individual into multiple supposed replicates, including
    specimens collected under different conditions from the same individual.
    """
    biological_columns = [name for name in DONOR_ALIASES if name in obs.columns]
    explicit = params.get("donor_col")
    limitation = "Biological identity labels do not establish independent biological sampling."
    if explicit is not None:
        column = str(explicit)
        if column not in obs.columns:
            raise ValueError(f"Donor metadata column {column!r} is not present")
        for biological_col in biological_columns:
            if biological_col == column:
                continue
            multiplicity = obs.groupby(biological_col, observed=True, dropna=True)[column].nunique()
            if multiplicity.gt(1).any():
                raise ValueError(
                    f"Explicit donor_col {column!r} splits one biological individual in "
                    f"{biological_col!r} into multiple specimen/library identities; "
                    "these are not independent biological replicates"
                )
        if column not in DONOR_ALIASES:
            limitation = (
                f"User explicitly declared {column!r} as donor identity. Specimen/library "
                "identifiers alone do not establish biological identity or independence; "
                "the sampling design must support this role mapping."
            )
        return DonorRole(column, "user_explicit", limitation)
    if not biological_columns:
        return DonorRole(None, "absent", "Biological donor identity is unavailable; sample/library identifiers are not biological replicates.")
    selected = biological_columns[0]
    for other in biological_columns[1:]:
        if not obs[selected].astype("string").str.strip().equals(obs[other].astype("string").str.strip()):
            raise ValueError(f"Ambiguous donor metadata: {selected!r} and {other!r} disagree; provide donor_col explicitly")
    return DonorRole(selected, "observed_biological_alias", limitation)
