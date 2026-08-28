"""
4-Tier Scientific Language Protocol enforcer and formatter.
"""

import re
from typing import Tuple, Optional
from eacbp.schemas.evidence import LanguageTier, ClaimNode, ConfidenceScore


FORBIDDEN_CAUSAL_VERBS_OBSERVATIONAL = [
    r"\bdrives\b",
    r"\bdrive\b",
    r"\bcauses\b",
    r"\bcause\b",
    r"\bproves\b",
    r"\bdemonstrates causality\b",
    r"\bis responsible for\b",
]


class LanguageEnforcer:
    """Enforces strict epistemic calibration across the 4 Language Tiers."""

    @staticmethod
    def audit_statement(statement: str, tier: LanguageTier, causal_status: str) -> Tuple[bool, Optional[str]]:
        """Audits whether a statement overclaims causality in observational studies."""
        if causal_status == "observational" and tier in (LanguageTier.LEVEL_1_OBSERVATION, LanguageTier.LEVEL_2_STATISTICAL_INFERENCE):
            for pattern in FORBIDDEN_CAUSAL_VERBS_OBSERVATIONAL:
                if re.search(pattern, statement, re.IGNORECASE):
                    return False, (
                        f"EPISTEMIC VIOLATION: Statement contains strong causal verb matching '{pattern}' "
                        f"in an observational study under {tier.value}. Use associative or qualified terminology "
                        f"(e.g. 'is associated with', 'correlates with', 'may participate in')."
                    )
        return True, None

    @staticmethod
    def format_claim_sentence(
        observation: str,
        stat_inference: Optional[str] = None,
        interpretation: Optional[str] = None,
        hypothesis: Optional[str] = None,
    ) -> str:
        """Formats a scientific claim spanning the four calibrated tiers."""
        parts = [f"**[{LanguageTier.LEVEL_1_OBSERVATION.value}]**: {observation}"]
        if stat_inference:
            parts.append(f"**[{LanguageTier.LEVEL_2_STATISTICAL_INFERENCE.value}]**: {stat_inference}")
        if interpretation:
            parts.append(f"**[{LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION.value}]**: {interpretation}")
        if hypothesis:
            parts.append(f"**[{LanguageTier.LEVEL_4_HYPOTHESIS.value}]**: {hypothesis}")
        return "\n\n".join(parts)
