from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.schemas.evidence import EvidenceNode, EvidenceType


def node(name, score, roots, kind=EvidenceType.PSEUDOBULK_DEG):
    return EvidenceNode(evidence_id=name, type=kind, score=score, summary=name,
                        source_task_id=name, data_origin_uris=roots)


def test_more_same_matrix_methods_do_not_inflate_dimension():
    first = node("a", .5, ["adata://s/raw/v1"])
    stronger = node("b", .99, ["adata://s/raw/v1"])
    original = ConfidenceCalculator.calculate([first], [])
    repeated = ConfidenceCalculator.calculate([first, stronger, stronger], [])
    assert repeated.association == original.association
    assert repeated.overall == original.overall


def test_shared_roots_are_transitive_and_order_independent():
    a = node("a", .8, ["r1"])
    b = node("b", .6, ["r2"])
    bridge = node("c", .9, ["r1", "r2"])
    assert ConfidenceCalculator.calculate([a, b, bridge], []) == ConfidenceCalculator.calculate([bridge, b, a], [])


def test_sensitivity_does_not_add_independent_support():
    a = node("a", .7, ["r1"])
    sensitivity = node("loo", 1, ["r1"], EvidenceType.SENSITIVITY_ANALYSIS)
    assert ConfidenceCalculator.calculate([a], []) == ConfidenceCalculator.calculate([a, sensitivity], [])
