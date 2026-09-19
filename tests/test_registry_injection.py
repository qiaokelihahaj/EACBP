import pytest
from eacbp.capabilities import CapabilityRegistry, create_default_capability_registry
from eacbp.adapters import ChatCellAdapter
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry


def test_explicit_registry_is_preserved_without_default_replacement(tmp_path):
    registry = CapabilityRegistry()
    custom = ChatCellAdapter()
    registry.register(custom)
    orchestrator = ScientificOrchestrator(ArtifactRegistry(str(tmp_path)), registry)
    assert orchestrator.capability_registry is registry
    assert registry.get(custom.capability_name, custom.implementation_id) is custom
    assert list(registry.list_capabilities()) == [custom.capability_name]


def test_duplicate_registration_requires_explicit_overwrite():
    registry = CapabilityRegistry()
    original, replacement = ChatCellAdapter(), ChatCellAdapter()
    registry.register(original)
    with pytest.raises(ValueError):
        registry.register(replacement)
    assert registry.get(original.capability_name, original.implementation_id) is original
    registry.register(replacement, overwrite=True)
    assert registry.get(replacement.capability_name, replacement.implementation_id) is replacement


def test_default_factory_includes_optional_planes_without_orchestrator():
    registered = create_default_capability_registry().list_capabilities()
    assert {"knowledge_retrieval", "cell_cell_communication", "genetic_perturbation_simulation", "chatcell_dialogue_prediction"} <= registered.keys()
