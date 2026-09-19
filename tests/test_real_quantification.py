"""STARsolo process-boundary tests; generated matrices are fixtures, not aligned reads."""
from pathlib import Path
import subprocess
import numpy as np
import pytest
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.quantification import FASTQQuantificationCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def setup_star(tmp_path):
    reads = []
    for lane in (1, 2):
        pair = []
        for read in (1, 2):
            path = tmp_path / f"sample_L{lane}_R{read}.fastq"
            path.write_text("@r\nACGT\n+\nIIII\n")
            pair.append(str(path))
        reads.append(pair)
    genome = tmp_path / "genome"
    genome.mkdir()
    (genome / "Genome").write_bytes(b"test reference")
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("AAAA\nBBBB\n")
    reg = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "fastq://star/reads/v1"
    reg.register(uri, {"samples": {"s1": {"R1": [p[0] for p in reads], "R2": [p[1] for p in reads],
                       "metadata": {"donor": "donor1", "condition": "treated", "batch": "batch1"}}}}, ArtifactType.FASTQ, "star", "input", "input")
    task = TaskContract(task_id="quant", capability="quantification", method="starsolo_v1", input_artifacts=[uri],
                        expected_outputs=["adata://star/raw/v1"], parameters={"star_bin": "STAR", "genome_dir": str(genome),
                        "whitelist_path": str(whitelist), "chemistry": "10xv3", "work_dir": str(tmp_path / "work")})
    return reg, task, reads


def write_matrix(command, folder="filtered"):
    prefix = Path(command[command.index("--outFileNamePrefix") + 1])
    directory = prefix / "Solo.out" / "Gene" / folder
    directory.mkdir(parents=True)
    mmwrite(str(directory / "matrix.mtx"), csr_matrix(np.array([[1, 2], [0, 3], [4, 1]])))
    (directory / "features.tsv").write_text("ENSG1\tA\tGene Expression\nENSG2\tA\tGene Expression\nENSG3\tC\tGene Expression\n")
    (directory / "barcodes.tsv").write_text("AAAA\nBBBB\n")


def test_starsolo_multilane_command_matrix_and_hashes(tmp_path, monkeypatch):
    reg, task, reads = setup_star(tmp_path)
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        write_matrix(command)
        return subprocess.CompletedProcess(command, 0, "", "")
    monkeypatch.setattr("eacbp.capabilities.quantification.subprocess.run", run)
    result = FASTQQuantificationCapability(implementation_id="starsolo_v1").execute(task, reg)
    assert result.status == TaskStatus.SUCCESS, result.error_message
    command = commands[0]
    index = command.index("--readFilesIn")
    assert command[index + 1:index + 3] == [",".join(p[1] for p in reads), ",".join(p[0] for p in reads)]
    assert command[command.index("--soloUMIlen") + 1] == "12"
    _, data = reg.get(result.output_artifacts[0])
    assert data.shape == (2, 3) and list(data.var.index) == ["ENSG1", "ENSG2", "ENSG3"]
    assert data.obs.cell_id.tolist() == ["s1_AAAA", "s1_BBBB"]
    assert data.obs.donor.tolist() == ["donor1", "donor1"]
    assert result.metrics["input_file_hashes"] or result.metrics["input_hashes"]


@pytest.mark.parametrize("failure", ["exit", "only_raw", "changed_input"])
def test_starsolo_failure_never_publishes_counts(tmp_path, monkeypatch, failure):
    reg, task, reads = setup_star(tmp_path)
    def run(command, **kwargs):
        if failure == "exit":
            raise subprocess.CalledProcessError(1, command, stderr="STAR failed")
        write_matrix(command, "raw" if failure == "only_raw" else "filtered")
        if failure == "changed_input":
            Path(reads[0][0]).write_text("changed")
        return subprocess.CompletedProcess(command, 0, "", "")
    monkeypatch.setattr("eacbp.capabilities.quantification.subprocess.run", run)
    result = FASTQQuantificationCapability(implementation_id="starsolo_v1").execute(task, reg)
    assert result.status != TaskStatus.SUCCESS
    assert not reg.exists("adata://star/raw/v1")
