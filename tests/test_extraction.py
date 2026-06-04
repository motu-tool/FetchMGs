import gzip
import pathlib
import pytest
from fetchmgs.fetchmgs import extraction_genes, extraction_genomes, load_input_files_from_file

DATA_DIR = pathlib.Path(__file__).parent / "data"


def read_gz(path):
    with gzip.open(path, "rt") as f:
        return f.read()


def assert_dir_matches(tmp_path, expected_dir):
    """Compare every .gz file in expected_dir against its decompressed counterpart in tmp_path."""
    expected_files = sorted(expected_dir.glob("*.gz"))
    assert expected_files, f"No expected files found in {expected_dir}"
    for expected_gz in expected_files:
        actual_name = expected_gz.name[:-3]  # strip trailing .gz
        actual_file = tmp_path / actual_name
        assert actual_file.exists(), f"Missing output file: {actual_name}"
        assert actual_file.read_text() == read_gz(expected_gz), f"Mismatch in {actual_name}"


def _single_cases(data_dir):
    """Return (input_file, expected_dir) pairs from per-sample output directories."""
    cases = []
    for output_dir in sorted(
        d for d in data_dir.iterdir()
        if d.is_dir() and d.name.startswith("output_") and d.name != "output_multi"
    ):
        for gz_file in output_dir.glob("*.fetchMGs.*.gz"):
            input_name = gz_file.name.split(".fetchMGs.")[0]
            input_file = data_dir / "input" / input_name
            if input_file.exists():
                cases.append((input_file, output_dir))
                break
    return cases


# ── genome mode ───────────────────────────────────────────────────────────────

_GENOME_DIR = DATA_DIR / "genomes"
_GENOME_SINGLE = _single_cases(_GENOME_DIR)


@pytest.mark.parametrize("input_fa,expected_dir", _GENOME_SINGLE,
                         ids=[c[1].name for c in _GENOME_SINGLE])
def test_genome_single(input_fa, expected_dir, tmp_path):
    extraction_genomes([input_fa], tmp_path, mode="genome", threads=1, very_best=True)
    assert_dir_matches(tmp_path, expected_dir)


def test_genome_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_GENOME_DIR)
    input_files = load_input_files_from_file(_GENOME_DIR / "input" / "genomes")
    extraction_genomes(input_files, tmp_path, mode="genome", threads=1, very_best=True)
    assert_dir_matches(tmp_path, _GENOME_DIR / "output_multi")


# ── gene mode (pre-called genome genes) ──────────────────────────────────────

_GENOME_GENES_DIR = DATA_DIR / "genomes_genes"
_GENOME_GENES_SINGLE = _single_cases(_GENOME_GENES_DIR)


@pytest.mark.parametrize("input_faa,expected_dir", _GENOME_GENES_SINGLE,
                         ids=[c[1].name for c in _GENOME_GENES_SINGLE])
def test_gene_from_genomes_single(input_faa, expected_dir, tmp_path):
    input_fna = input_faa.with_name(input_faa.name.replace(".faa.gz", ".fna.gz"))
    extraction_genes([input_faa], [input_fna], tmp_path, threads=1, very_best=True)
    assert_dir_matches(tmp_path, expected_dir)


def test_gene_from_genomes_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_GENOME_GENES_DIR)
    input_faa_files = load_input_files_from_file(_GENOME_GENES_DIR / "input" / "genome_genes_faa")
    input_fna_files = load_input_files_from_file(_GENOME_GENES_DIR / "input" / "genome_genes_fna")
    extraction_genes(input_faa_files, input_fna_files, tmp_path, threads=1, very_best=True)
    assert_dir_matches(tmp_path, _GENOME_GENES_DIR / "output_multi")


# ── metagenome mode ───────────────────────────────────────────────────────────

_METAGENOME_DIR = DATA_DIR / "metagenome"
_METAGENOME_SINGLE = _single_cases(_METAGENOME_DIR)


@pytest.mark.parametrize("input_fa,expected_dir", _METAGENOME_SINGLE,
                         ids=[c[1].name for c in _METAGENOME_SINGLE])
def test_metagenome_single(input_fa, expected_dir, tmp_path):
    extraction_genomes([input_fa], tmp_path, mode="metagenome", threads=1, very_best=False)
    assert_dir_matches(tmp_path, expected_dir)


def test_metagenome_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_METAGENOME_DIR)
    input_files = load_input_files_from_file(_METAGENOME_DIR / "input" / "metagenomes")
    extraction_genomes(input_files, tmp_path, mode="metagenome", threads=1, very_best=False)
    assert_dir_matches(tmp_path, _METAGENOME_DIR / "output_multi")


# ── gene mode (pre-called metagenome genes) ───────────────────────────────────

_METAGENOME_GENES_DIR = DATA_DIR / "metagenome_genes"
_METAGENOME_GENES_SINGLE = _single_cases(_METAGENOME_GENES_DIR)


@pytest.mark.parametrize("input_faa,expected_dir", _METAGENOME_GENES_SINGLE,
                         ids=[c[1].name for c in _METAGENOME_GENES_SINGLE])
def test_gene_from_metagenomes_single(input_faa, expected_dir, tmp_path):
    input_fna = input_faa.with_name(input_faa.name.replace(".faa.gz", ".fna.gz"))
    extraction_genes([input_faa], [input_fna], tmp_path, threads=1, very_best=False)
    assert_dir_matches(tmp_path, expected_dir)


def test_gene_from_metagenomes_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_METAGENOME_GENES_DIR)
    input_faa_files = load_input_files_from_file(_METAGENOME_GENES_DIR / "input" / "metagenome_genes_faa")
    input_fna_files = load_input_files_from_file(_METAGENOME_GENES_DIR / "input" / "metagenome_genes_fna")
    extraction_genes(input_faa_files, input_fna_files, tmp_path, threads=1, very_best=False)
    assert_dir_matches(tmp_path, _METAGENOME_GENES_DIR / "output_multi")
