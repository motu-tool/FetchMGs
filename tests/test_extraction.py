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
    extraction_genes([input_faa], [input_fna], tmp_path, threads=1, very_best=True, chunk_size=100)
    assert_dir_matches(tmp_path, expected_dir)


def test_gene_from_genomes_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_GENOME_GENES_DIR)
    input_faa_files = load_input_files_from_file(_GENOME_GENES_DIR / "input" / "genome_genes_faa")
    input_fna_files = load_input_files_from_file(_GENOME_GENES_DIR / "input" / "genome_genes_fna")
    extraction_genes(input_faa_files, input_fna_files, tmp_path, threads=1, very_best=True, chunk_size=100)
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
    extraction_genes([input_faa], [input_fna], tmp_path, threads=1, very_best=False, chunk_size=100)
    assert_dir_matches(tmp_path, expected_dir)


def test_gene_from_metagenomes_multi(tmp_path, monkeypatch):
    monkeypatch.chdir(_METAGENOME_GENES_DIR)
    input_faa_files = load_input_files_from_file(_METAGENOME_GENES_DIR / "input" / "metagenome_genes_faa")
    input_fna_files = load_input_files_from_file(_METAGENOME_GENES_DIR / "input" / "metagenome_genes_fna")
    extraction_genes(input_faa_files, input_fna_files, tmp_path, threads=1, very_best=False, chunk_size=100)
    assert_dir_matches(tmp_path, _METAGENOME_GENES_DIR / "output_multi")


# ── load_input_files_from_file ────────────────────────────────────────────────

def test_map_file_with_empty_lines(tmp_path):
    input_file = _GENOME_DIR / "input" / "ACIN21-1_SAMN05422111_MAG_00000002.fa.gz"
    map_file = tmp_path / "map.txt"
    map_file.write_text(f"{input_file}\n\n")
    result = load_input_files_from_file(map_file)
    assert result == [input_file]


# ── 2.2.0: threads, batches across files, sequence types, empty input, ties ──

from fetchmgs.fetchmgs import classify_sequences, check_sequence_type, main


def write_gz(path, text):
    with gzip.open(path, "wt") as f:
        f.write(text)


def parse_fasta(text):
    records = []
    for line in text.splitlines():
        if line.startswith(">"):
            records.append([line, ""])
        else:
            records[-1][1] += line
    return records


def read_fasta_gz(path):
    return parse_fasta(read_gz(path))


def test_genome_multi_threads(tmp_path, monkeypatch):
    monkeypatch.chdir(_GENOME_DIR)
    input_files = load_input_files_from_file(_GENOME_DIR / "input" / "genomes")
    extraction_genomes(input_files, tmp_path, mode="genome", threads=4, very_best=True)
    assert_dir_matches(tmp_path, _GENOME_DIR / "output_multi")


def test_metagenome_multi_threads(tmp_path, monkeypatch):
    monkeypatch.chdir(_METAGENOME_DIR)
    input_files = load_input_files_from_file(_METAGENOME_DIR / "input" / "metagenomes")
    extraction_genomes(input_files, tmp_path, mode="metagenome", threads=4, very_best=False)
    assert_dir_matches(tmp_path, _METAGENOME_DIR / "output_multi")


@pytest.mark.parametrize("chunk_size", [7, 1000, 500000])
def test_gene_batches_across_files(chunk_size, tmp_path):
    """Two copies of one gene set (same gene names) in one run: batches cross the file boundary, hits must go back
    to the right file, and each output must equal the single-file one."""
    name = "ACIN21-1_SAMN05422152_METAG.scaffolds.min500.sub.genes"
    source = _METAGENOME_GENES_DIR / "input"
    expected_dir = _METAGENOME_GENES_DIR / f"output_{name}"
    inputs = tmp_path / "in"
    inputs.mkdir()
    faa, fna = [], []
    for copy in ("a", "b"):
        for kind, files in (("faa", faa), ("fna", fna)):
            target = inputs / f"{copy}.{kind}.gz"
            target.write_bytes((source / f"{name}.{kind}.gz").read_bytes())
            files.append(target)
    out = tmp_path / "out"
    extraction_genes(faa, fna, out, threads=2, very_best=False, chunk_size=chunk_size)
    for copy in ("a", "b"):
        for kind in ("faa", "fna", "scores"):
            expected = read_gz(expected_dir / f"{name}.faa.gz.fetchMGs.{kind}.gz")
            assert (out / f"{copy}.faa.gz.fetchMGs.{kind}").read_text() == expected, f"{copy} {kind}"


@pytest.mark.parametrize("sequences,expected", [
    (["ACGTACGTNNACGT", "acgtu"], "NT"),
    (["ACGTRYKMACGTACGTACGTACGTACGTACGTACGTACGTSWBDHV" + "ACGT" * 30], "NT"),
    (["MKVLAAGSTEEFIQP*", "MSTUOBZJX"], "AA"),
    (["RYKMSWDHVRYKM"], "AA"),  # nucleotide letters, but less than 90% plain bases: a peptide
    (["ACGT-ACGT"], "UK"),
    ([""], None),
])
def test_classify_sequences(sequences, expected):
    assert classify_sequences(sequences) == expected


def _with_codes(sequence, codes):
    """Every 97th base replaced, cycling through codes."""
    chars = list(sequence)
    for n, i in enumerate(range(50, len(chars), 97)):
        chars[i] = codes[n % len(codes)]
    return "".join(chars)


def test_genome_ambiguous_bases(tmp_path):
    """IUPAC codes in a genome are accepted and called like N (pyrodigal reads every non-ACGT letter as unknown)."""
    source = _GENOME_DIR / "input" / "ACIN21-1_SAMN05422111_MAG_00000002.fa.gz"
    records = read_fasta_gz(source)
    for label, codes in (("iupac", "RYKMSWBDHV"), ("n", "N")):
        (tmp_path / label).mkdir()
        write_gz(tmp_path / label / "genome.fa.gz",
                 "".join(f"{header}\n{_with_codes(sequence, codes)}\n" for header, sequence in records))
    check_sequence_type([tmp_path / "iupac" / "genome.fa.gz"], "NT")
    for label in ("iupac", "n"):
        extraction_genomes([tmp_path / label / "genome.fa.gz"], tmp_path / f"out_{label}", mode="genome", threads=2,
                           very_best=True)
    outputs = sorted(p.name for p in (tmp_path / "out_n").iterdir())
    assert outputs == sorted(p.name for p in (tmp_path / "out_iupac").iterdir())
    for name in outputs:
        assert (tmp_path / "out_iupac" / name).read_text() == (tmp_path / "out_n" / name).read_text(), name
    assert (tmp_path / "out_n" / "genome.fa.gz.fetchMGs.faa").read_text(), "no marker genes found"


def test_gene_ambiguous_bases_in_nucleotide_file(tmp_path):
    """-d with IUPAC codes in the genes passes the check; the marker genes keep the input's letters."""
    name = "ACIN21-1_SAMN05422111_MAG_00000002.genes"
    source = _GENOME_GENES_DIR / "input"
    records = read_fasta_gz(source / f"{name}.fna.gz")
    fna = tmp_path / "genes.fna.gz"
    write_gz(fna, "".join(f"{header}\n{_with_codes(sequence, 'RYKMSWBDHV')}\n" for header, sequence in records))
    check_sequence_type([fna], "NT")
    extraction_genes([source / f"{name}.faa.gz"], [fna], tmp_path / "out", threads=1, very_best=True)
    expected_dir = _GENOME_GENES_DIR / f"output_{name}"
    for kind in ("faa", "scores"):
        assert (tmp_path / "out" / f"{name}.faa.gz.fetchMGs.{kind}").read_text() == \
            read_gz(expected_dir / f"{name}.faa.gz.fetchMGs.{kind}.gz")
    by_name = {header[1:].split()[0]: sequence for header, sequence in records}
    markers = parse_fasta((tmp_path / "out" / f"{name}.faa.gz.fetchMGs.fna").read_text())
    assert markers
    for header, sequence in markers:
        gene = header[1:].rsplit(".", 1)[0]
        assert sequence == _with_codes(by_name[gene], "RYKMSWBDHV")


@pytest.mark.parametrize("mode", ["gene", "genome", "metagenome"])
def test_empty_input(mode, tmp_path, monkeypatch):
    """An input without sequences gives empty outputs and exit code 0."""
    empty = tmp_path / "empty.fa.gz"
    write_gz(empty, "")
    args = ["fetchMGs", "extraction", str(empty), mode, str(tmp_path / "out"), "-t", "2", "-b", "10"]
    if mode == "gene":
        empty_fna = tmp_path / "empty.fna.gz"
        write_gz(empty_fna, "")
        args += ["-d", str(empty_fna)]
    monkeypatch.setattr("sys.argv", args)
    with pytest.raises(SystemExit) as exit_info:
        main()
    assert exit_info.value.code == 0
    out = tmp_path / "out"
    assert (out / "empty.fa.gz.fetchMGs.faa").read_text() == ""
    assert (out / "empty.fa.gz.fetchMGs.scores").read_text() == "#protein_sequence_id\tHMM bit score\tCOG\n"
    assert (out / "empty.fa.gz.fetchMGs.fna").read_text() == ""


@pytest.mark.parametrize("copy_first", [False, True])
def test_very_best_tie_goes_to_first_in_file(copy_first, tmp_path):
    """-v: of two genes with the same score for a COG, the one earlier in the input file is reported."""
    name = "ACIN21-1_SAMN05422111_MAG_00000002.genes"
    source = _GENOME_GENES_DIR / "input"
    expected_scores = read_gz(_GENOME_GENES_DIR / f"output_{name}" / f"{name}.faa.gz.fetchMGs.scores.gz")
    marker = expected_scores.splitlines()[1].split("\t")[0]
    records = read_fasta_gz(source / f"{name}.faa.gz")
    original = next(r for r in records if r[0][1:].split()[0] == marker)
    copy = [f">zz_copy_of_{marker}", original[1]]
    records = [copy] + records if copy_first else records + [copy]
    faa = tmp_path / "genes.faa.gz"
    write_gz(faa, "".join(f"{header}\n{sequence}\n" for header, sequence in records))
    extraction_genes([faa], [None], tmp_path / "out", threads=2, very_best=True, chunk_size=50)
    scores = (tmp_path / "out" / "genes.faa.gz.fetchMGs.scores").read_text()
    winner, loser = (f"zz_copy_of_{marker}", marker) if copy_first else (marker, f"zz_copy_of_{marker}")
    assert f"\n{winner}\t" in scores and f"\n{loser}\t" not in scores
