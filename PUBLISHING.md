# Publishing fetchMGs to PyPI

## Prerequisites

```bash
pip install build twine
```

## Steps

### 1. Bump the version

Edit `pyproject.toml` — the `version` field is the single source of truth:

```toml
[project]
version = "2.1.X"
```

Also update `__date__` in [fetchmgs/fetchmgs.py](fetchmgs/fetchmgs.py) if desired.

### 2. Clean previous build artifacts

```bash
rm -rf dist/
```

### 3. Build the source distribution and wheel

```bash
python -m build
```

This produces `dist/fetchMGs-<version>.tar.gz` (sdist) and `dist/fetchMGs-<version>-py3-none-any.whl` (wheel).

### 4. Verify the distributions

```bash
twine check dist/*
```

Fix any warnings before uploading.

### 5. (Optional) Test on TestPyPI first

```bash
twine upload --repository testpypi dist/*
```

Install from TestPyPI to verify:

```bash
pip install --index-url https://test.pypi.org/simple/ fetchMGs
```

### 6. Upload to PyPI

```bash
twine upload dist/*
```

Enter your PyPI username (`__token__`) and API token when prompted, or configure `~/.pypirc`:

```ini
[pypi]
username = __token__
password = pypi-<your-token>
```
