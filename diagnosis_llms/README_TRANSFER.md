# Notebook Transfer Bundle

This bundle contains 3 notebooks moved from `diagnosis_llms` plus required support files.

## Included notebooks
- `notebooks/display_reasoning.ipynb`
- `notebooks/multi_class_cont.ipynb`
- `notebooks/test_notebook.ipynb`

## Included support code
- `src/diagnosis_llms/paths.py`
- `src/diagnosis_llms/__init__.py`

These notebooks currently import `from diagnosis_llms.paths import get_paths` in their first cell.
If your target repo does not include this module, either:
1. copy `src/diagnosis_llms/` from this bundle, or
2. edit notebook first cells to use direct `Path` logic without that import.

## Environment setup options
- `env/pyproject.notebook-transfer-snippet.toml`: merge this dependency group into target `pyproject.toml`.
- `env/requirements.notebook-transfer.txt`: pip/uv requirements alternative.
- `env/kernel_setup.sh`: installs deps + registers Jupyter kernel.

## Recommended setup in target repo
1. Merge the snippet from `env/pyproject.notebook-transfer-snippet.toml` into your target `pyproject.toml`.
2. Run:
   - `uv sync --group notebook-transfer`
   - `uv run --group notebook-transfer python -m ipykernel install --user --name notebook-transfer --display-name "Python (notebook-transfer)"`
3. In VS Code, select kernel `Python (notebook-transfer)`.

## Runtime caveats
- `test_notebook.ipynb` is exploratory and pulls external resources (Hugging Face, GitHub, web downloads).
- It also expects shell tools in some cells (`wget`, `grep`) and may create local files/directories.
- Cells using OpenAI require `OPENAI_API_KEY` in environment.
