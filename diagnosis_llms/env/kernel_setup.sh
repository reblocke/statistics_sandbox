#!/usr/bin/env bash
set -euo pipefail

# Option A: using dependency-group in pyproject.toml
uv sync --group notebook-transfer

# Option B: using requirements file
# uv pip install -r env/requirements.notebook-transfer.txt

uv run --group notebook-transfer python -m ipykernel install --user \
  --name notebook-transfer \
  --display-name "Python (notebook-transfer)"

echo "Kernel installed: Python (notebook-transfer)"
