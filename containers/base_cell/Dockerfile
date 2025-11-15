FROM python:3.12-slim

WORKDIR /app

# Copy pyproject.toml and optional lockfile
COPY pyproject.toml uv.lock* /app/

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential git procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Verify procps is installed
RUN which ps || (echo "ps not found!" && exit 1)

# Install uv and Python dependencies
RUN pip install --no-cache-dir uv
RUN uv lock
RUN uv sync

# Add the virtual environment to PATH
ENV PATH="/app/.venv/bin:$PATH"

# Copy the bin directory with base_cell.py
COPY bin/ /app/bin/