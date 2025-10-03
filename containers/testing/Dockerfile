FROM python:3.12-slim

WORKDIR /app

# Copy pyproject.toml and optional lockfile
COPY pyproject.toml poetry.lock* /app/

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential git && \
    rm -rf /var/lib/apt/lists/*

# Install uv and Python dependencies
RUN pip install --no-cache-dir uv
RUN uv sync