# Dockerfile for running dlt pipelines

# Dockerfile is based heavily on the example uv dockerfile:
# https://github.com/astral-sh/uv-docker-example

# Use a Python image with uv pre-installed
FROM ghcr.io/astral-sh/uv:python3.13-trixie-slim

# Set environment variable to noninteractive to prevent prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive

# add tini
RUN apt-get update -y && apt-get install -y --no-install-recommends tini git

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Omit development dependencies
# ENV UV_NO_DEV=1
ENV UV_NO_SYNC=1

# Ensure installed tools can be executed out of the box
ENV UV_TOOL_BIN_DIR=/usr/local/bin

# Install the project into `/app`
WORKDIR /app

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-editable

# Then, add the rest of the project source code and install it
# Installing separately from its dependencies allows optimal layer caching
COPY --chown=nonroot:nonroot . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

COPY --chmod=+x ./scripts/entrypoint.sh /app/
# Use the non-root user to run our application
USER nonroot
ENTRYPOINT ["./entrypoint.sh"]
