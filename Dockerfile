# Dockerfile for running dlt pipelines

# Dockerfile is based heavily on the example uv dockerfile:
# https://github.com/astral-sh/uv-docker-example

# Use a Python image with uv pre-installed
FROM ghcr.io/astral-sh/uv:python3.14-trixie

ARG QSV_VERSION="25.0.0"
ARG XML_FILE_SPLITTER_VERSION="v0.1.4"

# Set environment variable to noninteractive to prevent prompts during apt operations
ENV DEBIAN_FRONTEND=noninteractive

# add tini and git
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    libwayland-client0 \
    libxml2-utils \
    parallel \
    rclone \
    tini \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# download and install the xml_file_splitter and qsv binaries, and copy them to /usr/local/bin
RUN ARCH=$(uname -m) && \
    wget https://github.com/ialarmedalien/xml_file_splitter/releases/download/${XML_FILE_SPLITTER_VERSION}/xml_file_splitter-${ARCH}-unknown-linux-gnu.tar.gz && \
    tar -xvf xml_file_splitter-${ARCH}-unknown-linux-gnu.tar.gz && \
    mv xml_file_splitter-${ARCH}-unknown-linux-gnu/xml_file_splitter /usr/local/bin/ && \
    # qsv release -- only need the `qsv` binary from it
    wget https://github.com/dathere/qsv/releases/download/${QSV_VERSION}/qsv-${QSV_VERSION}-${ARCH}-unknown-linux-gnu.zip && \
    unzip qsv-${QSV_VERSION}-${ARCH}-unknown-linux-gnu.zip -d /tmp/qsv && \
    mv /tmp/qsv/qsv /usr/local/bin/ && \
    rm -rf /tmp/*

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
    && useradd --system --gid 999 --uid 999 --create-home nonroot

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Omit development dependencies
# ENV UV_NO_DEV=1

# don't try to synchronise each time uv is executed
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
COPY --chown=nonroot:nonroot .dlt /app/.dlt
COPY --chown=nonroot:nonroot docs /app/docs
COPY --chown=nonroot:nonroot notebooks /app/notebooks
COPY --chown=nonroot:nonroot scripts /app/scripts
COPY --chown=nonroot:nonroot src /app/src
COPY --chown=nonroot:nonroot tests /app/tests
COPY --chown=nonroot:nonroot README.md /app/README.md
COPY --chown=nonroot:nonroot pyproject.toml /app/pyproject.toml
COPY --chown=nonroot:nonroot uv.lock /app/uv.lock

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-editable

RUN chown -R nonroot:nonroot /app
RUN chmod +x ./scripts/*.sh

# create a home dir for the user for dlt config files
# files in $HOME override those in the execution directory
# Note: set HOME after `uv sync` so root build steps use /root/.cache/uv
# rather than creating root-owned cache files in /home/nonroot/.cache/uv
RUN mkdir -p /home/nonroot && chown -R nonroot:nonroot /home/nonroot
ENV HOME="/home/nonroot"

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

USER nonroot
ENTRYPOINT ["./scripts/entrypoint.sh"]
