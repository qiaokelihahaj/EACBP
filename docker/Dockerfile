ARG PYTHON_BASE=python:3.12.14-slim-bookworm@sha256:392307d22300de8b5986851a12d9176dfc0fc073e65bf6523ebd7dcbeb23564e

FROM ${PYTHON_BASE} AS wheel-builder
WORKDIR /build
COPY docker/requirements-linux-py312-build.lock ./requirements-build.lock
RUN python -m pip install --no-cache-dir --timeout=120 --index-url https://pypi.org/simple \
        --require-hashes -r requirements-build.lock
COPY pyproject.toml README.md LICENSE ./
COPY eacbp/ ./eacbp/
COPY slurm/run_eacbp_container.sbatch ./slurm/run_eacbp_container.sbatch
RUN python -m pip wheel --no-deps --no-build-isolation --wheel-dir=/wheels .

FROM ${PYTHON_BASE} AS runtime

ENV VIRTUAL_ENV=/opt/venv \
    PATH="/opt/venv/bin:${PATH}" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONNOUSERSITE=1 \
    CUDA_VISIBLE_DEVICES="" \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    NUMBA_NUM_THREADS=1 \
    HOME=/scratch/home \
    XDG_CACHE_HOME=/scratch/cache \
    NUMBA_CACHE_DIR=/scratch/cache/numba \
    MPLCONFIGDIR=/scratch/cache/matplotlib \
    TMPDIR=/scratch/tmp

RUN python -m venv /opt/venv \
    && mkdir -p /opt/eacbp-build /outputs /scratch \
    && chown -R 10001:10001 /outputs /scratch
COPY docker/requirements-linux-py312-scientific.lock /opt/eacbp-build/requirements-scientific.lock
RUN /opt/venv/bin/python -m pip install --no-cache-dir --timeout=120 \
        --index-url https://pypi.org/simple --require-hashes \
        -r /opt/eacbp-build/requirements-scientific.lock
# Keep the locked scientific environment cached when only EACBP source changes.
COPY --from=wheel-builder /wheels/ /wheels/
RUN /opt/venv/bin/python -m pip install --no-cache-dir --no-deps /wheels/eacbp-0.1.0-py3-none-any.whl \
    && /opt/venv/bin/python -m pip check \
    && /opt/venv/bin/python -c "import anndata, cellrank, celltypist, decoupler, eacbp, liana, pydeseq2, scanpy, skimage" \
    && /opt/venv/bin/python -m pip freeze --all > /opt/eacbp-build/pip-freeze.txt \
    && rm -rf /wheels

# Scientific import checks create caches as root during the build. The default
# runtime user must also be able to use these paths without a scratch mount.
RUN chown -R 10001:10001 /scratch

COPY --chmod=0555 docker/entrypoint.sh /usr/local/bin/eacbp-entrypoint
WORKDIR /outputs
USER 10001:10001
ENTRYPOINT ["/usr/local/bin/eacbp-entrypoint"]
CMD ["/opt/venv/bin/eacbp", "--help"]

# Revision metadata must not invalidate the cached scientific environment.
ARG EACBP_SOURCE_REVISION=unknown
LABEL org.opencontainers.image.title="EACBP CPU analysis" \
      org.opencontainers.image.revision="${EACBP_SOURCE_REVISION}"
