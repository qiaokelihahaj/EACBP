ARG PYTHON_BASE=python:3.12-slim-bookworm

FROM ${PYTHON_BASE} AS wheel-builder
WORKDIR /build
COPY pyproject.toml README.md LICENSE ./
COPY eacbp/ ./eacbp/
COPY slurm/run_eacbp_container.sbatch ./slurm/run_eacbp_container.sbatch
RUN python -m pip install --no-cache-dir "setuptools>=61" wheel \
    && python -m pip wheel --no-deps --no-build-isolation --wheel-dir=/wheels .

FROM ${PYTHON_BASE} AS runtime
ARG EACBP_SOURCE_REVISION=unknown
LABEL org.opencontainers.image.title="EACBP CPU analysis" \
      org.opencontainers.image.revision="${EACBP_SOURCE_REVISION}"

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
COPY --from=wheel-builder /wheels/ /wheels/
RUN /opt/venv/bin/python -m pip install --no-cache-dir \
        "eacbp[standard,fate,advanced-statistics] @ file:///wheels/eacbp-0.1.0-py3-none-any.whl" \
    && /opt/venv/bin/python -m pip check \
    && /opt/venv/bin/python -c "import anndata, cellrank, decoupler, eacbp, pydeseq2, scanpy" \
    && /opt/venv/bin/python -m pip freeze --all > /opt/eacbp-build/pip-freeze.txt \
    && rm -rf /wheels

COPY --chmod=0555 docker/entrypoint.sh /usr/local/bin/eacbp-entrypoint
WORKDIR /outputs
USER 10001:10001
ENTRYPOINT ["/usr/local/bin/eacbp-entrypoint"]
CMD ["/opt/venv/bin/eacbp", "--help"]
