FROM nfcore/base:latest

LABEL author="Nextflow Pipeline Maintainers" \
      description="Docker image containing dependencies for Nextflow pipeline testing"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    fastqc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    fastq-generator \
    nf-test==0.9.2

# Install fastp
RUN wget -q http://opengene.org/fastp/fastp_0.23.2 -O /usr/local/bin/fastp && \
    chmod a+x /usr/local/bin/fastp

# Set workdir
WORKDIR /pipeline

# Add test entry point
COPY .github/workflows/test-entrypoint.sh /usr/local/bin/test-entrypoint.sh
RUN chmod +x /usr/local/bin/test-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/test-entrypoint.sh"]
