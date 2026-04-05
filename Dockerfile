FROM rust:1.80-slim-bookworm AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src/ src/
RUN cargo build --release

# ---

FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    gzip \
    bc \
    python3 \
    python3-pip \
    hmmer \
    && rm -rf /var/lib/apt/lists/*

# Install Prodigal (baseline gene finder)
RUN curl -sL https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux -o /usr/local/bin/prodigal \
    && chmod +x /usr/local/bin/prodigal

# Install barrnap
RUN curl -sL https://github.com/tseemann/barrnap/archive/refs/tags/0.9.tar.gz | tar xz -C /opt \
    && ln -s /opt/barrnap-0.9/bin/barrnap /usr/local/bin/barrnap

# Copy binary, model, and scripts
COPY --from=builder /build/target/release/prokrustes /usr/local/bin/prokrustes
COPY models/ /opt/prokrustes/models/
COPY scripts/ /opt/prokrustes/scripts/

# Download all 10 reference genomes
RUN mkdir -p /data && \
    echo "Downloading reference genomes..." && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" | gunzip > /data/ecoli_k12.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz" | gunzip > /data/ecoli_k12.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz" | gunzip > /data/salmonella_lt2.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz" | gunzip > /data/salmonella_lt2.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" | gunzip > /data/bacillus_subtilis.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz" | gunzip > /data/bacillus_subtilis.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz" | gunzip > /data/pseudomonas_aeruginosa.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz" | gunzip > /data/pseudomonas_aeruginosa.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz" | gunzip > /data/saureus_nctc8325.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz" | gunzip > /data/saureus_nctc8325.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz" | gunzip > /data/mycobacterium_tuberculosis.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz" | gunzip > /data/mycobacterium_tuberculosis.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz" | gunzip > /data/mycoplasma_genitalium.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz" | gunzip > /data/mycoplasma_genitalium.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/725/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna.gz" | gunzip > /data/synechocystis_pcc6803.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/725/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.gff.gz" | gunzip > /data/synechocystis_pcc6803.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.fna.gz" | gunzip > /data/bacteroides_fragilis.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.gff.gz" | gunzip > /data/bacteroides_fragilis.gff && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/685/GCF_000008685.2_ASM868v2/GCF_000008685.2_ASM868v2_genomic.fna.gz" | gunzip > /data/borrelia_burgdorferi.fasta && \
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/685/GCF_000008685.2_ASM868v2/GCF_000008685.2_ASM868v2_genomic.gff.gz" | gunzip > /data/borrelia_burgdorferi.gff && \
    echo "All genomes downloaded."

WORKDIR /workspace

ENTRYPOINT ["prokrustes"]
CMD ["--help"]
