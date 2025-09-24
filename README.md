# CRISPR Guide Design Toolkit

A high-performance toolkit for designing and analyzing CRISPR guide RNAs with advanced features including degenerate PAM support, genome-scale processing, and comprehensive off-target analysis.

## 🔬 Key Features

- **Degenerate PAM Support**: Use IUPAC nucleotide codes (V, R, Y, W, S, K, M, B, D, H, N) for precise PAM matching
- **Circular Genome Support**: Seamlessly handles circular chromosomes, plasmids, and organellar genomes
- **Genome-Scale Processing**: Handle millions of guides efficiently with memory-optimized streaming
- **Advanced Off-Target Detection**: Bowtie-based alignment with configurable mismatch tolerance
- **Multiple CRISPR Systems**: Support for Cas9, Cas12a, and custom PAM sequences
- **Genomic Coordinate Sorting**: Results ordered by chromosome and position for easy analysis
- **Gene Annotation**: Automatic mapping of guides to genes with locus tags and feature types
- **High-Performance Computing**: Polars-based data processing for large datasets

## 🚀 Quick Start

### Installation

1. **Install Miniforge** (if not already installed):
   ```bash
   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
   bash Miniforge3-$(uname)-$(uname -m).sh
   ```

2. **Clone and setup**:
   ```bash
   git clone https://github.com/ryandward/crispr_experiment.git
   cd crispr_experiment
   mamba env create -f environment.yml -n crispr_experiment
   mamba activate crispr_experiment
   ```

### Basic Usage

**Design guides for Cas12a (TTTV PAM)**:
```bash
python assets/find_guides.py \
    --genome-file your_genome.gb \
    --pam TTTV \
    --barcode-length 24 \
    --mismatches 0 \
    --pam-direction upstream > cas12a_guides.tsv
```

**Design guides for Cas9 (NGG PAM)**:
```bash
python assets/find_guides.py \
    --genome-file your_genome.gb \
    --pam NGG \
    --barcode-length 20 \
    --mismatches 2 \
    --pam-direction downstream > cas9_guides.tsv
```

## 🧬 Degenerate PAM Examples

This toolkit supports IUPAC degenerate nucleotide codes for precise PAM matching:

| Code | Matches | Example | Use Case |
|------|---------|---------|----------|
| `V` | A, C, G | `TTTV` → TTTA, TTTC, TTTG | Cas12a (excludes TTTT) |
| `R` | A, G | `TTTR` → TTTA, TTTG | Purine-specific |
| `Y` | C, T | `TTTY` → TTTC, TTTT | Pyrimidine-specific |
| `W` | A, T | `TTTW` → TTTA, TTTT | Weak bonds |
| `S` | C, G | `TTTS` → TTTC, TTTG | Strong bonds |
| `N` | A, C, G, T | `TTTN` → All combinations | Any nucleotide |

**Example**: `TTTV` matches TTTA, TTTC, TTTG but excludes TTTT - perfect for Cas12a specificity.

## 📊 Output Format

Results are provided in tab-separated format with genomic coordinate sorting:

```
spacer                   locus_tag      gene           chr          pam   mismatches  tar_start   tar_end
ATCGATCGATCGATCGATCG     LOC101252303   gene_name      NC_015444.3  TTTC  0          12345       12365
GCTAGCTAGCTAGCTAGCTA     LOC101252304   another_gene   NC_015444.3  TTTA  0          12400       12420
```

## ⚡ Performance Features

- **Memory Efficient**: Process 30+ million guides without crashes
- **Streaming Processing**: Results written incrementally to handle large datasets
- **Parallel Alignment**: Multi-threaded Bowtie alignment for speed
- **Smart Deduplication**: Genomic coordinate-based duplicate removal
- **Optimized I/O**: Direct file output bypasses memory bottlenecks

## 🧮 Analysis Modes

### Unique Mapping Only (`-k 1`)
```bash
# Only guides with unique genome targets
--mismatches 0  # Perfect matches only
```

### Off-Target Detection (`-k 2+`)
```bash
# Include potential off-targets for analysis
--mismatches 2  # Allow up to 2 mismatches
```

### Gene-Specific Targeting
```bash
# Target specific genomic regions
--regulatory 500  # Include 500bp regulatory regions
```

## 🛠 Advanced Usage

### Large Genome Processing
For genomes with millions of potential guides:

```bash
# Use minimal mismatches and unique mapping
python assets/find_guides.py \
    --genome-file large_genome.gb \
    --pam TTTV \
    --barcode-length 24 \
    --mismatches 0 \
    --pam-direction upstream > results.tsv
```

### Custom PAM Sequences
```bash
# Complex degenerate PAMs
python assets/find_guides.py \
    --pam NNGRRT \  # N=any, G=G, R=A/G, T=T
    --pam-direction downstream
```

### Circular Genomes
```bash
# Plasmids, organellar genomes, bacterial chromosomes
python assets/find_guides.py \
    --genome-file plasmid.gb \     # Automatically detects circular topology
    --pam TTTV \
    --barcode-length 24 > plasmid_guides.tsv

# Handles coordinate wrapping at origin automatically
```

## 📦 Available Tools

### Command Line Interface
- `find_guides.py` - Main guide design tool
- `targets.py` - Alignment and annotation engine

### Graphical Interface
- `crispr_experiment.py` - Unified interface for all tools
- `find_guides_gui.py` - Interactive guide design
- `design_mismatches_gui.py` - Mismatch analysis
- `assembly_finder_gui.py` - Genome database search

## 🔧 Dependencies

Core computational dependencies:
- **Bowtie**: Fast genome alignment
- **Polars**: High-performance data processing  
- **BioPython**: Sequence analysis
- **PyArrow**: Efficient data formats
- **Rich**: Enhanced terminal output

## 📈 Scalability

Successfully tested with:
- **Tomato genome**: 900+ Mb, 31K+ genes, linear chromosomes
- **Bacterial plasmids**: Circular topology, seamless coordinate wrapping
- **Guide discovery**: 33+ million potential sites
- **Final output**: Millions of unique guides
- **Processing time**: ~30 minutes on standard hardware

## 🆘 Troubleshooting

### Memory Issues
```bash
# Reduce batch size for very large genomes
export BATCH_SIZE=25000

# Use faster storage for temp files
export TMPDIR=/path/to/fast/storage
```

### Performance Optimization
```bash
# Increase thread count
--threads 32

# Use unique mapping for faster processing
--mismatches 0
```

## 📚 Citation

If you use this toolkit in your research, please cite:
```
[Citation information to be added]
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md).

## 📄 License

[License information to be added]

## 🐛 Support

- [Issues](https://github.com/ryandward/crispr_experiment/issues)
- [Discussions](https://github.com/ryandward/crispr_experiment/discussions)
