# Alliance of Genome Resources on AWS

Welcome to the Alliance of Genome Resources (Alliance) open data repository on AWS. This documentation helps you access and use comprehensive genomic, genetic, and molecular data from multiple model organisms.

---

## What is the Alliance of Genome Resources?

The Alliance of Genome Resources is a consortium integrating data from leading model organism databases:

- 🪰 **Drosophila melanogaster** and other Drosophila species
- 🪱 **Caenorhabditis elegans**
- 🐟 **Danio rerio** (zebrafish)
- 🐭 **Mus musculus** (mouse)
- 🐀 **Rattus norvegicus** (rat)
- 🍺 **Saccharomyces cerevisiae** (yeast)
- 🐸 **Xenopus laevis** and **Xenopus tropicalis** (frogs)
- 👤 **Homo sapiens** (human reference data)

**Mission:** Provide unified, high-quality genomic data to accelerate biological research and human disease understanding.

---

## Quick Start

### 1. Browse Available Data

Visit the Alliance downloads page or explore the S3 buckets:

**FlyBase Data (Public HTTPS):**
```bash
# Browse via wget
wget -qO- https://s3ftp.flybase.org/releases/current/precomputed_files/ | grep '.tsv.gz'
```

**Alliance Data (S3):**
```bash
# List Alliance releases
aws s3 ls s3://mod-datadumps/ --no-sign-request

# List disease data in latest release
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/ --no-sign-request
```

### 2. Download Your First File

```bash
# Download gene annotation IDs
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz

# Decompress
gunzip fbgn_annotation_ID_current.tsv.gz

# View
head fbgn_annotation_ID_current.tsv
```

### 3. Explore Tutorials

Check out [TUTORIAL.md](./TUTORIAL.md) for step-by-step guides on:
- Downloading and accessing data
- Gene annotation lookups
- Disease gene discovery
- Expression analysis
- Building interaction networks
- Cross-species orthology mapping

---

## Documentation

### 📚 [Data Organization Documentation](./DATA_DOCUMENTATION.md)

Complete reference guide covering:
- **Dataset Overview** - Scale, update frequency, data categories
- **S3 Bucket Structure** - Directory organization and file naming
- **Data Categories** - Detailed descriptions of all data types:
  - Gene annotations
  - Expression data (bulk RNA-Seq, single-cell RNA-Seq)
  - Disease associations
  - Molecular and genetic interactions
  - Orthology relationships
  - Variants and alleles
  - Genome sequences and annotations
- **Access Methods** - S3, HTTPS, Python boto3, PostgreSQL, API
- **File Formats** - TSV, JSON, FASTA, GFF3, GTF, VCF, MITAB
- **Common Use Cases** - Real-world examples

### 🎓 [Tutorial Guide](./TUTORIAL.md)

Hands-on tutorials with working code examples:

1. **Getting Started** - Download your first dataset (10 min)
2. **Gene Annotation Lookups** - ID conversion and batch processing (15 min)
3. **Disease Gene Discovery** - Find disease-associated genes (20 min)
4. **RNA-Seq Expression Analysis** - Analyze developmental expression (25 min)
5. **Protein Interaction Networks** - Build and visualize PPI networks (20 min)
6. **Cross-Species Orthology** - Map genes between species (15 min)

---

## Data Categories

### Core Annotation Data
- **Gene Annotations** - IDs, symbols, descriptions, mappings
- **Gene Ontology** - Molecular function, biological process, cellular component
- **Gene Groups** - Pathway and functional groupings

### Functional Genomics
- **Bulk RNA-Seq** - RPKM/FPKM expression matrices across tissues and stages
- **Single-Cell RNA-Seq** - Cell cluster expression from multiple datasets
- **Curated Expression** - Spatiotemporal expression patterns with ontology terms

### Disease and Phenotypes
- **Disease Associations** - Links to human diseases (Disease Ontology)
- **Phenotypes** - Mutant and variant phenotype annotations
- **Human Disease Models** - Model organism connections to human disease

### Interactions
- **Physical Interactions** - Protein-protein, protein-RNA, RNA-RNA (PSI-MI TAB format)
- **Genetic Interactions** - Suppression, enhancement, synthetic lethality

### Comparative Genomics
- **Orthology** - Cross-species gene relationships with DIOPT scores
- **Paralogy** - Within-species gene duplications

### Variants and Alleles
- **VCF Files** - Genomic variants in standard format
- **Allele Annotations** - Detailed allele and variant descriptions
- **Genotype-Phenotype** - Links between genetic changes and phenotypes

### Genome Sequences
- **FASTA** - Chromosomes, genes, transcripts, proteins
- **GFF3/GTF** - Genome annotations for analysis pipelines
- **Transposable Elements** - TE sequences and insertion sites

---

## Access Methods

### 🌐 Web Browser
- **Alliance Portal:** https://www.alliancegenome.org
- **Downloads Page:** https://www.alliancegenome.org/downloads
- **FTP Browser:** https://s3ftp.flybase.org/releases/current/

### ☁️ AWS S3 (Recommended for Large Downloads)
```bash
# Anonymous access - no AWS account needed
# FlyBase data
aws s3 ls s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/ --no-sign-request
aws s3 cp s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz . --no-sign-request

# Alliance data
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/ --no-sign-request
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --no-sign-request
```

### 📥 Direct Download (wget/curl)
```bash
# wget
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz

# curl
curl -O https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz
```

### 🐍 Python (boto3)
```python
import boto3
from botocore import UNSIGNED
from botocore.client import Config

# Create S3 client with anonymous access
s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

# Download FlyBase data
s3.download_file('s3ftp.flybase.org',
                'releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz',
                'fbgn_annotation_ID_current.tsv.gz')

# Download Alliance data
s3.download_file('mod-datadumps',
                '8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz',
                'DISEASE-ALLIANCE_COMBINED_2.tsv.gz')
```

### 🗄️ PostgreSQL Database
```bash
# Public read-only access
psql -h chado.flybase.org -U flybase flybase
```

### 🔌 REST API
```bash
# Get gene information
curl https://www.alliancegenome.org/api/gene/FBgn0000001

# Search genes
curl https://www.alliancegenome.org/api/search?category=gene&q=white
```

---

## File Naming Conventions

Alliance files follow predictable patterns:

### Current Release (Always Latest)
```
<data_type>_<details>_current.tsv.gz
dmel-all-<feature>-current.fasta.gz
```

Examples:
- `fbgn_annotation_ID_current.tsv.gz`
- `gene_rpkm_matrix_current.tsv.gz`

### Versioned Releases
```
<data_type>_<details>_fb_YYYY_MM.tsv.gz
dmel-all-<feature>-rX.YY.fasta.gz
```

Examples:
- `fbgn_annotation_ID_fb_2023_06.tsv.gz` (FB2023_06 release)
- `dmel-all-chromosome-r6.55.fasta.gz` (genome release 6.55)

---

## Data Updates

**Release Schedule:**
- **Major Releases:** Quarterly (every ~3 months)
- **Hot Fixes:** As needed for critical corrections
- **Continuous Updates:** Daily for time-sensitive annotations

**Versioning:**
- Alliance releases: `FB[YEAR]_[MONTH]` (e.g., FB2023_06)
- Genome releases: `r[MAJOR].[MINOR]` (e.g., r6.55)

**Current vs. Archived:**
- `/releases/current/` - Always points to latest release
- `/releases/FB2023_06/` - Specific archived release

---

## Common Use Cases

### 🧬 Gene Research
- Convert gene symbols to database IDs
- Retrieve gene descriptions and annotations
- Find genes in specific pathways or GO terms

### 🏥 Disease Studies
- Identify genes associated with human diseases
- Find model organism disease models
- Map disease genes to orthologs

### 📊 Expression Analysis
- Compare gene expression across developmental stages
- Analyze tissue-specific expression
- Explore single-cell expression patterns

### 🔗 Interaction Networks
- Build protein-protein interaction networks
- Analyze genetic interactions
- Find interaction partners for proteins of interest

### 🧫 Comparative Genomics
- Map orthologs between species
- Find conserved genes and pathways
- Compare genomic features across model organisms

### 🧪 Variant Analysis
- Access genomic variants in VCF format
- Link variants to phenotypes
- Study allele effects

---

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| **TSV** | Tab-separated values | General data tables |
| **JSON** | JavaScript Object Notation | Structured data, API responses |
| **FASTA** | Sequence data | Genomic/protein sequences |
| **GFF3** | Genome annotations | Genome browsers, analysis |
| **GTF** | Gene Transfer Format | RNA-Seq pipelines |
| **VCF** | Variant Call Format | Variant analysis |
| **MITAB** | PSI-MI TAB | Protein interactions |
| **XML** | Chado XML | Complete database dumps |

All compressed files use **gzip** (`.gz` extension).

---

## System Requirements

### Minimal Setup (Browser Only)
- Web browser
- No special software needed

### Command Line Access
- Linux, Mac, or Windows (WSL)
- `wget` or `curl`
- `gunzip` (usually pre-installed)

### Data Analysis
**Python:**
- Python 3.7+
- pandas, boto3, matplotlib, biopython

**R:**
- R 4.0+
- tidyverse, data.table

**Tools:**
- AWS CLI (optional, for S3 access)
- IGV, JBrowse (genome visualization)
- Cytoscape (network analysis)

---

## Citation and License

### Citing Alliance Data

**Primary Citation:**
> Alliance of Genome Resources Consortium. Alliance of Genome Resources Portal: unified model organism research platform. *Nucleic Acids Research* (2023). https://doi.org/10.1093/nar/gkac1003

### Data License

Most Alliance data is available under **CC0 1.0 Universal (Public Domain Dedication)**. Some datasets may use **CC-BY 4.0** (attribution required).

**License Details:** https://www.alliancegenome.org/terms-of-use

### Attribution Requirements

When publishing research using Alliance data:
1. Cite the Alliance consortium paper (above)
2. Include release version numbers for reproducibility
3. Acknowledge specific data sources when applicable
4. Link to https://www.alliancegenome.org in web applications

---

## Support and Help

### Documentation
- **Data Documentation:** [DATA_DOCUMENTATION.md](./DATA_DOCUMENTATION.md)
- **Tutorials:** [TUTORIAL.md](./TUTORIAL.md)
- **Alliance Homepage:** https://www.alliancegenome.org
- **API Docs:** https://www.alliancegenome.org/api/swagger-ui

### Get Help
- **Email:** help@alliancegenome.org
- **GitHub Issues:** https://github.com/alliance-genome/agr_java_software/issues

### Community
- **Twitter:** @alliancegenome
- **Mailing Lists:** Check Alliance website for announcements

---

## AWS Open Data Sponsorship Program

This dataset is part of the [AWS Open Data Sponsorship Program](https://registry.opendata.aws), which provides free hosting for publicly available high-value datasets.

**Benefits:**
- ✓ No AWS account required for downloads
- ✓ No data egress fees
- ✓ High-speed S3 access
- ✓ Global availability
- ✓ Automatic backups and archiving

**Registry Entry:** https://registry.opendata.aws/alliance-genome-resources/

---

## Quick Reference

### Essential Files

| Data Type | File | Location |
|-----------|------|----------|
| Gene IDs | `fbgn_annotation_ID_*.tsv.gz` | `precomputed_files/genes/` |
| Expression Matrix | `gene_rpkm_matrix_*.tsv.gz` | `precomputed_files/genes/` |
| Disease Data | `disease_model_annotations_*.tsv.gz` | `precomputed_files/disease/` |
| Interactions | `physical_interactions_mitab_*.tsv.gz` | `precomputed_files/interactions/` |
| Orthologs | `dmel_human_orthologs_disease_*.tsv.gz` | `precomputed_files/orthologs/` |
| Genome Sequence | `dmel-all-chromosome-*.fasta.gz` | `genomes/.../fasta/` |
| Genome Annotation | `dmel-all-*.gff.gz` | `genomes/.../gff/` |

### Example Commands

```bash
# Download gene annotations via wget
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz

# Download FlyBase data using AWS CLI
aws s3 cp s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz . --no-sign-request

# Download Alliance data using AWS CLI
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --no-sign-request

# List files in FlyBase
aws s3 ls s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/ --no-sign-request

# List files in Alliance
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/ --no-sign-request

# Decompress
gunzip fbgn_annotation_ID_current.tsv.gz

# View first 10 data rows (skip comment lines)
grep -v '^#' fbgn_annotation_ID_current.tsv | head -10
```

---

## Changelog

### Version 1.0 (2025-10-17)
- Initial documentation release
- Comprehensive data organization guide
- Six hands-on tutorials
- AWS S3 access instructions
- Registry of Open Data submission

---

## Contributing

We welcome feedback and contributions:

1. **Report Issues:** Use GitHub issues for bugs or documentation improvements
2. **Suggest Tutorials:** Email help@alliancegenome.org with tutorial ideas
3. **Share Use Cases:** Tell us how you're using Alliance data
4. **Contribute Code:** Submit pull requests for example scripts

---

## Additional Resources

### Alliance Components
- **WormBase:** https://www.wormbase.org
- **ZFIN:** https://zfin.org
- **MGI:** http://www.informatics.jax.org
- **RGD:** https://rgd.mcw.edu
- **SGD:** https://www.yeastgenome.org
- **Xenbase:** http://www.xenbase.org

### Related Projects
- **Gene Ontology:** http://geneontology.org
- **Disease Ontology:** https://disease-ontology.org
- **UniProt:** https://www.uniprot.org
- **NCBI:** https://www.ncbi.nlm.nih.gov

### Tools and Browsers
- **JBrowse:** https://jbrowse.org
- **IGV:** https://software.broadinstitute.org/software/igv/
- **Cytoscape:** https://cytoscape.org

---

## Contact

**Alliance of Genome Resources**
- Website: https://www.alliancegenome.org
- Email: help@alliancegenome.org
- Twitter: @alliancegenome
- GitHub: https://github.com/alliance-genome

---

**Documentation Version:** 1.0
**Last Updated:** 2025-10-17
**Maintained by:** Alliance of Genome Resources Consortium
