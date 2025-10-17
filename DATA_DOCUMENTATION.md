# Alliance of Genome Resources - Data Organization Documentation

## Dataset Overview

The Alliance of Genome Resources (Alliance) is a comprehensive knowledge base providing integrated genomic, genetic, and molecular data across multiple model organisms. The Alliance aggregates and harmonizes data from model organism database components including:

- *Drosophila melanogaster* and other Drosophila species (FlyBase component)
- *Caenorhabditis elegans* (WormBase component)
- *Danio rerio* (ZFIN component)
- *Mus musculus* (MGI component)
- *Rattus norvegicus* (RGD component)
- *Saccharomyces cerevisiae* (SGD component)
- *Xenopus laevis* and *Xenopus tropicalis* (Xenbase component)
- *Homo sapiens* reference data

### Data Scale and Update Frequency

**Total Dataset Size:** Multiple terabytes of genomic annotations, sequences, and experimental data

**Update Schedule:**
- Major releases: Quarterly (approximately every 3 months)
- Alliance versioning: Numbered releases (e.g., 8.3.0, 8.2.0)
- FlyBase component versioning: FB[YEAR]_[MONTH] (e.g., FB2025_04)

**Data Access:**
1. Alliance-wide integrated data (cross-species)
2. Model organism-specific data (including FlyBase for Drosophila)

---

## Part 1: Alliance-Wide Integrated Data

### Alliance S3 Repository Structure

**Primary Bucket:** `s3://mod-datadumps/`

**Organization:** Data is organized by release version, then by data type, then by model organism

```
s3://mod-datadumps/
├── 8.3.0/                                    # Latest release
│   ├── DISEASE-ALLIANCE/                     # Disease associations
│   │   ├── COMBINED/                         # All species combined
│   │   ├── FB/                              # Drosophila (FlyBase)
│   │   ├── MGI/                             # Mouse
│   │   ├── RGD/                             # Rat
│   │   ├── SGD/                             # Yeast
│   │   ├── WB/                              # Worm
│   │   ├── XBXL/                            # Xenopus laevis
│   │   ├── XBXT/                            # Xenopus tropicalis
│   │   ├── ZFIN/                            # Zebrafish
│   │   └── HUMAN/                           # Human
│   ├── DISEASE-ALLIANCE-JSON/               # Same structure, JSON format
│   ├── EXPRESSION-ALLIANCE/                  # Expression annotations
│   ├── EXPRESSION-ALLIANCE-JSON/            # Expression in JSON
│   ├── GENE-DESCRIPTION-JSON/               # Gene descriptions (JSON)
│   ├── GENE-DESCRIPTION-TSV/                # Gene descriptions (TSV)
│   ├── GENE-DESCRIPTION-TXT/                # Gene descriptions (plain text)
│   ├── INTERACTION-GEN/                      # Genetic interactions
│   ├── INTERACTION-MOL/                      # Molecular interactions
│   ├── ORTHOLOGY-ALLIANCE/                   # Orthology data (TSV)
│   └── ORTHOLOGY-ALLIANCE-JSON/             # Orthology data (JSON)
├── 8.2.0/                                    # Previous release
├── 8.1.0/
├── 8.0.0/
└── variants/                                  # Separate variants directory
    └── 8.3.0/
        ├── FB/
        ├── MGI/
        ├── RGD/
        ├── WB/
        ├── XBXL/
        ├── XBXT/
        └── ZFIN/
```

### Alliance File Naming Conventions

Files follow these patterns:

**Disease Data:**
```
DISEASE-ALLIANCE_COMBINED_[number].tsv.gz
DISEASE-ALLIANCE_[MOD]_[number].tsv.gz
```

Examples:
- `DISEASE-ALLIANCE_COMBINED_2.tsv.gz` - All species disease data
- `DISEASE-ALLIANCE_FB_0.tsv.gz` - Drosophila disease data

**Cross-Species Files:**
- Expression: `EXPRESSION-ALLIANCE_COMBINED_*.tsv.gz`
- Gene Descriptions: Available in JSON, TSV, and TXT formats
- Interactions: `INTERACTION-MOL_COMBINED_*.tsv` and `INTERACTION-GEN_COMBINED_*.tsv`
- Orthology: `ORTHOLOGY-ALLIANCE_COMBINED_*.tsv`

---

## Alliance Data Categories

### 1. Disease Associations

**Location:** `s3://mod-datadumps/[VERSION]/DISEASE-ALLIANCE/`

**Available Formats:**
- TSV: `DISEASE-ALLIANCE/`
- JSON: `DISEASE-ALLIANCE-JSON/`

**Species Coverage:**
All Alliance organisms plus combined cross-species file

**Typical File Content:**
- Gene identifiers
- Disease Ontology (DO) IDs and terms
- Association types (is_model_of, etc.)
- Evidence codes
- Publications

**Use Cases:**
- Find all genes associated with a specific disease
- Identify model organism disease models
- Cross-species disease gene comparison

### 2. Expression Data

**Location:** `s3://mod-datadumps/[VERSION]/EXPRESSION-ALLIANCE/`

**Available Formats:**
- TSV: `EXPRESSION-ALLIANCE/`
- JSON: `EXPRESSION-ALLIANCE-JSON/`

**Data Types:**
- Spatiotemporal expression patterns
- Developmental stage annotations
- Anatomical location annotations
- Expression assay types

**Species Coverage:**
All Alliance model organisms

### 3. Gene Descriptions

**Location:** `s3://mod-datadumps/[VERSION]/GENE-DESCRIPTION-*/`

**Three Format Options:**
- **JSON:** `GENE-DESCRIPTION-JSON/` - Structured data
- **TSV:** `GENE-DESCRIPTION-TSV/` - Tab-delimited tables
- **TXT:** `GENE-DESCRIPTION-TXT/` - Plain text summaries

**Content:**
- Automated gene summaries
- Curated gene descriptions
- Gene symbols and IDs
- Cross-references

### 4. Molecular Interactions

**Location:** `s3://mod-datadumps/[VERSION]/INTERACTION-MOL/`

**Interaction Types:**
- Protein-protein interactions
- Protein-RNA interactions
- RNA-RNA interactions

**Format:** PSI-MI TAB compatible

**Species Coverage:**
All Alliance organisms plus SARS-CoV-2 interactions

### 5. Genetic Interactions

**Location:** `s3://mod-datadumps/[VERSION]/INTERACTION-GEN/`

**Interaction Types:**
- Suppression
- Enhancement
- Synthetic lethality/rescue

**Format:** TSV

### 6. Orthology Data

**Location:** `s3://mod-datadumps/[VERSION]/ORTHOLOGY-ALLIANCE/`

**Available Formats:**
- TSV: `ORTHOLOGY-ALLIANCE/`
- JSON: `ORTHOLOGY-ALLIANCE-JSON/`

**Content:**
- Cross-species ortholog predictions
- DIOPT scores
- Multiple algorithm consensus data

### 7. Variants (VCF Format)

**Location:** `s3://mod-datadumps/variants/[VERSION]/`

**Organisms with VCF Data:**
- FB (Drosophila melanogaster)
- MGI (Mus musculus)
- RGD (Rattus norvegicus)
- WB (Caenorhabditis elegans)
- ZFIN (Danio rerio)
- XBXL (Xenopus laevis)
- XBXT (Xenopus tropicalis)

**Format:** Standard VCF 4.2

---

## Alliance Data Access Methods

### AWS CLI Access

**List available releases:**
```bash
aws s3 ls s3://mod-datadumps/ --profile ctabone
```

**List disease data types:**
```bash
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/ --profile ctabone
```

**Download combined disease data:**
```bash
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --profile ctabone
```

**Download orthology data:**
```bash
aws s3 cp s3://mod-datadumps/8.3.0/ORTHOLOGY-ALLIANCE/COMBINED/ . --recursive --profile ctabone
```

### Alliance API Access

**REST API Base:** `https://www.alliancegenome.org/api/`

**Example API Calls:**
```bash
# Get gene information
curl https://www.alliancegenome.org/api/gene/FBgn0000001

# Search genes
curl "https://www.alliancegenome.org/api/search?category=gene&q=white"

# Get disease information
curl https://www.alliancegenome.org/api/disease/DOID:162
```

---

## Part 2: FlyBase-Specific Data (Drosophila Component)

### FlyBase S3 Repository Structure

**Primary Bucket:** `s3://s3ftp.flybase.org/`

**Public HTTPS Access:** `https://s3ftp.flybase.org/`

```
s3://s3ftp.flybase.org/
├── releases/
│   ├── current/                              # Always latest release
│   │   ├── precomputed_files/
│   │   │   ├── genes/
│   │   │   ├── alleles/
│   │   │   ├── synonyms/
│   │   │   ├── expression/
│   │   │   ├── interactions/
│   │   │   ├── human_disease/
│   │   │   ├── orthologs/
│   │   │   ├── references/
│   │   │   ├── ontologies/
│   │   │   ├── stocks/
│   │   │   └── (15+ other categories)
│   │   ├── chado-xml/
│   │   ├── psql/
│   │   ├── dmel_r6.65/                      # Current genome
│   │   └── collaborators/
│   ├── FB2025_04/                            # Latest versioned
│   ├── FB2025_03/
│   ├── FB2025_02/
│   ├── FB2025_01/
│   ├── FB2024_06/
│   └── (earlier releases)
└── genomes/
    └── Drosophila_melanogaster/
        ├── current/
        │   ├── fasta/
        │   ├── gff/
        │   └── gtf/
        └── dmel_r6.65_FB2025_04/
```

### FlyBase File Naming Conventions

**Current Release Files:**
```
<data_type>_<details>_current.tsv.gz
```

Examples:
- `fbgn_annotation_ID_current.tsv.gz`
- `gene_rpkm_matrix_current.tsv.gz`

**Versioned Release Files:**
```
<data_type>_<details>_fb_YYYY_MM.tsv.gz
```

Examples:
- `fbgn_annotation_ID_fb_2025_04.tsv.gz`
- `disease_model_annotations_fb_2025_04.tsv.gz`

**Genome Files:**
```
dmel-all-<feature_type>-r<version>.fasta.gz
dmel-all-r<version>.gff.gz
```

Examples:
- `dmel-all-chromosome-r6.65.fasta.gz`
- `dmel-all-CDS-r6.65.fasta.gz`
- `dmel-all-r6.65.gff.gz`

---

## FlyBase Data Categories

### 1. Gene Annotations

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/`

**Key Files:**
- `fbgn_annotation_ID_*.tsv.gz` - Gene ID to annotation ID (CG numbers)
- `fbgn_NAseq_Uniprot_*.tsv.gz` - Database cross-references
- `gene_map_table_*.tsv` - Genetic and cytological positions
- `best_gene_summary_*.tsv` - Curated gene descriptions

**Use Cases:**
- Convert gene symbols to FBgn IDs
- Get annotation IDs (CG numbers)
- Retrieve gene descriptions

### 2. Gene Expression Data

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/`

**Expression Files:**
- `gene_rpkm_matrix_*.tsv.gz` - RNA-Seq expression matrix
- `gene_rpkm_report_*.tsv.gz` - Detailed RPKM values
- `scRNA-Seq_gene_expression_*.tsv.gz` - Single-cell expression
- `curated_expression_*.tsv.gz` - Spatiotemporal patterns

**Data Types:**
- Bulk RNA-Seq (RPKM/FPKM values)
- Single-cell RNA-Seq (cell cluster expression)
- FlyAtlas2 data (FPKM values)
- Curated developmental expression

### 3. Disease and Phenotype Data

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/human_disease/`

**Key Files:**
- `disease_model_annotations_*.tsv.gz` - Disease Ontology annotations
- `human_disease_models_*.tsv.gz` - Detailed disease models
- `genotype_phenotype_data_*.tsv` - Phenotypic annotations

### 4. Physical Interactions

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/`

**File:** `physical_interactions_mitab_*.tsv.gz`

**Format:** PSI-MI TAB 2.7

**Content:**
- Protein-protein interactions
- Protein-RNA interactions
- Detection methods
- Publications

### 5. Genetic Interactions

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/`

**Files:**
- `gene_genetic_interactions_*.tsv` - Gene-level interactions
- `allele_genetic_interactions_*.tsv` - Allele-level interactions

### 6. Orthology Data

**Location:** `s3://s3ftp.flybase.org/releases/current/precomputed_files/orthologs/`

**Key Files:**
- `dmel_human_orthologs_disease_*.tsv.gz` - Human orthologs with disease info
- `dmel_paralogs_*.tsv.gz` - Drosophila paralogs

### 7. Genome Sequences

**Location:** `s3://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta/`

**Available Sequences:**
- Chromosomes: `dmel-all-chromosome-*.fasta.gz`
- CDS: `dmel-all-CDS-*.fasta.gz`
- Transcripts: `dmel-all-transcript-*.fasta.gz`
- Proteins: `dmel-all-translation-*.fasta.gz`
- Exons, introns, ncRNAs, tRNAs, miRNAs, etc.

### 8. Genome Annotations

**Location:** `s3://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/`

**GFF3 Files:**
- `gff/dmel-all-r*.gff.gz` - Complete genome annotation
- `gff/dmel-all-filtered-r*.gff.gz` - High-confidence annotations

**GTF Files:**
- `gtf/dmel-all-r*.gtf.gz` - Gene Transfer Format

---

## FlyBase Data Access Methods

### Public HTTPS Access (No Authentication)

**wget Examples:**
```bash
# Download gene annotations
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz

# Download expression matrix
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/gene_rpkm_matrix_current.tsv.gz

# Download genome sequence
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.65.fasta.gz
```

### AWS CLI Access (With Profile)

```bash
# List FlyBase releases
aws s3 ls s3://s3ftp.flybase.org/releases/ --profile chris

# List current precomputed files
aws s3 ls s3://s3ftp.flybase.org/releases/current/precomputed_files/ --profile chris

# Download specific file
aws s3 cp s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz . --profile chris
```

### PostgreSQL Database Access

**Public Read-Only Database:**
```bash
psql -h chado.flybase.org -U flybase flybase
```

**Connection Details:**
- Host: `chado.flybase.org`
- Database: `flybase`
- User: `flybase`
- Password: (none - read-only)
- Schema: Chado

---

## File Format Details

### Compressed Files

All data files use **GNU gzip** compression (`.gz` extension)

**Decompression:**
```bash
gunzip file.tsv.gz              # Decompress
gunzip -k file.tsv.gz           # Keep original
zcat file.tsv.gz | head         # View without decompressing
```

### TSV File Structure

**Header Format:**
```
# Comment lines start with #
# Multiple metadata lines
# Column descriptions
#
column1	column2	column3
data1	data2	data3
```

**Field Separators:**
- Tab-delimited (`\t`)
- Pipe-separated (`|`) for multiple values within fields

### JSON Format

Standard JSON or JSON Lines format depending on file

### VCF Format

Standard Variant Call Format 4.2

### PSI-MI TAB Format

15-column tab-delimited format for molecular interactions

---

## Common Use Cases

### Cross-Species Disease Analysis

**Scenario:** Find all model organism genes associated with Alzheimer's disease

**Alliance Files:**
```bash
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --profile ctabone
```

Filter for DOID:10652 (Alzheimer's disease)

### Drosophila Gene Expression Analysis

**Scenario:** Analyze developmental expression patterns

**FlyBase Files:**
```bash
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/gene_rpkm_matrix_current.tsv.gz
```

### Orthology Mapping

**Alliance Level:**
```bash
aws s3 cp s3://mod-datadumps/8.3.0/ORTHOLOGY-ALLIANCE/COMBINED/ . --recursive --profile ctabone
```

**FlyBase Level:**
```bash
wget https://s3ftp.flybase.org/releases/current/precomputed_files/orthologs/dmel_human_orthologs_disease_current.tsv.gz
```

---

## License and Citation

### Data License

Alliance data is available under **CC0 1.0 Universal (Public Domain Dedication)** and **CC-BY 4.0**.

**License Details:** https://www.alliancegenome.org/terms-of-use

### Citation

**Alliance Consortium:**
> Alliance of Genome Resources Consortium. Alliance of Genome Resources Portal: unified model organism research platform. *Nucleic Acids Research* (2023). https://doi.org/10.1093/nar/gkac1003

### Attribution

When using Alliance data:
1. Cite the Alliance consortium paper
2. Include release version numbers
3. Acknowledge data sources
4. Link to https://www.alliancegenome.org

---

## Support and Documentation

### Resources

- **Alliance Homepage:** https://www.alliancegenome.org
- **Downloads Page:** https://www.alliancegenome.org/downloads
- **API Documentation:** https://www.alliancegenome.org/api/swagger-ui
- **FlyBase:** https://flybase.org
- **GitHub:** https://github.com/alliance-genome

### Help

- **Email:** help@alliancegenome.org
- **GitHub Issues:** https://github.com/alliance-genome/agr_java_software/issues

### Tutorials

See [TUTORIAL.md](./TUTORIAL.md) for step-by-step guides.

---

## Quick Reference

### Alliance Files

| Data Type | Location | Pattern |
|-----------|----------|---------|
| Disease | `mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/` | `DISEASE-ALLIANCE_COMBINED_*.tsv.gz` |
| Expression | `mod-datadumps/8.3.0/EXPRESSION-ALLIANCE/COMBINED/` | `EXPRESSION-ALLIANCE_COMBINED_*.tsv.gz` |
| Orthology | `mod-datadumps/8.3.0/ORTHOLOGY-ALLIANCE/COMBINED/` | `ORTHOLOGY-ALLIANCE_COMBINED_*.tsv` |
| Variants | `mod-datadumps/variants/8.3.0/FB/` | VCF files |

### FlyBase Files

| Data Type | Location | Pattern |
|-----------|----------|---------|
| Gene IDs | `releases/current/precomputed_files/genes/` | `fbgn_annotation_ID_current.tsv.gz` |
| Expression | `releases/current/precomputed_files/genes/` | `gene_rpkm_matrix_current.tsv.gz` |
| Disease | `releases/current/precomputed_files/human_disease/` | `disease_model_annotations_*.tsv.gz` |
| Genome | `genomes/Drosophila_melanogaster/current/fasta/` | `dmel-all-chromosome-*.fasta.gz` |

---

**Document Version:** 2.0 (Revised)
**Last Updated:** 2025-10-17
**Verified Against:** AWS S3 buckets (mod-datadumps and s3ftp.flybase.org)
