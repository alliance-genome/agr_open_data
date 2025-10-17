# AWS S3 Structure Summary - Alliance and FlyBase

**Date:** 2025-10-17
**Purpose:** Document actual S3 bucket structures for documentation updates

---

## Alliance S3 Buckets (Profile: ctabone)

### Primary Data Bucket: `s3://mod-datadumps/`

**Versioned Releases:**
- Latest: `8.3.0/`
- Previous: `8.2.0/`, `8.1.0/`, `8.0.0/`, `7.5.0/`, etc.

**Directory Structure within each version (e.g., `8.3.0/`):**

```
s3://mod-datadumps/8.3.0/
├── DISEASE-ALLIANCE/
│   ├── COMBINED/
│   ├── FB/          (FlyBase/Drosophila)
│   ├── HUMAN/
│   ├── MGI/         (Mouse)
│   ├── RGD/         (Rat)
│   ├── SGD/         (Yeast)
│   ├── WB/          (Worm)
│   ├── XBXL/        (Xenopus laevis)
│   ├── XBXT/        (Xenopus tropicalis)
│   └── ZFIN/        (Zebrafish)
├── DISEASE-ALLIANCE-JSON/
│   └── (same subdirectories as above)
├── EXPRESSION-ALLIANCE/
│   └── (same subdirectories)
├── EXPRESSION-ALLIANCE-JSON/
│   └── (same subdirectories)
├── GENE-DESCRIPTION-JSON/
│   └── (same subdirectories)
├── GENE-DESCRIPTION-TSV/
│   └── (same subdirectories)
├── GENE-DESCRIPTION-TXT/
│   └── (same subdirectories)
├── INTERACTION-GEN/         (Genetic interactions)
│   └── (same subdirectories)
├── INTERACTION-MOL/         (Molecular interactions)
│   └── (same subdirectories)
├── ORTHOLOGY-ALLIANCE/
│   └── (combined orthology data)
└── ORTHOLOGY-ALLIANCE-JSON/
    └── (same data in JSON format)
```

**Additional directories (partial list from 49 total):**
- AGM/
- ALLELE-GFF/
- ALLELE/
- BGI/
- BIOGRID-ORCS/
- CONSTRUCT/
- CROSSREFERENCEUNIPROT/
- DAF/
- GAF/
- GENECROSSREFERENCE/
- GFF/
- HTPDATASAMPLE/
- HTPDATASET/
- HTPOSTVEPVCF/
- HTVCF/
- Human-genes-interacting-with-JSON/
- Human-genes-interacting-with/
- ... (and more)

### Variants Bucket Structure: `s3://mod-datadumps/variants/`

```
s3://mod-datadumps/variants/
├── 8.3.0/
│   ├── FB/
│   ├── HUMAN/
│   ├── MGI/
│   ├── RGD/
│   ├── WB/
│   ├── XBXL/
│   ├── XBXT/
│   └── ZFIN/
├── 8.2.0/
└── (earlier versions)
```

**File Naming Pattern (Alliance):**
- `DISEASE-ALLIANCE_COMBINED_0.tsv.gz`
- `DISEASE-ALLIANCE_COMBINED_1.tsv.gz`
- (Files appear to be numbered, possibly split by size)

---

## FlyBase S3 Buckets (Profile: chris)

### Primary FTP Bucket: `s3://s3ftp.flybase.org/`

**Top-Level Structure:**
```
s3://s3ftp.flybase.org/
├── releases/
│   ├── current/              # Always points to latest release
│   ├── FB2025_04/            # Latest versioned release
│   ├── FB2025_03/
│   ├── FB2025_02/
│   ├── FB2025_01/
│   ├── FB2024_06/
│   ├── FB2024_05/
│   ├── FB2024_04/
│   ├── FB2024_03/
│   ├── FB2023_05/
│   ├── FB2017_05/
│   ├── FB2014_03/
│   └── dmel_REFSEQ/
├── genomes/
│   ├── Drosophila_melanogaster/
│   └── dmel/
├── alliance/
├── flybase/
├── flybase-data/
├── docs/
├── maps/
├── news/
├── ReadMe/
└── icons/
```

**Current Release Structure (`releases/current/`):**
```
s3://s3ftp.flybase.org/releases/current/
├── precomputed_files/
│   ├── genes/
│   ├── alleles/
│   ├── aberrations/
│   ├── chemicals/
│   ├── experimental_tools/
│   ├── go/
│   ├── human_disease/
│   ├── insertions/
│   ├── map_conversion/
│   ├── metadata/
│   ├── ontologies/
│   ├── orthologs/
│   ├── reagents/
│   ├── references/
│   ├── species/
│   ├── stocks/
│   ├── synonyms/
│   ├── transposons/
│   └── collaborators/
├── chado-xml/
├── psql/
├── dmel_r6.65/           # Current genome release
└── collaborators/
```

**Genomes Structure:**
```
s3://s3ftp.flybase.org/genomes/
├── Drosophila_melanogaster/
│   ├── current/
│   ├── dmel_r6.65_FB2025_04/
│   ├── dmel_r6.64_FB2025_03/
│   └── (earlier genome releases)
└── dmel/
```

**File Naming Patterns (FlyBase):**
- Current: `fbgn_annotation_ID_current.tsv.gz`
- Versioned: `fbgn_annotation_ID_fb_2025_04.tsv.gz`
- Genome: `dmel-all-chromosome-r6.65.fasta.gz`

---

## Key Differences from Initial Documentation

### Incorrect Assumptions:
1. - Bucket name was `alliance-genome-data` → - Actually `mod-datadumps`
2. - Alliance structure mirrored FlyBase → - Alliance has completely different structure
3. - FlyBase was at `s3://alliance-genome-data/` → - Actually at `s3://s3ftp.flybase.org/`
4. - Combined Alliance+FlyBase structure → - Separate buckets with different organizations

### Correct Structure:
1. - Alliance: `s3://mod-datadumps/` with version directories containing data type subdirectories
2. - FlyBase: `s3://s3ftp.flybase.org/releases/` with release directories containing category subdirectories

---

## Access Patterns

### Alliance Data Access:
```bash
# List Alliance releases
aws s3 ls s3://mod-datadumps/ --profile ctabone

# List disease data types
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/ --profile ctabone

# Download combined disease data
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --profile ctabone
```

### FlyBase Data Access:
```bash
# List FlyBase releases
aws s3 ls s3://s3ftp.flybase.org/releases/ --profile chris

# List current precomputed files
aws s3 ls s3://s3ftp.flybase.org/releases/current/precomputed_files/ --profile chris

# Download gene annotations
aws s3 cp s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz . --profile chris
```

---

## Public Access Notes

**Important:** The documentation should note that:
1. Alliance data in `mod-datadumps` may require AWS credentials (not anonymous)
2. FlyBase data at `s3ftp.flybase.org` is publicly accessible via HTTPS
3. Public URL pattern for FlyBase: `https://s3ftp.flybase.org/releases/current/...`

---

## Documentation Updates Required

1. **DATA_DOCUMENTATION.md:**
   - Section 1: Alliance data structure (mod-datadumps)
   - Section 2: FlyBase data structure (s3ftp.flybase.org)
   - Fix all S3 bucket references
   - Update file path examples

2. **TUTORIAL.md:**
   - Update all download examples
   - Use correct bucket names
   - Provide both Alliance and FlyBase examples

3. **README.md:**
   - Update quick start examples
   - Fix S3 access commands
   - Correct file path references

4. **Registry YAML:**
   - Use correct bucket names
   - Provide accurate data access patterns
