# Alliance of Genome Resources - Data Usage Tutorials

This guide provides practical, step-by-step tutorials for accessing and using Alliance of Genome Resources data on AWS. Each tutorial demonstrates a specific use case with complete, working examples.

---

## Table of Contents

1. [Tutorial 1: Getting Started - Downloading Your First Dataset](#tutorial-1-getting-started---downloading-your-first-dataset)
2. [Tutorial 2: Gene Annotation Lookups and ID Conversion](#tutorial-2-gene-annotation-lookups-and-id-conversion)
3. [Tutorial 3: Disease Gene Discovery and Analysis](#tutorial-3-disease-gene-discovery-and-analysis)
4. [Tutorial 4: RNA-Seq Expression Analysis Across Developmental Stages](#tutorial-4-rna-seq-expression-analysis-across-developmental-stages)
5. [Tutorial 5: Building Protein Interaction Networks](#tutorial-5-building-protein-interaction-networks)
6. [Tutorial 6: Cross-Species Orthology Mapping](#tutorial-6-cross-species-orthology-mapping)

---

## Understanding Alliance Data Organization

The Alliance of Genome Resources provides data through two main S3 buckets:

### 1. Alliance-Wide Integrated Data
- **Location:** `s3://mod-datadumps/`
- **Content:** Cross-species integrated datasets (disease associations, expression, orthology, interactions)
- **Access:** Public (use `--no-sign-request` with AWS CLI)
- **Organization:** By release version (8.3.0/, 8.2.0/, etc.)

### 2. FlyBase-Specific Data (Drosophila Component)
- **Location:** `s3://s3ftp.flybase.org/`
- **Content:** Comprehensive Drosophila annotations, genomes, and precomputed files
- **Access:** Public HTTPS or AWS CLI (use `--no-sign-request`)
- **Organization:** By release (current/, FB2025_04/, etc.)

Most tutorials will use **FlyBase data** for organism-specific examples, with **Alliance data** for cross-species analyses.

---

## Prerequisites

### Tools and Setup

**Required:**
- Terminal/command line access
- `wget` or `curl` (for downloading files)
- `gunzip` (for decompressing .gz files)
- Text editor or spreadsheet software

**Optional but Recommended:**
- Python 3.7+ with pandas, boto3, matplotlib
- AWS CLI (for S3 access)
- R with tidyverse (for data analysis)
- Cytoscape (for network visualization)

### Installation Commands

**Python Dependencies:**
```bash
pip install pandas boto3 matplotlib seaborn biopython
```

**AWS CLI:**
```bash
# Installation varies by platform
# See: https://aws.amazon.com/cli/

# No AWS account or configuration needed - all data is publicly accessible
```

---

## Tutorial 1: Getting Started - Downloading Your First Dataset

**Goal:** Learn how to navigate both Alliance and FlyBase data structures and download gene annotation data.

**Time:** 10 minutes

### Part A: FlyBase-Specific Data (Public Access)

FlyBase data is publicly accessible via HTTPS - no AWS credentials needed!

#### Step 1: Explore Available FlyBase Data

```bash
# List top-level directories (via wget HTML parsing)
wget -qO- https://s3ftp.flybase.org/releases/current/ | grep -o 'href="[^"]*"' | grep -o '[^"]*/$'
```

Expected output shows directories like:
- `precomputed_files/`
- `chado-xml/`
- `psql/`
- `dmel_r6.65/` (genome release)

#### Step 2: Browse FlyBase Gene Data Files

```bash
# List gene-related files
wget -qO- https://s3ftp.flybase.org/releases/current/precomputed_files/genes/ | grep '.tsv.gz'
```

#### Step 3: Download FlyBase Gene Annotation File

```bash
# Download gene annotation IDs
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz

# Decompress the file
gunzip fbgn_annotation_ID_current.tsv.gz
```

#### Step 4: Inspect the Data

```bash
# View the first 20 lines (including headers)
head -n 20 fbgn_annotation_ID_current.tsv

# Count total genes
grep -v '^#' fbgn_annotation_ID_current.tsv | wc -l
```

#### Step 5: Alternative - Using AWS CLI for FlyBase

```bash
# List files using AWS CLI (no credentials needed)
aws s3 ls s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/ --no-sign-request

# Download using AWS CLI
aws s3 cp s3://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz . --no-sign-request
```

### Part B: Alliance-Wide Data (Public Access)

Alliance integrated data is also publicly accessible.

#### Step 6: Explore Alliance Data Structure

```bash
# List available Alliance releases
aws s3 ls s3://mod-datadumps/ --no-sign-request

# Expected: 8.3.0/, 8.2.0/, 8.1.0/, etc.
```

#### Step 7: Browse Alliance Disease Data

```bash
# List disease data categories in latest release
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/ --no-sign-request

# Expected: COMBINED/, FB/, MGI/, RGD/, SGD/, WB/, XBXL/, XBXT/, ZFIN/, HUMAN/
```

#### Step 8: Download Alliance Combined Disease Data

```bash
# Download combined disease associations across all species
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --no-sign-request

# Decompress
gunzip DISEASE-ALLIANCE_COMBINED_2.tsv.gz

# Inspect
head -n 20 DISEASE-ALLIANCE_COMBINED_2.tsv
```

### What You Learned

- How to access FlyBase data via public HTTPS
- How to access Alliance integrated data via AWS CLI (no credentials needed)
- How to download and decompress data files
- File naming conventions (`_current` for FlyBase, numbered files for Alliance)
- Understanding the two-tier data organization

---

## Tutorial 2: Gene Annotation Lookups and ID Conversion

**Goal:** Convert gene symbols to FlyBase IDs and retrieve gene annotations using FlyBase data.

**Time:** 15 minutes

### Scenario

You have a list of gene symbols (e.g., `white`, `ebony`, `yellow`) and need to:
1. Find their FlyBase IDs (FBgn numbers)
2. Get their annotation IDs (CG numbers)
3. Retrieve full gene names

### Step 1: Download Required FlyBase Files

```bash
# Download gene ID mappings
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz
gunzip fbgn_annotation_ID_current.tsv.gz

# Download gene synonyms (includes full names)
wget https://s3ftp.flybase.org/releases/current/precomputed_files/synonyms/fb_synonym_current.tsv.gz
gunzip fb_synonym_current.tsv.gz
```

### Step 2: Look Up a Single Gene (Command Line)

```bash
# Find the gene 'white' in the synonyms file
grep -i "^FBgn[0-9]*\tDmel\twhite\t" fb_synonym_current.tsv

# Expected output:
# FBgn0003996	Dmel	w	white	...
```

### Step 3: Batch Conversion with Python

Save this as `gene_lookup.py`:

```python
import pandas as pd
import sys

# Load data files
print("Loading gene annotation data...")
annotations = pd.read_csv('fbgn_annotation_ID_current.tsv', sep='\t', comment='#')
synonyms = pd.read_csv('fb_synonym_current.tsv', sep='\t', comment='#')

# List of genes to look up
genes_to_find = ['white', 'ebony', 'yellow', 'sevenless', 'Notch']

print(f"\nLooking up {len(genes_to_find)} genes...\n")

results = []
for gene_symbol in genes_to_find:
    # Find in synonyms file
    match = synonyms[synonyms['current_symbol'] == gene_symbol]
    
    if not match.empty:
        fbgn = match.iloc[0]['primary_FBid']
        fullname = match.iloc[0]['current_fullname']
        
        # Get annotation ID
        annot_match = annotations[annotations['gene_symbol'] == gene_symbol]
        cg_id = annot_match.iloc[0]['annotation_ID'] if not annot_match.empty else 'N/A'
        
        results.append({
            'Symbol': gene_symbol,
            'FBgn': fbgn,
            'CG_ID': cg_id,
            'Full_Name': fullname
        })
        
        print(f"- {gene_symbol:12} → {fbgn:12} → {cg_id:10} → {fullname}")
    else:
        print(f"- {gene_symbol:12} → Not found")

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv('gene_lookup_results.csv', index=False)
print(f"\nResults saved to gene_lookup_results.csv")
```

Run the script:

```bash
python gene_lookup.py
```

**Expected Output:**
```
Loading gene annotation data...

Looking up 5 genes...

- white        → FBgn0003996  → CG2759     → white
- ebony        → FBgn0000527  → CG3331     → ebony
- yellow       → FBgn0004034  → CG3879     → yellow
- sevenless    → FBgn0003366  → CG18085    → sevenless
- Notch        → FBgn0004647  → CG3936     → Notch

Results saved to gene_lookup_results.csv
```

### Step 4: Reverse Lookup (FBgn to Symbol)

```python
# Reverse lookup: FBgn → Symbol
fbgn_list = ['FBgn0003996', 'FBgn0000527', 'FBgn0004034']

for fbgn in fbgn_list:
    match = synonyms[synonyms['primary_FBid'] == fbgn]
    if not match.empty:
        symbol = match.iloc[0]['current_symbol']
        print(f"{fbgn} → {symbol}")
```

### What You Learned

- How to perform gene symbol → FBgn ID conversion
- How to retrieve annotation IDs (CG numbers)
- Batch lookup techniques using Python pandas
- Working with multiple FlyBase data files together

---

## Tutorial 3: Disease Gene Discovery and Analysis

**Goal:** Find all genes associated with a specific human disease using both Alliance-wide and FlyBase-specific data.

**Time:** 20 minutes

### Scenario

You're researching Alzheimer's disease (DOID:10652) and want to:
1. Find all Alliance genes associated with this disease across all species
2. Focus on FlyBase (Drosophila) disease models
3. Get human orthologs for cross-species comparison

### Part A: Alliance-Wide Disease Data

#### Step 1: Download Alliance Disease Association Data

```bash
# Download combined disease data across all Alliance species
aws s3 cp s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/COMBINED/DISEASE-ALLIANCE_COMBINED_2.tsv.gz . --no-sign-request
gunzip DISEASE-ALLIANCE_COMBINED_2.tsv.gz

# Browse FlyBase-specific disease data from Alliance
aws s3 ls s3://mod-datadumps/8.3.0/DISEASE-ALLIANCE/FB/ --no-sign-request
```

#### Step 2: Search Alliance Data for Disease-Associated Genes

Save as `alliance_disease_search.py`:

```python
import pandas as pd

# Load Alliance combined disease annotations
print("Loading Alliance disease association data...")
disease_data = pd.read_csv('DISEASE-ALLIANCE_COMBINED_2.tsv',
                          sep='\t', comment='#', low_memory=False)

print(f"Loaded {len(disease_data)} disease associations")
print(f"\nColumn names: {list(disease_data.columns)}")

# Search for Alzheimer's disease across all Alliance species
alzheimers = disease_data[disease_data.apply(
    lambda row: 'DOID:10652' in str(row.values) or 
                'Alzheimer' in str(row.values), axis=1)]

print(f"\nFound {len(alzheimers)} gene associations with Alzheimer's disease across Alliance\n")

# Group by species if species column exists
if 'species' in disease_data.columns or 'taxon' in disease_data.columns:
    species_col = 'species' if 'species' in disease_data.columns else 'taxon'
    species_counts = alzheimers[species_col].value_counts()
    print("Associations by species:")
    print(species_counts)

# Show first 10 associations
print("\nFirst 10 associated genes:")
print(alzheimers.head(10))

# Save results
alzheimers.to_csv('alliance_alzheimers_genes.csv', index=False)
print(f"\nFull results saved to alliance_alzheimers_genes.csv")
```

### Part B: FlyBase-Specific Disease Data

#### Step 3: Download FlyBase Disease Data

```bash
# Download FlyBase disease model annotations
wget https://s3ftp.flybase.org/releases/current/precomputed_files/human_disease/disease_model_annotations_fb_current.tsv.gz
gunzip disease_model_annotations_fb_current.tsv.gz
```

#### Step 4: Search FlyBase Data

Save as `flybase_disease_search.py`:

```python
import pandas as pd

# Load FlyBase disease annotations
print("Loading FlyBase disease association data...")
disease_data = pd.read_csv('disease_model_annotations_fb_current.tsv',
                          sep='\t', comment='#')

# Search for Alzheimer's disease (DOID:10652)
alzheimers = disease_data[disease_data['DOid'].str.contains('DOID:10652', na=False) |
                         disease_data['DOterm'].str.contains('Alzheimer', case=False, na=False)]

print(f"\nFound {len(alzheimers)} gene associations with Alzheimer's disease in FlyBase\n")

# Show gene details
print("Associated Drosophila genes:")
print(alzheimers[['DBobjectSymbol', 'DBobjectID', 'DOterm',
                  'EvidenceCode', 'AssociationType']].to_string(index=False))

# Save results
alzheimers.to_csv('flybase_alzheimers_genes.csv', index=False)
print(f"\nFull results saved to flybase_alzheimers_genes.csv")
```

Run it:

```bash
python flybase_disease_search.py
```

#### Step 5: Find Human Orthologs for Disease Genes

```bash
# Download FlyBase human orthology data
wget https://s3ftp.flybase.org/releases/current/precomputed_files/orthologs/dmel_human_orthologs_disease_fb_current.tsv.gz
gunzip dmel_human_orthologs_disease_fb_current.tsv.gz
```

```python
import pandas as pd

# Load orthology data
orthologs = pd.read_csv('dmel_human_orthologs_disease_fb_current.tsv',
                       sep='\t', comment='#')

# Load our Alzheimer's genes (Drosophila)
alzheimers_dmel = pd.read_csv('flybase_alzheimers_genes.csv')
dmel_genes = alzheimers_dmel['DBobjectID'].tolist()

# Find human orthologs
human_orthologs = orthologs[orthologs['FBgn_ID'].isin(dmel_genes)]

print(f"\nFound {len(human_orthologs)} Drosophila genes with human orthologs\n")
print(human_orthologs[['Dmel_gene_symbol', 'Human_gene_symbol',
                       'DIOPT_score', 'OMIM_Phenotype_IDs']].to_string(index=False))

human_orthologs.to_csv('alzheimers_human_orthologs.csv', index=False)
```

#### Step 6: Visualize Disease Gene Distribution

```python
import matplotlib.pyplot as plt
import pandas as pd

# Load data
disease_data = pd.read_csv('disease_model_annotations_fb_current.tsv',
                          sep='\t', comment='#')

# Count top 10 diseases by gene associations
top_diseases = disease_data['DOterm'].value_counts().head(10)

# Create bar plot
plt.figure(figsize=(12, 6))
top_diseases.plot(kind='barh')
plt.xlabel('Number of Gene Associations')
plt.ylabel('Disease')
plt.title('Top 10 Diseases by Gene Associations in FlyBase')
plt.tight_layout()
plt.savefig('top_diseases_flybase.png', dpi=300)
print("Chart saved to top_diseases_flybase.png")
```

### What You Learned

- How to search Alliance-wide disease associations across all species
- How to use FlyBase-specific disease model data
- Cross-species disease model analysis
- Mapping model organism genes to human orthologs
- Visualization of disease-gene relationships

---

## Tutorial 4: RNA-Seq Expression Analysis Across Developmental Stages

**Goal:** Analyze gene expression patterns across Drosophila developmental stages using FlyBase RNA-Seq data.

**Time:** 25 minutes

### Scenario

You want to understand how a set of genes (e.g., developmental regulators) are expressed across different life stages.

### Step 1: Download FlyBase Expression Data

```bash
# Download RNA-Seq expression matrix
wget https://s3ftp.flybase.org/releases/current/precomputed_files/genes/gene_rpkm_matrix_fb_current.tsv.gz
gunzip gene_rpkm_matrix_fb_current.tsv.gz
```

### Step 2: Load and Explore the Expression Matrix

Save as `expression_analysis.py`:

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load expression matrix
print("Loading expression matrix (this may take a minute)...")
expr = pd.read_csv('gene_rpkm_matrix_fb_current.tsv',
                  sep='\t', comment='#', low_memory=False)

print(f"Matrix shape: {expr.shape[0]} genes × {expr.shape[1]} columns")
print(f"\nFirst few columns: {list(expr.columns[:10])}")

# Show sample names (dataset columns)
sample_cols = [col for col in expr.columns if not col.startswith(('gene_', 'FBgn'))]
print(f"\nTotal expression samples: {len(sample_cols)}")
print("\nFirst 5 samples:")
for col in sample_cols[:5]:
    print(f"  - {col}")
```

### Step 3: Filter for Developmental Stage Samples

```python
# Define genes of interest (developmental regulators)
genes_of_interest = ['even-skipped', 'fushi tarazu', 'engrailed', 'wingless', 'hedgehog']

# Filter matrix for these genes
expr_subset = expr[expr['gene_symbol'].isin(genes_of_interest)].copy()

print(f"\nFound {len(expr_subset)} of {len(genes_of_interest)} genes in dataset")

# Select only embryonic and larval stage samples (example pattern)
# Note: Actual column names will vary - inspect your data first!
stage_columns = [col for col in expr.columns if 'embryo' in col.lower() or
                'larva' in col.lower() or 'adult' in col.lower()]

if len(stage_columns) == 0:
    print("No stage-specific columns found. Using first 20 sample columns.")
    stage_columns = sample_cols[:20]

print(f"\nAnalyzing {len(stage_columns)} developmental samples")
```

### Step 4: Create Expression Heatmap

```python
# Prepare data for heatmap
heatmap_data = expr_subset[['gene_symbol'] + stage_columns].set_index('gene_symbol')

# Convert to numeric, handling any non-numeric values
heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')

# Log-transform expression values (adding pseudocount)
heatmap_data_log = np.log2(heatmap_data + 1)

# Create heatmap
plt.figure(figsize=(16, 6))
sns.heatmap(heatmap_data_log,
            cmap='YlOrRd',
            cbar_kws={'label': 'log2(RPKM + 1)'},
            yticklabels=True,
            xticklabels=True)
plt.title('Developmental Gene Expression Across Stages', fontsize=14)
plt.xlabel('Developmental Stage / Sample', fontsize=10)
plt.ylabel('Gene', fontsize=12)
plt.xticks(rotation=90, fontsize=8)
plt.tight_layout()
plt.savefig('expression_heatmap.png', dpi=300)
print("\nHeatmap saved to expression_heatmap.png")
```

### Step 5: Plot Individual Gene Expression Profiles

```python
# Create line plots for each gene
fig, axes = plt.subplots(len(genes_of_interest), 1, figsize=(14, 10), sharex=True)

for idx, gene in enumerate(genes_of_interest):
    gene_data = expr_subset[expr_subset['gene_symbol'] == gene][stage_columns]
    
    if not gene_data.empty:
        values = gene_data.iloc[0].values
        values_numeric = pd.to_numeric(values, errors='coerce')
        
        axes[idx].plot(range(len(values_numeric)), values_numeric, marker='o', linewidth=2)
        axes[idx].set_ylabel('RPKM', fontsize=10)
        axes[idx].set_title(f'{gene} Expression', fontsize=11)
        axes[idx].grid(True, alpha=0.3)

axes[-1].set_xlabel('Sample Index', fontsize=10)
plt.tight_layout()
plt.savefig('gene_expression_profiles.png', dpi=300)
print("Gene profiles saved to gene_expression_profiles.png")
```

### Step 6: Find Highly Expressed Genes in Specific Stages

```python
# Find top expressed genes in embryonic stage (example)
embryo_cols = [col for col in stage_columns if 'embryo' in col.lower()]

if len(embryo_cols) > 0:
    # Calculate mean expression across embryonic samples
    expr['embryo_mean'] = expr[embryo_cols].apply(pd.to_numeric, errors='coerce').mean(axis=1)
    
    # Sort and get top 20
    top_embryo = expr.nlargest(20, 'embryo_mean')[['gene_symbol', 'gene_fullname', 'embryo_mean']]
    
    print("\nTop 20 genes expressed in embryonic stages:")
    print(top_embryo.to_string(index=False))
    
    top_embryo.to_csv('top_embryo_genes.csv', index=False)
```

### What You Learned

- How to work with large RNA-Seq expression matrices from FlyBase
- Filtering expression data by developmental stage
- Creating expression heatmaps for visualization
- Identifying stage-specific highly expressed genes
- Log-transformation and normalization techniques

---

## Tutorial 5: Building Protein Interaction Networks

**Goal:** Download physical interaction data from both Alliance and FlyBase and build protein-protein interaction networks.

**Time:** 20 minutes

### Scenario

You want to:
1. Download protein interaction data in PSI-MI TAB format
2. Build an interaction network for a protein of interest
3. Visualize the network
4. Export for Cytoscape

### Part A: Alliance Molecular Interactions (Cross-Species)

#### Step 1: Download Alliance Interaction Data

```bash
# Download Alliance molecular interactions for Drosophila
aws s3 ls s3://mod-datadumps/8.3.0/INTERACTION-MOL/FB/ --no-sign-request
aws s3 cp s3://mod-datadumps/8.3.0/INTERACTION-MOL/FB/INTERACTION-MOL_FB_0.tsv.gz . --no-sign-request
gunzip INTERACTION-MOL_FB_0.tsv.gz
```

### Part B: FlyBase-Specific Interactions

#### Step 2: Download FlyBase Interaction Data

```bash
# Download physical interactions in MITAB format
wget https://s3ftp.flybase.org/releases/current/precomputed_files/interactions/physical_interactions_mitab_fb_current.tsv.gz
gunzip physical_interactions_mitab_fb_current.tsv.gz
```

#### Step 3: Parse MITAB Format

Save as `ppi_network.py`:

```python
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load interaction data
print("Loading physical interactions...")
interactions = pd.read_csv('physical_interactions_mitab_fb_current.tsv',
                          sep='\t', comment='#',
                          header=None)  # MITAB format has no header

# Assign column names (PSI-MI TAB 2.7 format)
col_names = ['ID_A', 'ID_B', 'Alt_ID_A', 'Alt_ID_B', 'Alias_A', 'Alias_B',
             'Interaction_Detection_Method', 'Publication_Author', 'Publication_ID',
             'Taxon_A', 'Taxon_B', 'Interaction_Type', 'Source_Database',
             'Interaction_ID', 'Confidence_Score']

interactions.columns = col_names[:len(interactions.columns)]

print(f"Loaded {len(interactions)} interactions")

# Extract gene symbols from Alias columns (format: flybase:gene(gene name))
def extract_gene_symbol(alias_str):
    if pd.isna(alias_str):
        return None
    # Parse format like "flybase:eve(gene name)"
    if '(' in alias_str:
        return alias_str.split('(')[0].split(':')[-1]
    return None

interactions['Gene_A'] = interactions['Alias_A'].apply(extract_gene_symbol)
interactions['Gene_B'] = interactions['Alias_B'].apply(extract_gene_symbol)

# Remove rows where gene symbols couldn't be extracted
interactions_clean = interactions.dropna(subset=['Gene_A', 'Gene_B'])
print(f"Cleaned to {len(interactions_clean)} interactions with gene symbols")
```

#### Step 4: Build Network for Gene of Interest

```python
# Choose a gene to analyze (e.g., 'Notch')
gene_of_interest = 'Notch'

# Find all interactions involving this gene
gene_interactions = interactions_clean[
    (interactions_clean['Gene_A'] == gene_of_interest) |
    (interactions_clean['Gene_B'] == gene_of_interest)
]

print(f"\nFound {len(gene_interactions)} interactions for {gene_of_interest}")

# Build network
G = nx.Graph()

for _, row in gene_interactions.iterrows():
    G.add_edge(row['Gene_A'], row['Gene_B'])

print(f"Network has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

# Calculate network statistics
print(f"\nNetwork statistics for {gene_of_interest}:")
print(f"  Direct interactors: {G.degree(gene_of_interest)}")
print(f"  Network density: {nx.density(G):.3f}")

# Find most connected proteins
degree_centrality = nx.degree_centrality(G)
top_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:10]
print(f"\nTop 10 most connected proteins in this network:")
for gene, centrality in top_nodes:
    print(f"  {gene:20} {centrality:.3f}")
```

#### Step 5: Visualize the Network

```python
# Visualize network
plt.figure(figsize=(14, 14))

# Position nodes using spring layout
pos = nx.spring_layout(G, k=0.5, iterations=50)

# Draw network
nx.draw_networkx_nodes(G, pos,
                       node_size=[3000 if node == gene_of_interest else 800
                                for node in G.nodes()],
                       node_color=['red' if node == gene_of_interest else 'lightblue'
                                  for node in G.nodes()],
                       alpha=0.7)

nx.draw_networkx_edges(G, pos, alpha=0.3)

nx.draw_networkx_labels(G, pos,
                        font_size=8,
                        font_family='sans-serif')

plt.title(f'Protein Interaction Network for {gene_of_interest}', fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.savefig(f'{gene_of_interest}_network.png', dpi=300, bbox_inches='tight')
print(f"\nNetwork visualization saved to {gene_of_interest}_network.png")
```

#### Step 6: Export for Cytoscape

```python
# Export edge list for Cytoscape
edge_list = gene_interactions[['Gene_A', 'Gene_B', 'Interaction_Detection_Method',
                               'Publication_ID']].copy()
edge_list.columns = ['Source', 'Target', 'Detection_Method', 'Publication']

edge_list.to_csv(f'{gene_of_interest}_network_cytoscape.csv', index=False)
print(f"Cytoscape import file saved to {gene_of_interest}_network_cytoscape.csv")

# Export node attributes
node_attributes = pd.DataFrame({
    'Gene': list(G.nodes()),
    'Degree': [G.degree(n) for n in G.nodes()],
    'Betweenness_Centrality': [nx.betweenness_centrality(G)[n] for n in G.nodes()],
    'Is_Query_Gene': [n == gene_of_interest for n in G.nodes()]
})

node_attributes.to_csv(f'{gene_of_interest}_node_attributes.csv', index=False)
print(f"Node attributes saved to {gene_of_interest}_node_attributes.csv")
```

#### Step 7: Build Second-Degree Network

```python
# Expand to second-degree interactions
second_degree_genes = set()
for neighbor in G.neighbors(gene_of_interest):
    second_degree_genes.update(G.neighbors(neighbor))

# Get all interactions among these genes
expanded_interactions = interactions_clean[
    interactions_clean['Gene_A'].isin(second_degree_genes) &
    interactions_clean['Gene_B'].isin(second_degree_genes)
]

print(f"\nExpanded network:")
print(f"  Second-degree interactors: {len(second_degree_genes)}")
print(f"  Total interactions: {len(expanded_interactions)}")
```

### What You Learned

- How to access Alliance molecular interaction data
- How to parse PSI-MI TAB format interaction data from FlyBase
- Building protein-protein interaction networks
- Network analysis (centrality, degree, density)
- Network visualization with matplotlib
- Exporting data for Cytoscape analysis

---

## Tutorial 6: Cross-Species Orthology Mapping

**Goal:** Map genes between species using both Alliance-wide and FlyBase orthology data.

**Time:** 15 minutes

### Scenario

You have identified interesting genes in one species and want to find their orthologs in other model organisms for comparative studies.

### Part A: Alliance Orthology Data (Cross-Species)

#### Step 1: Download Alliance Orthology Data

```bash
# Download Alliance combined orthology data
aws s3 ls s3://mod-datadumps/8.3.0/ORTHOLOGY-ALLIANCE/ --no-sign-request
aws s3 cp s3://mod-datadumps/8.3.0/ORTHOLOGY-ALLIANCE/ORTHOLOGY-ALLIANCE_COMBINED_0.tsv.gz . --no-sign-request
gunzip ORTHOLOGY-ALLIANCE_COMBINED_0.tsv.gz
```

#### Step 2: Explore Alliance Orthology

Save as `alliance_orthology.py`:

```python
import pandas as pd

# Load Alliance orthology data
print("Loading Alliance orthology data...")
ortho = pd.read_csv('ORTHOLOGY-ALLIANCE_COMBINED_0.tsv',
                   sep='\t', comment='#', low_memory=False)

print(f"Loaded {len(ortho)} orthology relationships")
print(f"\nColumns: {list(ortho.columns)}")

# Show species pairs
print("\nOrthology relationships summary:")
print(ortho.head(20))
```

### Part B: FlyBase-Specific Orthology (Drosophila-Human)

#### Step 3: Download FlyBase Orthology Data

```bash
# Download Drosophila-human orthology
wget https://s3ftp.flybase.org/releases/current/precomputed_files/orthologs/dmel_human_orthologs_disease_fb_current.tsv.gz
gunzip dmel_human_orthologs_disease_fb_current.tsv.gz

# Note: FlyBase provides other species pairs as well
# Check the orthologs/ directory for more options
```

#### Step 4: Simple Orthology Lookup

Save as `orthology_mapping.py`:

```python
import pandas as pd

# Load orthology data
orthologs = pd.read_csv('dmel_human_orthologs_disease_fb_current.tsv',
                       sep='\t', comment='#')

print(f"Loaded {len(orthologs)} Drosophila-Human orthology relationships")

# Look up human orthologs for Drosophila genes
dmel_genes = ['white', 'Notch', 'Toll', 'nanos', 'oskar']

print("\nDrosophila → Human Orthology:")
print("=" * 80)

for gene in dmel_genes:
    matches = orthologs[orthologs['Dmel_gene_symbol'] == gene]
    
    if not matches.empty:
        for _, row in matches.iterrows():
            print(f"\n{gene} (FBgn: {row['FBgn_ID']})")
            print(f"  → Human ortholog: {row['Human_gene_symbol']} ({row['Human_gene_ID']})")
            print(f"  → DIOPT score: {row['DIOPT_score']}")
            if pd.notna(row['OMIM_Phenotype_IDs']):
                print(f"  → Associated OMIM phenotypes: {row['OMIM_Phenotype_IDs']}")
    else:
        print(f"\n{gene}: No human ortholog found")
```

Run it:

```bash
python orthology_mapping.py
```

#### Step 5: Reverse Mapping (Human to Drosophila)

```python
# Reverse lookup: Human → Drosophila
human_genes = ['TP53', 'BRCA1', 'APP', 'HTT']

print("\n\nHuman → Drosophila Orthology:")
print("=" * 80)

for gene in human_genes:
    matches = orthologs[orthologs['Human_gene_symbol'] == gene]
    
    if not matches.empty:
        for _, row in matches.iterrows():
            print(f"\n{gene} ({row['Human_gene_ID']})")
            print(f"  → Drosophila ortholog: {row['Dmel_gene_symbol']} ({row['FBgn_ID']})")
            print(f"  → DIOPT score: {row['DIOPT_score']}")
            print(f"  → Best score: {row['Best_Score']}")
            print(f"  → Best reverse score: {row['Best_Score_Reverse']}")
    else:
        print(f"\n{gene}: No Drosophila ortholog found")
```

#### Step 6: Filter by Orthology Confidence

```python
# Filter for high-confidence orthologs (DIOPT score >= 8)
high_conf = orthologs[orthologs['DIOPT_score'] >= 8]

print(f"\n\nHigh-confidence orthologs (DIOPT >= 8): {len(high_conf)} pairs")
print(f"Percentage of total: {len(high_conf)/len(orthologs)*100:.1f}%")

# Distribution of DIOPT scores
score_dist = orthologs['DIOPT_score'].value_counts().sort_index()
print("\nDIOPT Score Distribution:")
for score, count in score_dist.items():
    print(f"  Score {score}: {count} pairs")
```

#### Step 7: Disease-Associated Orthologs

```python
# Find orthologs associated with specific diseases
# Filter for genes with OMIM phenotype associations
disease_orthologs = orthologs[orthologs['OMIM_Phenotype_IDs'].notna()]

print(f"\n\nOrthologs with disease associations: {len(disease_orthologs)}")

# Show examples
print("\nExamples of disease-associated orthologs:")
print(disease_orthologs[['Dmel_gene_symbol', 'Human_gene_symbol',
                        'OMIM_Phenotype_IDs', 'DIOPT_score']].head(10).to_string(index=False))

# Save for further analysis
disease_orthologs.to_csv('disease_associated_orthologs.csv', index=False)
print("\nSaved to disease_associated_orthologs.csv")
```

#### Step 8: Batch Orthology Conversion

```python
def batch_orthology_lookup(gene_list, species_from='Dmel', species_to='Human'):
    """
    Convert a list of genes to their orthologs
    """
    results = []
    
    for gene in gene_list:
        if species_from == 'Dmel' and species_to == 'Human':
            matches = orthologs[orthologs['Dmel_gene_symbol'] == gene]
            if not matches.empty:
                for _, row in matches.iterrows():
                    results.append({
                        'Input_Gene': gene,
                        'Input_Species': 'Drosophila',
                        'Ortholog': row['Human_gene_symbol'],
                        'Ortholog_Species': 'Human',
                        'Ortholog_ID': row['Human_gene_ID'],
                        'DIOPT_Score': row['DIOPT_score'],
                        'Confidence': 'High' if row['DIOPT_score'] >= 8 else 'Medium' if row['DIOPT_score'] >= 5 else 'Low'
                    })
        # Add other species pairs as needed
    
    return pd.DataFrame(results)

# Example usage
my_genes = ['white', 'ebony', 'yellow', 'Notch', 'Delta', 'wingless', 'hedgehog']
ortholog_table = batch_orthology_lookup(my_genes)

print("\nBatch Orthology Results:")
print(ortholog_table.to_string(index=False))

ortholog_table.to_csv('batch_orthology_results.csv', index=False)
```

### What You Learned

- How to access Alliance-wide orthology data across all species
- How to map genes between Drosophila and human using FlyBase data
- Understanding DIOPT orthology scores
- Filtering orthologs by confidence level
- Identifying disease-associated orthologs
- Batch orthology conversion workflows

---

## Next Steps and Advanced Usage

### Additional Resources

**Alliance API:**
- https://www.alliancegenome.org/api/swagger-ui
- Programmatic access to all Alliance data
- Real-time queries without downloading files

**JBrowse Genome Browser:**
- https://jbrowse.alliancegenome.org
- Visual exploration of genomic features
- Integration with expression and variant data

**Complete Data Documentation:**
- See [DATA_DOCUMENTATION.md](./DATA_DOCUMENTATION.md) for full file descriptions

### Combining Multiple Data Types

Many powerful analyses combine different data types:

1. **Expression + Disease:** Find genes with tissue-specific expression that are associated with diseases
2. **Interactions + Orthology:** Build conserved interaction networks across species
3. **Variants + Expression:** Correlate genetic variants with expression changes
4. **GO Terms + Expression:** Analyze functional enrichment in expression clusters

### Getting Help

- **Alliance Help:** help@alliancegenome.org
- **GitHub Issues:** https://github.com/alliance-genome/agr_java_software/issues
- **Documentation:** https://www.alliancegenome.org/documentation

---

## Appendix: Common Troubleshooting

### File Download Issues

**Problem:** `wget` command not found
**Solution:** Use `curl` instead:
```bash
curl -O https://s3ftp.flybase.org/releases/current/precomputed_files/genes/fbgn_annotation_ID_current.tsv.gz
```

**Problem:** Download is very slow
**Solution:** For FlyBase data, HTTPS is usually fastest. For Alliance data, use AWS CLI with `--no-sign-request`.

**Problem:** AWS CLI not configured
**Solution:** No configuration needed! Just use `--no-sign-request` flag:
```bash
aws s3 ls s3://mod-datadumps/ --no-sign-request
aws s3 ls s3://s3ftp.flybase.org/releases/ --no-sign-request
```

### Data Parsing Issues

**Problem:** File has extra comment lines
**Solution:** Use `comment='#'` parameter in pandas:
```python
pd.read_csv('file.tsv', sep='\t', comment='#')
```

**Problem:** Memory error loading large files
**Solution:** Use chunking:
```python
chunks = pd.read_csv('large_file.tsv', sep='\t', comment='#', chunksize=10000)
for chunk in chunks:
    # Process each chunk
    process(chunk)
```

### Python Package Issues

**Problem:** Module not found
**Solution:** Install missing packages:
```bash
pip install pandas matplotlib seaborn networkx biopython boto3
```

**Problem:** Version conflicts
**Solution:** Create a virtual environment:
```bash
python3 -m venv alliance_env
source alliance_env/bin/activate  # On Windows: alliance_env\Scripts\activate
pip install -r requirements.txt
```

---

**Document Version:** 2.0
**Last Updated:** 2025-10-17
**Changes:** Updated all S3 bucket references to use correct Alliance (mod-datadumps) and FlyBase (s3ftp.flybase.org) buckets. Removed all credential requirements - all data is publicly accessible with `--no-sign-request`.
**Feedback:** Please report issues or suggestions to help@alliancegenome.org
