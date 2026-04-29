#!/bin/bash
#SBATCH --job-name=download_resources_mm39
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/download_resources_mm39_%j.out
#SBATCH --error=logs/download_resources_mm39_%j.err

# =============================================================================
# Download Genomic Resources for Allele-Specific PRO-seq Analysis (mm39/GRCm39)
# =============================================================================
# This script downloads:
# 1. mm39 reference genome (with indices)
# 2. Mouse Genomes Project VCF files (strain-specific variants)
#    - Already in GRCm39 coordinates (Ensembl style: 1,2,3... not chr1,chr2,chr3)
# 3. Chromosome sizes file
#
# For F121-9 cells (129S1/SvImJ × CAST/EiJ F1 hybrid)
# For eventual BL6 × CAST analysis
#
# NOTE: The MGP REL-2112-v8 VCF files are ALREADY in GRCm39 coordinates!
# The VCF uses Ensembl-style chromosome names (1, 2, 3...) not UCSC (chr1, chr2...)
# We handle this by providing both naming conventions.
# =============================================================================

set -euo pipefail

# Configuration
RESOURCE_DIR="${1:-/n/scratch/users/a/alb1273/pausing_phase_project_intermediates/resources}"
THREADS=4

# Create directory structure
mkdir -p "${RESOURCE_DIR}"/{genome,vcf,indices}
mkdir -p logs

cd "${RESOURCE_DIR}"

echo "=============================================="
echo "Starting resource download (mm39/GRCm39): $(date)"
echo "Resource directory: ${RESOURCE_DIR}"
echo "=============================================="

# =============================================================================
# 1. Download mm39 Reference Genome
# =============================================================================
echo ""
echo "[1/5] Downloading mm39 reference genome..."

GENOME_DIR="${RESOURCE_DIR}/genome"
cd "${GENOME_DIR}"

# We need BOTH UCSC (chr1, chr2...) and Ensembl (1, 2...) versions
# because the VCF uses Ensembl naming but many tools expect UCSC naming

# UCSC version (chr1, chr2, chr3...)
if [[ ! -f "mm39.fa" ]]; then
    echo "  Downloading mm39.fa.gz from UCSC (chr1, chr2... naming)..."
    wget -q --show-progress \
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" \
        -O mm39.fa.gz
    
    echo "  Decompressing..."
    gunzip mm39.fa.gz
    
    echo "  Creating samtools index..."
    samtools faidx mm39.fa
else
    echo "  mm39.fa already exists, skipping download"
fi

# Ensembl version (1, 2, 3...) - needed for VCF compatibility
if [[ ! -f "mm39.ensembl.fa" ]]; then
    echo "  Downloading Ensembl GRCm39 (1, 2... naming for VCF compatibility)..."
    wget -q --show-progress \
        "https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
        -O mm39.ensembl.fa.gz
    
    echo "  Decompressing..."
    gunzip mm39.ensembl.fa.gz
    
    # Rename to standard name
    if [[ -f "Mus_musculus.GRCm39.dna.primary_assembly.fa" ]]; then
        mv Mus_musculus.GRCm39.dna.primary_assembly.fa mm39.ensembl.fa
    fi
    
    echo "  Creating samtools index..."
    samtools faidx mm39.ensembl.fa
else
    echo "  mm39.ensembl.fa already exists, skipping download"
fi

# Download chromosome sizes (UCSC style)
if [[ ! -f "mm39.chrom.sizes" ]]; then
    echo "  Downloading chromosome sizes..."
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes" \
        -O mm39.chrom.sizes
    
    # Create a version with only standard chromosomes (no random, Un, etc.)
    grep -E "^chr[0-9]+|^chr[XYM]" mm39.chrom.sizes | sort -V > mm39.chrom.sizes.standard
    
    # Create Ensembl-style version (for VCF work)
    sed 's/^chr//' mm39.chrom.sizes.standard | sed 's/^M$/MT/' > mm39.chrom.sizes.ensembl
fi

echo "  Reference genome ready: ${GENOME_DIR}/mm39.fa (UCSC) and mm39.ensembl.fa (Ensembl)"

# =============================================================================
# 2. Download Mouse Genomes Project VCF Files
# =============================================================================
echo ""
echo "[2/5] Downloading Mouse Genomes Project VCF files..."
echo "  NOTE: These are already in GRCm39 coordinates with Ensembl-style chr names"

VCF_DIR="${RESOURCE_DIR}/vcf"
cd "${VCF_DIR}"

# The Mouse Genomes Project (Sanger) REL-2112-v8 is in GRCm39 coordinates
# VCF uses Ensembl chromosome naming (1, 2, 3...) not UCSC (chr1, chr2, chr3...)
# 
# Key strains for our analysis:
# - 129S1_SvImJ: One parent of F121-9
# - CAST_EiJ: Other parent of F121-9 (also used in BL6 x CAST)
# - C57BL_6NJ: Close to reference (for BL6 x CAST analysis)

MGP_FTP="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels"

if [[ ! -f "mgp_REL2021_snps.vcf.gz" ]]; then
    echo "  Downloading Mouse Genomes Project SNP VCF (REL-2112-v8, GRCm39)..."
    echo "  This is ~22GB and may take a while..."
    
    wget -q --show-progress \
        "${MGP_FTP}/mgp_REL2021_snps.vcf.gz" \
        -O mgp_REL2021_snps.vcf.gz
    
    # Download the CSI index (not TBI - file is too large for standard tabix)
    wget -q "${MGP_FTP}/mgp_REL2021_snps.vcf.gz.csi" \
        -O mgp_REL2021_snps.vcf.gz.csi
else
    echo "  MGP SNP VCF already exists, skipping download"
fi

# Verify this is indeed GRCm39
echo "  Verifying VCF genome build..."
zcat mgp_REL2021_snps.vcf.gz 2>/dev/null | head -100 | grep -E "##reference|##contig=<ID=1," | head -2 || \
    echo "  (Could not verify - check manually)"

# =============================================================================
# 3. Extract Strain-Specific VCFs
# =============================================================================
echo ""
echo "[3/5] Extracting strain-specific variant files..."

# Load bcftools if in module system
module load bcftools 2>/dev/null || true

# List available strains
echo "  Checking sample names in VCF..."
bcftools query -l mgp_REL2021_snps.vcf.gz | head -20

# --- F121-9: 129S1_SvImJ vs CAST_EiJ ---
echo ""
echo "  Extracting 129S1_SvImJ and CAST_EiJ variants for F121-9..."

if [[ ! -f "F121-9_het_snps.vcf.gz" ]]; then
    # Extract only sites where 129S1 and CAST differ (heterozygous in F1)
    # Filter for sites that are:
    # 1. Biallelic SNPs (not indels)
    # 2. Heterozygous in F1 (different genotypes between parents)
    
    bcftools view \
        -s 129S1_SvImJ,CAST_EiJ \
        -v snps \
        -m2 -M2 \
        mgp_REL2021_snps.vcf.gz \
    | bcftools view \
        -i 'GT[0]!=GT[1]' \
        -Oz -o F121-9_het_snps.vcf.gz
    
    tabix -p vcf F121-9_het_snps.vcf.gz
    
    # Count variants
    N_VARS=$(bcftools view -H F121-9_het_snps.vcf.gz | wc -l)
    echo "  Found ${N_VARS} heterozygous SNPs between 129S1 and CAST"
fi

# --- BL6 x CAST ---
echo "  Extracting C57BL_6NJ and CAST_EiJ variants for BL6xCAST..."

if [[ ! -f "BL6xCAST_het_snps.vcf.gz" ]]; then
    bcftools view \
        -s C57BL_6NJ,CAST_EiJ \
        -v snps \
        -m2 -M2 \
        mgp_REL2021_snps.vcf.gz \
    | bcftools view \
        -i 'GT[0]!=GT[1]' \
        -Oz -o BL6xCAST_het_snps.vcf.gz
    
    tabix -p vcf BL6xCAST_het_snps.vcf.gz
    
    N_VARS=$(bcftools view -H BL6xCAST_het_snps.vcf.gz | wc -l)
    echo "  Found ${N_VARS} heterozygous SNPs between BL6 and CAST"
fi

# =============================================================================
# 4. Create UCSC-style VCFs (add "chr" prefix)
# =============================================================================
echo ""
echo "[4/5] Creating UCSC-style VCFs (adding 'chr' prefix to chromosome names)..."

# Many downstream tools expect UCSC-style chromosome names
# Create versions with "chr" prefix

for VCF in F121-9_het_snps.vcf.gz BL6xCAST_het_snps.vcf.gz; do
    BASE=$(basename ${VCF} .vcf.gz)
    UCSC_VCF="${BASE}.ucsc.vcf.gz"
    
    if [[ ! -f "${UCSC_VCF}" ]]; then
        echo "  Converting ${BASE} to UCSC style..."
        
        # Add "chr" prefix to chromosome names
        # Use bcftools annotate with a rename file for clean conversion
        
        # Create chromosome rename map
        echo "1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
X chrX
Y chrY
MT chrM" > chr_rename.txt
        
        bcftools annotate \
            --rename-chrs chr_rename.txt \
            -Oz -o ${UCSC_VCF} \
            ${VCF}
        
        tabix -p vcf ${UCSC_VCF}
        echo "  Created ${UCSC_VCF}"
    else
        echo "  ${UCSC_VCF} already exists, skipping"
    fi
done

# Clean up
rm -f chr_rename.txt

# =============================================================================
# 5. Generate Summary Statistics
# =============================================================================
echo ""
echo "[5/5] Generating summary statistics..."

# Create a summary of SNP density per chromosome
echo "  Calculating SNP density per chromosome..."

for VCF in F121-9_het_snps.vcf.gz BL6xCAST_het_snps.vcf.gz; do
    BASE=$(basename ${VCF} .vcf.gz)
    
    echo "  Processing ${BASE}..."
    
    # Count SNPs per chromosome
    bcftools view -H ${VCF} \
        | cut -f1 \
        | sort -V \
        | uniq -c \
        | awk '{print $2"\t"$1}' \
        > ${BASE}_snp_counts.txt
    
    # Calculate density (SNPs per kb) using Ensembl chrom sizes
    echo -e "chr\tsnps\tlength\tsnps_per_kb" > ${BASE}_snp_density.txt
    while read chr count; do
        len=$(grep -w "^${chr}" ${RESOURCE_DIR}/genome/mm39.chrom.sizes.ensembl | cut -f2)
        if [[ -n "$len" ]]; then
            density=$(echo "scale=4; ${count} / (${len}/1000)" | bc)
            echo -e "${chr}\t${count}\t${len}\t${density}"
        fi
    done < ${BASE}_snp_counts.txt >> ${BASE}_snp_density.txt
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Resource download complete: $(date)"
echo "=============================================="
echo ""
echo "IMPORTANT NOTES:"
echo "  - VCF files use Ensembl chromosome naming (1, 2, 3... not chr1, chr2, chr3)"
echo "  - UCSC-style VCFs (.ucsc.vcf.gz) are also provided with 'chr' prefix"
echo "  - Use mm39.ensembl.fa reference with original VCFs"
echo "  - Use mm39.fa reference with .ucsc.vcf.gz files"
echo ""
echo "Files created:"
echo "  Reference genomes:"
echo "    ${GENOME_DIR}/mm39.fa          (UCSC naming: chr1, chr2...)"
echo "    ${GENOME_DIR}/mm39.ensembl.fa  (Ensembl naming: 1, 2...)"
echo ""
echo "  Chromosome sizes:"
echo "    ${GENOME_DIR}/mm39.chrom.sizes          (UCSC)"
echo "    ${GENOME_DIR}/mm39.chrom.sizes.ensembl  (Ensembl)"
echo ""
echo "  VCF files (Ensembl chromosome naming):"
echo "    ${VCF_DIR}/mgp_REL2021_snps.vcf.gz      (Full MGP, GRCm39)"
echo "    ${VCF_DIR}/F121-9_het_snps.vcf.gz       (129S1 vs CAST hets)"
echo "    ${VCF_DIR}/BL6xCAST_het_snps.vcf.gz     (BL6 vs CAST hets)"
echo ""
echo "  VCF files (UCSC chromosome naming):"
echo "    ${VCF_DIR}/F121-9_het_snps.ucsc.vcf.gz"
echo "    ${VCF_DIR}/BL6xCAST_het_snps.ucsc.vcf.gz"
echo ""
echo "SNP Summary (F121-9):"
head -6 ${VCF_DIR}/F121-9_het_snps_snp_density.txt
echo "..."
echo ""
echo "Total heterozygous SNPs:"
echo "  F121-9:   $(bcftools view -H ${VCF_DIR}/F121-9_het_snps.vcf.gz | wc -l)"
echo "  BL6xCAST: $(bcftools view -H ${VCF_DIR}/BL6xCAST_het_snps.vcf.gz | wc -l)"
