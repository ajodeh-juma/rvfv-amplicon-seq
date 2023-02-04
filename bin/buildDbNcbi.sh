#!/bin/sh -e

#-------------------------------------------------------------------------------
#
# Download a GenBank file form NCBI and build the corresponding SnpEff database
#
# Note: It is assumed to be run form SnpEff's directory and that 
#       the 'data' directory is 'data'
#
# Usage example:
#     $ cd ~/snpEff
#     $ ./scripts/buildDbNcbi.sh CP000724.1
#
#															Pablo Cingolani 2015
#-------------------------------------------------------------------------------

#---
# Command line arguments
#---
ID=$1

if [ -z "$ID" ]
then
	echo "Usage: $0 ncbi_genbank_accession"
	exit 1
fi

DIR="data/$ID"
GENE_FILE="$DIR/genes.gbk"

#---
# Download data
#---

# Create db directory
mkdir -p "$DIR" >/dev/null 2>&1 || true

# Download GenBank file
echo "Downloading genome $ID"
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$ID&rettype=gbwithparts&retmode=text" > $GENE_FILE

#---
# Build database
#---

# Add entry to config file
echo "$ID.genome : $ID" >> snpEff.config

# Build database
# java -jar snpEff.jar build -v $ID
snpEff build -v $ID

# snpEff \
# 	NC_014395.1 \
# 	-config /Users/jjuma/PhD_RVF2019/pipelines/rvfvampliconseq/snpEff.config \
# 	/Users/jjuma/PhD_RVF2019/pipelines/rvfvampliconseq/work/02/6d16bf9eb13811b01c027fb186cb0b/08HAB.vcf \
# 	-csvStats 08HAB.snpeff.csv | bgzip -c > 08HAB.snpeff.vcf.gz


# tabix -p vcf -f 08HAB.snpeff.vcf.gz

# bcftools stats 08HAB.snpeff.vcf.gz > 08HAB.bcftools.stats.vcf

# SnpSift extractFields -s "," -e "." 08HAB.snpeff.vcf.gz CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" \
#                 "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
#                 "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
#                 "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" \
#                 "EFF[*].AA_LEN" > 08HAB.snpeff.snpSift.table.txt

