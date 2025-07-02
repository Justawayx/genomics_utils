SNPEFF_DIR=/tscc/projects/ps-winzelerlab/TOOLS/snpEff
SNPEFF_SPECIES_ID=Pf3D7v3
VCF=$1

java -jar ${SNPEFF_DIR}/snpEff.jar -formatEFF -o vcf -ud 1000 -c ${SNPEFF_DIR}/snpEff.config ${SNPEFF_SPECIES_ID} $VCF > ${VCF}.ann.txt
