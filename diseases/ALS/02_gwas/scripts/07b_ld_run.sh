#!/bin/bash
# ALS-S2 step 07b (v2): LD expansion of focal ALS GWS variants (r2>=0.8, +/-500kb)
# against 1000G phase 3 GRCh38 (release 20190312), EUR and EAS panels.
# The panel VCFs use bare chromosome names ('9') and '.' IDs, so regions are
# chr-stripped and variant IDs are set to CHROM:POS:REF:ALT; focal variants are
# matched by position.
set -uo pipefail
source /etc/profile.d/modules.sh
module load bcftools/1.23 plink/1.9.0-beta4.4 2>/dev/null || module load bcftools plink

D=/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas
KG=/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data/1000genomes
cd $D/tmp

chroms=$(tail -n +2 ld_focal.tsv | cut -f2 | sort -u)
for c in $chroms; do
  vcf=$KG/ALL.chr${c}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
  sed 's/^chr//' regions_chr${c}.txt > regions_nochr_${c}.txt
  for pop in EUR EAS; do
    if [ "$pop" = "EAS" ]; then
      n=$(awk -F'\t' -v c=$c '$2==c && $4 ~ /EAS/' ld_focal.tsv | wc -l)
      [ "$n" -eq 0 ] && continue
    fi
    out=chr${c}_${pop}
    bcftools view -R regions_nochr_${c}.txt -S samples_${pop}.txt --force-samples $vcf 2> ${out}.bcftools.err \
      | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${out}.vcf.gz
    bcftools index -t -f ${out}.vcf.gz
    awk -F'\t' -v c=$c '$2==c {print $3}' ld_focal.tsv | while read pos; do
      bcftools query -r ${c}:${pos}-${pos} -f '%ID\n' ${out}.vcf.gz
    done | sort -u > focal_ids_${out}.txt || true
    nf=$(wc -l < focal_ids_${out}.txt)
    if [ "$nf" -eq 0 ]; then echo "chr${c} ${pop}: no focal variants in panel"; continue; fi
    plink --vcf ${out}.vcf.gz --keep-allele-order --r2 \
          --ld-snp-list focal_ids_${out}.txt \
          --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 \
          --out ld_${out} > /dev/null 2> ld_${out}.plink.err
    echo "chr${c} ${pop}: focal_in_panel=$nf, ld_pairs=$(tail -n +2 ld_${out}.ld 2>/dev/null | wc -l)"
  done
done
echo DONE
