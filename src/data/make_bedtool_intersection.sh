#!/bin/sh
#intersect variants with exons using bedtools, needs absolute paths, -a are exons, -b are variants and > XX is the new intersection bed file
#intersect_clinvar_gtex.py provides the data files for this

bedtools intersect -a /mnt/wartmann_pe940/git/disp/data/interim/from_intersect_clinvar_gtex.py/bed_gtex_exons.bed -b /mnt/wartmann_pe940/git/disp/data/interim/from_intersect_clinvar_gtex.py/bed_variant_summary.bed -wa -wb > /mnt/wartmann_pe940/git/disp/data/interim/bed_gtex_exons.bed
         
