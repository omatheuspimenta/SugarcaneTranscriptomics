#### 5. FEELnc
FEELnc: a tool for long non-coding RNA annotation and its application to the dog transcriptome. From: https://doi.org/10.1093/nar/gkw1306
Install and run following the official repository instructions with conda: https://github.com/tderrien/FEELnc 
```bash
FEELnc_codpot.pl \
    -i lncRNA_candidates.gtf \
    -a SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
    -g SofficinarumxspontaneumR570_771_v2.0.fasta \
    -m intergenic \
    --outdir feelnc_codpot
```  
```bash
FEELnc_classifier.pl \
	-i lncRNA_candidates.gtf.lncRNA.gtf \
	-a reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf > lncRNA_classes.txt
```
