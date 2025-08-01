# Retreiving Annotation via biomart

Annotation from [BioMart](https://useast.ensembl.org/biomart/martview/)  with Ensembl names is most flexible way to retrieve tabular annotation for an organism.

1. The [Biomart](https://useast.ensembl.org/biomart/martview/) start page should look like ...

<img src="annotation_mm_figures/annotation_figures1.png" alt="annotation_figures1" width="80%" style="border:5px solid #ADD8E6;"/>

1. First select the dataset, for gene expression experiment select Ensembl Genes 114 (version 114). The current version as of this workshop.

<img src="annotation_mm_figures/annotation_figures2.png" alt="annotation_figures2" width="80%" style="border:5px solid #ADD8E6;"/>

1. Then the Organism, Here Mouse genes which is based on the GRCm39 genome.

<img src="annotation_mm_figures/annotation_figures3.png" alt="annotation_figures3" width="80%" style="border:5px solid #ADD8E6;"/>

1. You can choose to filter to only a subset of genes. Or a chromosome, or regions. _We won't filter here_. BY default, all genes in the genome are selected.

<img src="annotation_mm_figures/annotation_figures4.png" alt="annotation_figures6" width="80%" style="border:5px solid #ADD8E6;"/>

<img src="annotation_mm_figures/annotation_figures7.png" alt="annotation_figures7" width="80%" style="border:5px solid #ADD8E6;"/>

1. Next select the attributes you want in the table.

<img src="annotation_mm_figures/annotation_figures8.png" alt="annotation_figures8" width="80%" style="border:5px solid #ADD8E6;"/>

1. Expand the 'GENE' tab, and select the attributes you want to retrieve. **HERE** recreate the list you see on the left side.

<img src="annotation_mm_figures/annotation_figures9.png" alt="annotation_figures9" width="80%" style="border:5px solid #ADD8E6;"/>

1. Click "Results" (Top left -ish).

<img src="annotation_mm_figures/annotation_figures10.png" alt="annotation_figures10" width="80%" style="border:5px solid #ADD8E6;"/>

1. Select "GO", to download a tab-separated value (tsv) file.

<img src="annotation_mm_figures/annotation_figures11.png" alt="annotation_figures11" width="80%" style="border:5px solid #ADD8E6;"/>

1. The file will save as "mart_export.txt", put the file into your rnaseq_example directory, rename to "ensembl_mm_114.txt" and open the file in a spreadsheet program to view the annotation.
