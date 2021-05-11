# Blastoid check



This repository contains the code used to check the Blastoids single-cel RNA-seq data

- Xiaodong Liu, Jia Ping Tan, Jan Schröder, Asma Aberkane, John F. Ouyang, Monika Mohenska, Sue Mei Lim, Yu B. Y. Sun, Joseph Chen, Guizhi Sun, Yichen Zhou, Daniel Poppe, Ryan Lister, Amander T. Clark, Owen J. L. Rackham, Jennifer Zenker# & Jose M. Polo&#. **Modelling human blastocysts by reprogramming fibroblasts into iBlastoids.** *Nature* (2021). https://doi.org/10.1038/s41586-021-03372-y.

- Leqian Yu\*, Yulei Wei\*, Jialei Duan\*, Daniel A. Schmitz, Masahiro Sakurai, Lei Wang, Kunhua Wang, Shuhua Zhao, Gary C. Hon# & Jun Wu#. **Blastocyst-like structures generated from human pluripotent stem cells.** *Nature* (2021). https://doi.org/10.1038/s41586-021-03356-y


***


- Combining all datasets: Combine.all.rawData.R
- [Restore cell identity of Zheng et al., 2019: Def_AMLC.R ](html_report/Def_AMLC.html)
- [Restore cell identity of Liu et al., 2021: Def_Iblastoids.R ](html_report/Def_Iblastoids.html)
- [Restore cell identity of Yu et al., 2021: Def_Eblastoids.R ](html_report/Def_Eblastoids.html)
- Quality control for cells and downsampling: QC.D2.downsampling.R
- normalization; Define lineage marker genes; UMAP projection of all datasets; Generating figures : [D2.norm.integration.figure.R](html_report/D2.norm.integration.figure.html)


