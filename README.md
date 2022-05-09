# Hi-C data

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

A (continuously updated) collection of references to Hi-C data and papers. Predominantly human/mouse Hi-C data, with replicates. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Large collections](#large-collections)
  - [Lieberman-Aiden lab](#lieberman-aiden-lab)
  - [Leonid Mirny lab](#leonid-mirny-lab)
  - [Bing Ren lab](#bing-ren-lab)
  - [Feng Yue lab](#feng-yue-lab)
  - [4D Nucleome](#4d-nucleome)
- [Cancer](#cancer)
  - [BRCA](#brca)
- [Tissue-specific](#tissue-specific)
  - [ENCODE](#encode)
  - [Brain](#brain)
  - [Cell lines](#cell-lines)
  - [Non-human data](#non-human-data)
- [Differential Hi-C](#differential-hi-c)
- [Timecourse Hi-C](#timecourse-hi-c)
- [Promoter-capture Hi-C](#promoter-capture-hi-c)
- [Integrative Hi-C](#integrative-hi-c)
- [Single-cell Hi-C](#single-cell-hi-c)
- [Micro-C](#micro-c)
- [GAM](#gam)
- [Imaging](#imaging)
- [CTCF](#ctcf)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Large collections

- [3DIV](http://kobic.kr/3div/download) - database of uniformly processed 315 Hi-C datasets, 80 human cell/tissue types. Bait-centric (SNP rsID, gene name, hg19 coordinates) visualization of long-range interactions in context of epigenomic (histone, enhancers) signals, numerical results. Custom BWA-MEM pipeline, Bias, distance effect removed. Coordinates of significant interactions, with annotations, are available for (FTP) download, http://kobic.kr/3div/download
    - Yang, Dongchan, Insu Jang, Jinhyuk Choi, Min-Seo Kim, Andrew J Lee, Hyunwoong Kim, Junghyun Eom, Dongsup Kim, Inkyung Jung, and Byungwook Lee. “[3DIV: A 3D-Genome Interaction Viewer and Database](https://doi.org/10.1093/nar/gkx1017).” Nucleic Acids Research 46, no. D1 (January 4, 2018)

- [Chorogenome](http://chorogenome.ie-freiburg.mpg.de:5003/) resource: Processed data (Hi-C, ChIP-seq) for Drosophila, Mouse, Human, http://chorogenome.ie-freiburg.mpg.de:5003/
    - Ramírez, Fidel, Vivek Bhardwaj, Laura Arrigoni, Kin Chung Lam, Björn A. Grüning, José Villaveces, Bianca Habermann, Asifa Akhtar, and Thomas Manke. “[High-Resolution TADs Reveal DNA Sequences Underlying Genome Organization in Flies](https://doi.org/10.1038/s41467-017-02525-w).” Nature Communications 9, no. 1 (December 2018). 

- [GITAR: An Open Source Tool for Analysis and Visualization of Hi-C Data](https://www.genomegitar.org/processed-data.html) - Includes a large collection of standardized processed data from 4D Nucleome. 20 hg38 and 2 mm10 datasets normalized by Yaffe-Tanay method, downloadable, include directionality index, HMM states, TAD analysis results. Text and HDF5 formats. https://www.genomegitar.org/processed-data.html

- [4DGenome](https://4dgenome.research.chop.edu/) - 3D significant interactions, from different literature sources
    - Teng, Li, Bing He, Jiahui Wang, and Kai Tan. “[4DGenome: A Comprehensive Database of Chromatin Interactions](https://doi.org/10.1093/bioinformatics/btv158).” Bioinformatics (Oxford, England) 31, no. 15 (August 1, 2015)


## Lieberman-Aiden lab

All HiC data released by Lieberman-Aiden group. Links to Amazon storage and GEO studies. http://aidenlab.org/data.html

- Vian, Laura, Aleksandra Pękowska, Suhas S.P. Rao, Kyong-Rim Kieffer-Kwon, Seolkyoung Jung, Laura Baranello, Su-Chen Huang, et al. “[The Energetics and Physiological Impact of Cohesin Extrusion](https://doi.org/10.1016/j.cell.2018.03.072).” Cell 173, no. 5 (May 2018) - Architectural stripes, created by extensive loading of cohesin near CTCF anchors, with Nipbl and Rad21 help. Little overlap between B cells and ESCs. Architectural stripes are sites for tumor-inducing TOP2beta DNA breaks. ATP is required for loop extrusion, cohesin translocation, but not required for maintenance, Replication of transcription is not important for loop extrusion. Zebra algorithm for detecting architectural stripes, image analysis, math in Methods. Human lymphoblastoid cells, mouse ESCs, mouse B-cells activated with LPS, CH12 B lymphoma cells, wild-type, treated with hydroxyurea (blocks DNA replication), flavopiridol (blocks transcription, PolII elongation), oligomycin (blocks ATP). Many other data types (e.g., ChIP-seq, ATAC-seq) [GSE82144](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82144), [GSE98119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98119)

- Lieberman-Aiden, Erez, Nynke L. van Berkum, Louise Williams, Maxim Imakaev, Tobias Ragoczy, Agnes Telling, Ido Amit, et al. “[Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome](https://doi.org/10.1126/science.1181369).” Science (New York, N.Y.) 326, no. 5950 (October 9, 2009) Gm12878, K562 cells. HindIII, NcoI enzymes. Two-three replicates. [GSE18199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199)

- Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. “[A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping](https://doi.org/10.1016/j.cell.2014.11.021).” Cell 159, no. 7 (December 18, 2014) - Human Gm12878, K562, IMR90, NHEC, HeLa cells, Mouse CH12 cells. Different digestion enzymes (HindIII, NcoI, Mbol, DpnII), different dilutions. Up to 35 biological replicates for Gm12878. [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), [Supplementary Table S1. Hi-C meta-data](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/1bcd7dea-7dbe-45af-8664-4ddd7fb54bc6/mmc2.xlsx)

- Sanborn, Adrian L., Suhas S. P. Rao, Su-Chen Huang, Neva C. Durand, Miriam H. Huntley, Andrew I. Jewett, Ivan D. Bochkov, et al. “[Chromatin Extrusion Explains Key Features of Loop and Domain Formation in Wild-Type and Engineered Genomes](https://doi.org/10.1073/pnas.1518552112).” Proceedings of the National Academy of Sciences of the United States of America 112, no. 47 (November 24, 2015).  HAP1, derived from chronic myelogenous leukemia cell line. Replicates. [GSE74072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74072)
    - [Sanborn_Aiden_2015_st01.xlsx](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4664323/bin/pnas.1518552112.st01.xlsx) - Sheet "S6 - Hi-C experiments" contains information about Hi-C experiments.

- Rao, Suhas S.P., Su-Chen Huang, Brian Glenn St Hilaire, Jesse M. Engreitz, Elizabeth M. Perez, Kyong-Rim Kieffer-Kwon, Adrian L. Sanborn, et al. “[Cohesin Loss Eliminates All Loop Domains](https://doi.org/10.1016/j.cell.2017.09.026).” Cell 171, no. 2 (2017) - HCT-116 human colorectal carcinoma cells. Timecourse, replicates under different conditions. [GSE104334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104334)
    - [Rao_Aiden_2017_mmc1.xlsx](http://www.cell.com/cms/attachment/2112693605/2084046365/mmc1.xlsx) - Table S1. Summary Statistics for the Datasets.


## Leonid Mirny lab

http://mirnylab.mit.edu/

- Data from multiple studies, in one place, in `.cool` format: ftp://cooler.csail.mit.edu/coolers
- Convert to any other format with `cooler` https://cooler.readthedocs.io/

- Depletion of the cohesin-loading factor Nipbl. Three conditions: wild-type, tamoxifen control and deltaNipbl mice liver. TADs disappear, A/B compartments reinforced, minimal nonspecific effect on gene expression. Disappearing TADs unmask finer level of chromatin organization that is better associated with epigenetic landscape. TADs and compartments are independent types of chromosomal organization, but overlapping. Ideas: Excluding low-coverage bins using MAD-max procedure (Methods). Global compartmentalization. [Lavaburst](https://github.com/nvictus/lavaburst) - TAD detection using Filippova method. [Tethered Hi-C, H3K4me3, H3K27ac, CTCF, Rad21, Smc3 ChIP-seq and RNA-seq data and visualization](http://mirnylab.mit.edu/projects/nipbl/), [GEO GSE93431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93431)
    - Schwarzer, Wibke, Nezar Abdennur, Anton Goloborodko, Aleksandra Pekowska, Geoffrey Fudenberg, Yann Loe-Mie, Nuno A. Fonseca, et al. “[Two Independent Modes of Chromatin Organization Revealed by Cohesin Removal](https://doi.org/10.1038/nature24281).” Nature, (02 2017)



## Bing Ren lab

http://chromosome.sdsc.edu/mouse/hi-c/download.html

Raw and normalized chromatin interaction matrices and TADs defined with DomainCaller. Mouse ES, cortex, Human ES, IMR90 fibroblasts. Two replicates per condition. GEO accession: [GSE35156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35156), [GSE43070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43070)

- 3D variability between 20 humans, lymphoblastoid cell lines, associated with variation in gene expression, histone modifications, transcription factor binding. Genetic variation (SNPs) is associated with loop strength, contact insulation, directionality, density of local contacts, SNPs in CTCF binding sites - QTLs for these. WASP approach to address allelic mapping biases, HiCNorm normalization to remove GC, mappability, fragment length biases, BNBC quantile normalization across samples. 40kb data, detecting A/B compartments (PC1), directionality index (DI), insulation score (INS), frequently interacting regions (FIRE score). Variability detected using limma:eBayes function. IWH for multiple testing correction. Power calculation for QTL detection in Hi-C data. Data and code: Hi-C BAM files, matrices, full QTL results, 3D variable regions, SNPs at http://renlab.sdsc.edu/renlab_website/download/iqtl/, http://renlab.sdsc.edu/iQTL/
    - Gorkin, David U., Yunjiang Qiu, Ming Hu, Kipper Fletez-Brant, Tristin Liu, Anthony D. Schmitt, Amina Noor, et al. “[Common DNA Sequence Variation Influences 3-Dimensional Conformation of the Human Genome](https://doi.org/10.1101/592741).” Preprint. Genomics, March 30, 2019. 

- Normal human cells, brain (dorsolateral prefrontal cortex, hippocampus), adrenal, bladder, lung, ovary, pancreas, etc. 21 human cell lines and primary tissues. Some replicates. [GSE87112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112). Used in [HiCDB paper](https://doi.org/10.1093/nar/gky789)
    - Schmitt, Anthony D., Ming Hu, Inkyung Jung, Zheng Xu, Yunjiang Qiu, Catherine L. Tan, Yun Li, et al. “[A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome](https://doi.org/10.1016/j.celrep.2016.10.061).” Cell Reports 17, no. 8 (November 2016) 

- Dixon, Jesse R., Siddarth Selvaraj, Feng Yue, Audrey Kim, Yan Li, Yin Shen, Ming Hu, Jun S. Liu, and Bing Ren. “[Topological Domains in Mammalian Genomes Identified by Analysis of Chromatin Interactions](https://doi.org/10.1038/nature11082).” Nature 485, no. 7398 (April 11, 2012)

- Jin, Fulai, Yan Li, Jesse R. Dixon, Siddarth Selvaraj, Zhen Ye, Ah Young Lee, Chia-An Yen, Anthony D. Schmitt, Celso A. Espinoza, and Bing Ren. “[A High-Resolution Map of the Three-Dimensional Chromatin Interactome in Human Cells](https://doi.org/10.1038/nature12644).” Nature 503, no. 7475 (November 14, 2013)


## Feng Yue lab

- [Data from Yue lab](http://3dgenome.fsm.northwestern.edu/publications.html)

- [NeoLoopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder), Wang, Xiaotao, Jie Xu, Baozhen Zhang, Ye Hou, Fan Song, Huijue Lyu, and Feng Yue. “[Genome-Wide Detection of Enhancer-Hijacking Events from Chromatin Interaction Data in Rearranged Genomes](https://doi.org/10.1038/s41592-021-01164-w).” Nature Methods, (June 2021)
    - [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01164-w/MediaObjects/41592_2021_1164_MOESM3_ESM.xlsx) - Details of the 50 cancer Hi-C datasets, references to GEO and 4DNucleome.
    - [Supplementary Table 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01164-w/MediaObjects/41592_2021_1164_MOESM4_ESM.xlsx) - Coordinates (hg38) of large SVs detected in each sample.
    - [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01164-w/MediaObjects/41592_2021_1164_MOESM5_ESM.xlsx) - Genomic coordinates of the detected neoloops in each sample.
    - [Supplementary Table 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01164-w/MediaObjects/41592_2021_1164_MOESM6_ESM.xlsx) - List of neoloop-involved genes identified in each sample.
    - [Supplementary Table 8](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01164-w/MediaObjects/41592_2021_1164_MOESM8_ESM.xlsx) - List of the annotated enhancer-hijacking events in 11 cancer cell lines: A549, K562, LNCaP, MCF7, T47D, HepG2, SK-MEL-5, NCI-H460, PANC-1, HT-1080 and C4-2B.

- [3D Genome Browser](http://3dgenome.fsm.northwestern.edu) - Classical datasets for TAD/loop identification, provided as raw and normalized matrices, genomic coordinates of TADs/loops, tools for various 3C data analysis.

- Iyyanki, Tejaswi. “[Subtype-Associated Epigenomic Landscape and 3D Genome Structure in Bladder Cancer](https://doi.org/10.1186/s13059-021-02325-y),” Genome Biology, 15 April 2021 - 3D genomics of bladder cancer. 4 cancer cell lines (luminal: RT4 and SW780; basal: SCABER and HT1376), 5 patients. H3K27ac ChIP-seq, RNA-seq (DESeq2), ATAC-seq (TCGA), Hi-C data (Arima, hg19). Peakachu for loop prediction, CNVs with HiNT and Hi-Cbreakfinder.
    - [Additional file 14: Table S13.](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-021-02325-y/MediaObjects/13059_2021_2325_MOESM14_ESM.xlsx) - Luminal-specific and Basal-specific loops.
    - [GEO GSE148079](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148079)
    - [GitHub](https://github.com/Qixuan771/Source-code-for-bladder-cancer-project)


## 4D Nucleome 

- [4D Nucleome Data Portal](https://data.4dnucleome.org/) - chromatin conformation-related datasets, uniformly prodessed, integrated with the HiGlass genome browser, comparison possible. Overview of first and second phases of the 4DN project. Other repositories that host Hi-C and similar datasets include the ENCODE portal, NCBI's GEO and EMBL-EBI’s ArrayExpress. 
    - Reiff, Sarah B, Andrew J Schroeder, Koray Kirli, Andrea Cosolo, Clara Bakker, Luisa Mercado, Soohyun B Lee, et al. “[The 4D Nucleome Data Portal: A Resource for Searching and Visualizing Curated Nucleomics Data](https://doi.org/10.1101/2021.10.14.464435).” Preprint. Genomics, October 15, 2021. 
    - Genomics datasets - [Hi-C](https://data.4dnucleome.org/resources/data-collections/high-resolution-hic-datasets) (in situ, dilution), Micro-C, DNase Hi-C Hi-C 3.0, Capture Hi-C, TCC, single-cell variants, SPRITE, GAM. ChIA-PET, ChIA-Drop, PLAC-seq, ChIP-seq, CUT&RUN, Repli-seq, others. .cool and .mcool formats, A/B compartments and TAD boundaries (insulation score) are called using cooltools. Data are in [three tiers](https://4dnucleome.org/cell-lines.html): Tier 1 (H1-ESC, GM12878, IMR90, HFF-hTERT (clone 6), and WTC-11), Tier 2 and untiered. 
    - [Microscopy datasets](https://data.4dnucleome.org/microscopy-data-overview) - standard FISH (DNA or RNA), multi-loci FISH, high-throughput FISH, dynamic single particle tracking. 
    - [Hi-C Processing Pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline) (BWA MEM with the -SP5M option)
    - [Domain Calling Pipelines](https://data.4dnucleome.org/resources/data-analysis/insulation_compartment_scores)
    - [4DN Visualization Workspace](https://data.4dnucleome.org/tools/visualization)
    - [4DN Software](https://www.4dnucleome.org/software.html) - alphabetical list of Hi-C software.

- [4dnucleome.eu](https://www.4dnucleome.eu/) - 4DNucleome Initiative in Europe

- [MuGVRE](https://www.multiscalegenomics.eu/MuGVRE/) - The MuG Virtual Research Environment supports the expanding 3D/4D genomics community by developing tools to integrate and visualize genomics data from sequence to 3D/4D chromatin dynamics data

# Cancer

- Changes in 3D genome are associated with CNVs in multiple myeloma cells (RPMI-8226 trt- and tetraploid, U266 nearly diploid). The number of TADs increases by \~25%, they become smaller, \~20% switch compartment. ICE normalization better accounts for CNVs than HiCNorm. CNV breakpoints overlap with TAD boundaries. 40kb resolution, replicates. [Code](https://github.com/ChengLiLab/myeloma), Hi-C, WGS, RNA-seq data [GSE87585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87585)
    - Wu, Pengze, Tingting Li, Ruifeng Li, Lumeng Jia, Ping Zhu, Yifang Liu, Qing Chen, Daiwei Tang, Yuezhou Yu, and Cheng Li. “[3D Genome of Multiple Myeloma Reveals Spatial Genome Disorganization Associated with Copy Number Variations](https://doi.org/10.1038/s41467-017-01793-w).” Nature Communications 8, no. 1 (December 2017)

- Curtaxins drugs affect 3D genome by DNA intercalation but without inducing DNA damage, compromise enhancer-promoter interactions, suppress oncogene expression, including MYC family genes, downregulates survival genes, partially disrupt TAD borders, decreases short-range interactions, the level of spatial segregation of the A/B compartments, depletes CTCF but not other factors. Hi-C in HT1080 fibrosarcoma cells. Data: Hi-C and CTCF ChIP-seq in duplicates [GSE122463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122463), gene expression in MM1.S and HeLa S3 cells [GSE117611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117611), H3K27ac [GSE117409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117409), nascent RNA transcription [GSE107633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107633)
    - Kantidze, Omar L., Artem V. Luzhin, Ekaterina V. Nizovtseva, Alfiya Safina, Maria E. Valieva, Arkadiy K. Golov, Artem K. Velichko, et al. “[The Anti-Cancer Drugs Curaxins Target Spatial Genome Organization](https://doi.org/10.1038/s41467-019-09500-7).” Nature Communications 10, no. 1 (December 2019). 

- 3D genomics of glioblastoma. Replicate samples from three patients. Sub-5kb-resolution Hi-C data, integration with ChIP- and RNA-seq. Data: Six Hi-C replicates, [EGAS00001003493](https://ega-archive.org/studies/EGAS00001003493), ChIP-seq [GSE121601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121601), RNA-seq data [EGAS00001003700](https://www.ebi.ac.uk/ega/studies/EGAS00001003700). [Processed data](https://wangftp.wustl.edu/hubs/johnston_gallo/)
    - Johnston, Michael J., Ana Nikolic, Nicoletta Ninkovic, Paul Guilhamon, Florence M.G. Cavalli, Steven Seaman, Franz J. Zemp, et al. “[High-Resolution Structural Genomics Reveals New Therapeutic Vulnerabilities in Glioblastoma](https://doi.org/10.1101/gr.246520.118).” Genome Research 29, no. 8 (August 2019) 

- Ten non-replicated Hi-C datasets. Two human lymphoblastoid cell lines with known chromosomal translocations (FY1199 and DD1618), transformed mouse cell line (EKLF), six human brain tumours: five glioblastomas ( GB176, GB180, GB182, GB183 and GB238) and one anaplastic astrocytoma (AA86), a normal human cell line control (GM07017). [GSE81879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81879)
- Harewood, Louise, Kamal Kishore, Matthew D. Eldridge, Steven Wingett, Danita Pearson, Stefan Schoenfelder, V. Peter Collins, and Peter Fraser. “[Hi-C as a Tool for Precise Detection and Characterisation of Chromosomal Rearrangements and Copy Number Variation in Human Tumours](https://doi.org/10.1186/s13059-017-1253-8).” Genome Biology 18, no. 1 (December 2017). 

- Prostate cancer, normal. RWPE1 prostate epithelial cells transfected with GFP or ERG oncogene. Two biological and up to four technical replicates. [GSE37752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37752)
    - Rickman, David S., T. David Soong, Benjamin Moss, Juan Miguel Mosquera, Jan Dlabal, Stéphane Terry, Theresa Y. MacDonald, et al. “[Oncogene-Mediated Alterations in Chromatin Conformation](https://doi.org/10.1073/pnas.1112570109).” Proceedings of the National Academy of Sciences of the United States of America 109, no. 23 (June 5, 2012)

- Taberlay, Phillippa C., Joanna Achinger-Kawecka, Aaron T. L. Lun, Fabian A. Buske, Kenneth Sabir, Cathryn M. Gould, Elena Zotenko, et al. “[Three-Dimensional Disorganization of the Cancer Genome Occurs Coincident with Long-Range Genetic and Epigenetic Alterations](https://doi.org/10.1101/gr.201517.115).” Genome Research 26, no. 6 (June 2016)
- Cancer, normal Hi-C. Prostate epithelial cells, PC3, LNCaP. Two-three replicates. [GSE73785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73785)

- Haplotype-resolved Hi-C of GM12878, integrated with RNA-seq and Bru-seq (nascent mRNA). Investigation of Monoallelic expression (MAE) and Allele-Biased expression (ABE). [GEO GSE159813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159813)
    - Lindsly, Stephen, Wenlong Jia, Haiming Chen, Sijia Liu, Scott Ronquist, Can Chen, Xingzhao Wen, et al. “[Functional Organization of the Maternal and Paternal Human 4D Nucleome](https://doi.org/10.1101/2020.03.15.992164),” bioRxiv, June 17, 2021.

## BRCA

- Comparative characterization of 3D genomics in TNBC. Cell lines (HMEC as normal and 5 BRCA subtypes, by the order of aggressiveness: T47D, ZR7530, HCC1954, HCC70, BT549). TNBC shows most dramatic changes, partially conserved across TNBC cell lines and TNBC tissues. TADs (CaTCH), loops (HiCCUPS), compartment (PC1) analyses. Local interactions are lost, "normal" TAD interactions weakened but TNBC TADs strenghtened; those changes are associated with CTCF loss/gain. 3В changes are associated with gene expression changes. Hi-C (replicates), ChIP-seq (CTCF, H3K27ac), RNA-seq, and ATAC-seq data are at [GSE167154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167150). <details>
    <summary>Paper</summary>
    Kim, Taemook, Sungwook Han, Yujin Chun, Hyeokjun Yang, Hyesung Min, Sook Young Jeon, Jang-il Kim, Hyeong-Gon Moon, and Daeyoup Lee. “Comparative Characterization of 3D Chromatin Organization in Triple-Negative Breast Cancers.” Experimental & Molecular Medicine, May 5, 2022. https://doi.org/10.1038/s12276-022-00768-2.
</details>

- BRCA gene targets regulated by SNPs - Capture-C of chromatin interactions centered on causal variants and promoters of causal genes (Variant- and Promoter Capture Hi-C) in six human mammary epithelial (B80T5, MCF10A) and breast cancer (MCF7, T47D, MDAMB231, Hs578T) cell lines. HindIII fragments, CHiCAGO and Peaky for significant interaction calling. PCA on interactions separates cell types, significant interactions enriched in epigenomic elements. 651 target genes at 139 independent breast cancer risk signals. Table 1 - top priority target genes. [HiCUP-processed capture Hi-C data (hg19)](https://osf.io/2cnw7/), [code](https://github.com/jmbeesley/Beesley_GenomeBiol2019), [Supplementary tables](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1877-y#Sec35), Tables S11 - 651 target genes. <details>
    <summary>Paper</summary>
    Beesley, Jonathan, Haran Sivakumaran, Mahdi Moradi Marjaneh, Luize G. Lima, Kristine M. Hillman, Susanne Kaufmann, Natasha Tuano, et al. “Chromatin Interactome Mapping at 139 Independent Breast Cancer Risk Signals.” Genome Biology 21, no. 1 (December 2020) https://doi.org/10.1186/s13059-019-1877-y
</details>

- Capture Hi-C (CHi-C) to annotate 63 breast cancer risk loci. 110 target genes at 33 loci, supported bu other evidence (eQTLs, disease-specific survival). Two ER+ breast cancer cell lines (T-47D, ZR-75-1), two ER− breast cancer cell lines (BT-20, MDA- MB-231), one “normal” breast epithelial cell line (Bre80-Q-TERT (Bre80)) and a non-breast lymphoblastoid cell line (GM06990). Approx 40% of interaction peaks are present in multiple cell lines. More interactions within TADs. [WashU session with all CHi-C interaction peaks](https://bit.ly/CHiC-BC). **Table 2** Risk loci which formed interaction peaks directly (N = 33) or via an adjacent risk locus (N = 3) with 110 target genes (locus, SNP, gene targets, nearest gene). **Table 3** Nine CHi-C putative target genes that were statistically significant eQTLs (FDR adjusted P < 0.1) (locus, SNP, gene, p-values in all, ER+/- cancers). **Table 4** Six CHi-C putative target genes for which there was orthogonal support for at least two additional data sources. [PRJEB23968](https://www.ebi.ac.uk/ena/browser/view/PRJEB23968?show=reads) - FASTQ files. <details>
    <summary>Supplementary material</summary>
    https://www.nature.com/articles/s41467-018-03411-9#Sec23
    - Supplementary Data 1: Captured genomic regions (Locus, SNP, hg19 coordinates, size, reference)
    - Supplementary Data 2: Numbers of statistically significant interaction peaks in six cell lines at 51 informative loci and 12 uninformative loci
    - Supplementary Data 3: Coordinates of interacting pairs detected in at least two cell lines (bedpe, -log10 FDR of interaction significance, cell line, numbed of cells)
    - Supplementary Data 4: Risk loci which formed interaction peaks with target genes in T-47D (T), ZR-75-1 (Z), Bre80 (Br), BT-20 (BT), MDA-MB-231 (M) and GM06990 (G) cell lines. (cytoband, SNP, gene targets).
    - Supplementary Data 5: Distances between published risk SNPs and putative CHi-C target genes (kb) at 36 informative risk loci (cytoband, SNP, hg19 coordinates, gene targets)
    - Supplementary Data 6: eQTL analysis of 69 protein coding target genes at 26 risk loci in TCGA breast cancer data
    - Supplementary Data 7: Disease-specific survival analysis of 97 target genes in Metabric data
    </details> <details>
    <summary>Paper</summary>
    Baxter, Joseph S., Olivia C. Leavy, Nicola H. Dryden, Sarah Maguire, Nichola Johnson, Vita Fedele, Nikiana Simigdala, et al. “Capture Hi-C Identifies Putative Target Genes at 33 Breast Cancer Risk Loci.” Nature Communications 9, no. 1 (December 2018): 1028. https://doi.org/10.1038/s41467-018-03411-9
</details>

- Breast cancer. Epithelial (MCF-10A) and breast cancer (MCF-7) cells. Tumor vs. normal comparison, replicate comparison. Two replicates for each. [GSE66733](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733). The data was reanalyzed in Fritz, Andrew J., Prachi N. Ghule, Joseph R. Boyd, Coralee E. Tye, Natalie A. Page, Deli Hong, David J. Shirley, et al. “[Intranuclear and Higher-Order Chromatin Organization of the Major Histone Gene Cluster in Breast Cancer](https://doi.org/10.1002/jcp.25996).” Journal of Cellular Physiology 233, no. 2 (February 2018) [GSE98552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98552)
    - Barutcu AR, Lajoie BR, McCord RP, Tye CE et al. [Chromatin interaction analysis reveals changes in small chromosome and telomere clustering between epithelial and breast cancer cells](https://doi.org/10.1186/s13059-015-0768-0). Genome Biol 2015 Sep 28;16:214. PMID: 26415882. 

- Breast cancer. T47D-MTLV cell line. 3D response to progesterone, integrative analysis, effect of cutting enzymes. Hi-C at 0h and 1h time points, with different enzymes. RNA-seq and ChIP-Seq available. No replicates. [GSE53463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53463)
    - Le Dily F, Baù D, Pohl A, Vicent GP et al. [Distinct structural transitions of chromatin topological domains correlate with coordinated hormone-induced gene regulation](http://www.genesdev.org/cgi/doi/10.1101/gad.241422.114). Genes Dev 2014 Oct 1;28(19):2151-62. PMID: 25274727. 

- Breast cancer. MCF-7 cell line. 3D response to estrogen, time course (0, 0.5h, 1h, 4h, 24h), replicate comparison. [GSE51687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51687)
    - Tordini, Fabio, Marco Aldinucci, Luciano Milanesi, Pietro Liò, and Ivan Merelli. “[The Genome Conformation As an Integrator of Multi-Omic Data: The Example of Damage Spreading in Cancer](https://doi.org/10.3389/fgene.2016.00194).” Frontiers in Genetics 7 (November 15, 2016). 



# Tissue-specific

- ChIA-PET loops and gene expression in 24 human cell types. RAD21, H3K27ac, RNA-seq. 28% of loops are variable, distinguish cells by tissue of origin, shorter, depleted of housekeeping genes, coincide with different chromatin states. Genes that have more interactions are depleted in housekeeping functions and enriched for pathogenic variants. [Supplementary material](https://www.nature.com/articles/s41586-020-2151-x#Sec59) has hg19 coordinates of [RAD21 peaks](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2151-x/MediaObjects/41586_2020_2151_MOESM4_ESM.xlsx), [Pan-cell type cohesin-mediated chromatin loops](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2151-x/MediaObjects/41586_2020_2151_MOESM5_ESM.xlsx), [H3K27ac peaks](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2151-x/MediaObjects/41586_2020_2151_MOESM6_ESM.xlsx), and more
    - Grubert, Fabian, Rohith Srivas, Damek  V Spacek, Maya Kasowski, Mariana Ruiz-Velasco, Nasa Sinnott-Armstrong, Peyton Greenside, et al. “[Landscape of Cohesin-Mediated Chromatin Loops in the Human Genome](https://doi.org/10.1038/s41586-020-2151-x).” Nature 583, no. 7818 (July 2020)

## ENCODE

Search query for any type of Hi-C data, e.g., [human brain Hi-C](https://www.encodeproject.org/search/?type=Experiment&assay_slims=3D+chromatin+structure&assay_title=Hi-C&organ_slims=brain)

- [Brain microvascular endothelial cell](https://www.encodeproject.org/experiments/ENCSR507AHE/), [GEO GSE105544](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105544)
- [Brain pericytes](https://www.encodeproject.org/experiments/ENCSR323QIP/), [GEO GSE105513](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105513)
- [Astrocyte of the cerebellum](https://www.encodeproject.org/experiments/ENCSR011GNI/), [GEO GSE105194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105194)
- [Neuroblastoma, derived from a bone marrow metastasis, SK-N-DZ treated and not treated with dimethyl sulfoxide](https://www.encodeproject.org/experiments/ENCSR105KFX/), [GEO GSE105275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105275)
- [Neuroepithelioma, derived from metastatis, SK-N-MC](https://www.encodeproject.org/experiments/ENCSR834DXR/), [GEO GSE105914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105914)

## Brain

- Won, Hyejung, Luis de la Torre-Ubieta, Jason L. Stein, Neelroop N. Parikshak, Jerry Huang, Carli K. Opland, Michael J. Gandal, et al. “[Chromosome Conformation Elucidates Regulatory Relationships in Developing Human Brain](https://www.ncbi.nlm.nih.gov/pubmed/27760116).” Nature, (October 27, 2016) - Two brain regions: the cortical and subcortical plate (CP), consisting primarily of post-mitotic neurons and the germinal zone (GZ), containing primarily mitotically active neural progenitors. Three replicates per condition. [GEO GSE77565](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77565). Controlled access.

- Bonev, Boyan, Netta Mendelson Cohen, Quentin Szabo, Lauriane Fritsch, Giorgio L. Papadopoulos, Yaniv Lubling, Xiaole Xu, et al. “[Multiscale 3D Genome Rewiring during Mouse Neural Development](https://doi.org/10.1016/j.cell.2017.09.043).” Cell, (October 2017)
    - Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107. Four HiC replicates in each condition. Mouse embryonic stem cells (ESs), neural progenitors (NPCs), and cortical neurons (CNs), purified NPC and CN populations from neocortex (ncx_NPC, ncx_CN). Replicated RNA-seq and ChIP-seq (H3K4me3, H4K9me3, H3K27ac, H3K36me3).
    - [Bonev-Cavalli_mmc1.xlsx] - Table S1. Summary Statistics for the Datasets, http://www.cell.com/cms/attachment/2111760282/2083800642/mmc1.xlsx

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “[Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation](https://doi.org/10.15252/msb.20156492).” Molecular Systems Biology, (December 23, 2015)
    - [GEO GSE59027](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59027) - mouse embryonic stem cells (ESC), neuronal progenitor cells (NPC) and neurons. Two datasets per cell type, digested using HindIII and NcoI enzymes. [Genomic coordinates for TADs identified from NcoI datasets](http://msb.embopress.org/content/msb/11/12/852/DC5/embed/inline-supplementary-material-5.xls?download=true)

- 5C libraries generated in Beagan et al. in pluripotent mouse ES cells and multipotent neural progenitor cells were downloaded from GEO accession numbers GSM1974095, GSM1974096, GSM1974099, and GSM1974100 (Beagan et al. 2016). [GEO GSE68582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68582)

## Cell lines

- Haarhuis, Judith H.I., Robin H. van der Weide, Vincent A. Blomen, J. Omar Yáñez-Cuna, Mario Amendola, Marjon S. van Ruiten, Peter H.L. Krijger, et al. “[The Cohesin Release Factor WAPL Restricts Chromatin Loop Extension](https://doi.org/10.1016/j.cell.2017.04.013).” Cell, (May 2017) - WAPL, cohesin's antagonist, DNA release factor, restricts loop length and prevents looping between incorrectly oriented CTCF sites. Together with SCC2/SCC4 complex, WAPL promotes correct assembly of chromosomal structures. WAPL WT and KO Hi-C, RNA-seq, ChIP-seq for CTCF and SMC1. Also, SCC4 KO and combined SCC4-WAPL KO Hi-C. Potential role of WAPL in mitosis chromosome condensation. Tools: HiC-Pro processing, HICCUPS, HiCseq, DI, SomaticSniper for variant calling. Data (Hi-C in custom paired BED format) : [GEO GSE95015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95015)

- Grubert, Fabian, Judith B. Zaugg, Maya Kasowski, Oana Ursu, Damek V. Spacek, Alicia R. Martin, Peyton Greenside, et al. “[Genetic Control of Chromatin States in Humans Involves Local and Distal Chromosomal Interactions](https://doi.org/10.1016/j.cell.2015.07.048).” Cell, (August 2015) - seven Hi-C replicates on Gm12878 cell line, [GEO GSE62742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62742)

- Naumova, Natalia, Maxim Imakaev, Geoffrey Fudenberg, Ye Zhan, Bryan R. Lajoie, Leonid A. Mirny, and Job Dekker. “[Organization of the Mitotic Chromosome](https://doi.org/10.1126/science.1236083).” Science (New York, N.Y.), (November 22, 2013) - E-MTAB-1948 - 5C and Hi-C chromosome conformation capture study on metaphase chromosomes from human HeLa, HFF1 and K562 cell lines across the cell cycle. Two biological and two technical replicates. [ArrayExpress E-MTAB-1948](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1948/samples/)

- Jessica Zuin et al., “[Cohesin and CTCF Differentially Affect Chromatin Architecture and Gene Expression in Human Cells](https://doi.org/10.1073/pnas.1317788111),” Proceedings of the National Academy of Sciences of the United States of America, (January 21, 2014) - CTCF and cohesin (RAD21 protein) are enriched in TAD boundaries. Depletion experiments. Different effect on inter- and intradomain interactions. Loss of cohesin leads to loss of local interactions, but TADs remained. Loss of CTCF leads to both loss of local and increase in inter-domain interactions. Different gene expression changes. TAD structures remain largely intact. Data: Hi-C, RNA-seq, RAD21 ChIP-seq for control and depleted RAD21 and CTCF in HEK293 hepatocytes. Two replicates in each condition. [GEO GSE44267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44267)

## Non-human data

- Erythrocytes 3D genome organization in ten species at the last nucleated stages of maturation (newly generated mouse erythroblasts data and previously generated public blood Hi-C data from other organisms). Lack loops and TADs, strong second diagonal pattern. [Raw data at SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA666472/). <details>
    <summary>Paper</summary>
    Ryzhkova, Anastasia, Alena Taskina, Anna Khabarova, Veniamin Fishman, and Nariman Battulin. “Erythrocytes 3D Genome Organization in Vertebrates.” Scientific Reports 11, no. 1 (December 2021): 4414. https://doi.org/10.1038/s41598-021-83903-9.
</details>

- RNA-seq, ATAC-seq, ChIP-seq, whole genome methylation (30X), Hi-C in 11 adult and two embryonic tissues on zebrafish. Comparison with human and mouse regulatory elements. Enrichment of evolutionary breakpoints at TAD boundaries, H3K4me3 and CCTF signal.De novo chr4 assembly (sex chromosome). scATAC-seq on zebrafish brain - 25 cell types. [GEO GSE134055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134055), [Tweet](https://twitter.com/yuefeng_1/status/1331664006881439744?s=20)
    - Yang, Hongbo, Yu Luan, Tingting Liu, Hyung Joo Lee, Li Fang, Yanli Wang, Xiaotao Wang, et al. “[A Map of Cis-Regulatory Elements and 3D Genome Structures in Zebrafish](https://doi.org/10.1038/s41586-020-2962-9).” Nature, November 25, 2020. 

- tagHi-C protocol for low-input tagmentation-based Hi-C. Applied to mouse hematopoiesis 10 major blood cell types. Changes in compartments and the Rabl configuration defining chromatin condensation. Gene-body-associating domains are a general property of highly-expressed genes. Spatial chromatin loops link GWAS SNPs to candidate blood-phenotype genes. HiC-Pro to Juicer. [GEO GSE142216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142216) - RNA-seq, replicates, [GEO GSE152918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152918) - tagHi-C data, replicates, combined .hic files 
    - Zhang C, Xu Z, Yang S, Sun G, Jia L, Zheng Z, Gu Q, Tao W, Cheng T, Li C, Cheng H. "[tagHi-C Reveals 3D Chromatin Architecture Dynamics during Mouse Hematopoiesis](https://doi.org/10.1016/j.celrep.2020.108206)." Cell Rep. 2020 Sep 29 

- Single-nucleus Hi-C data (scHi-C) of 88 Drosophila BG3 cells. 2-5M paired-end reads per cell, 10kb resolution. ORBITA pipeline to eliminate the effect of Phi29 DNA polymerase template switching. Chromatin compartments approx. 1Mb in size, non-hierarchical conserved TADs can be detected. Lots of biology, integration with other omics data. Raw and processed data in .cool format at [GEO GSE131811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131811)
    - Ulianov, Sergey V., Vlada V. Zakharova, Aleksandra A. Galitsyna, Pavel I. Kos, Kirill E. Polovnikov, Ilya M. Flyamer, Elena A. Mikhaleva, et al. “[Order and Stochasticity in the Folding of Individual Drosophila Genomes](https://doi.org/10.1038/s41467-020-20292-z).” Nature Communications 12, no. 1 (December 2021)

- 3D chromatin organization during spermatogenesis, mouse. Meyotic chromosomes in prophase have weak compartmentalization, TADs, loops. Enrichment in near inter-chromosomal interactions (close to diagonal). The X chromosome lacks domain organization during meiotic sex-chromosome inactivation. Concept and formula for evaluation of genomic compartment strength (Methods). [GEO ](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119805) - Hi-C of meiotic pachytene spermatocytes (PS; 2 biological replicates). Other public Hi-C, RNA-seq, ChIP-seq data.
    - Alavattam, Kris G. “[Attenuated Chromatin Compartmentalization in Meiosis and Its Maturation in Sperm Development](https://doi.org/10.1038/s41594-019-0189-y).” Molecular Biology 26 (2019): 17.

- 3D genome rearrangement is uncoupled from gene expression changes. Introduction, references for and against 3D genomics-gene expression links. Drosophila, a "balancer" line with highly rearranged chromosomes. Negligible association can be detected, but changes in genome topology are not predictive of changes in gene expression, loss of long-range interactions has little impact. [Processed data](http://furlonglab.embl.de/data), [GitHub](https://github.com/ajank/balancer-paper). Raw data: [Whole genome](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7510/), [Hi-C](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7512/), [Capture-C](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7513/), [RNA-seq](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7511/) <details>
    <summary>Paper</summary>

    Ghavi-Helm, Yad, Aleksander Jankowski, Sascha Meiers, Rebecca R. Viales, Jan O. Korbel, and Eileen E. M. Furlong. “[Highly Rearranged Chromosomes Reveal Uncoupling between Genome Topology and Gene Expression](https://doi.org/10.1038/s41588-019-0462-3).” Nature Genetics, July 15, 2019. 
</details>

- Global organization of the B cell genome throughout differentiation by the transcription factor Pax5. Mouse splenic CD4+ cells, B cells at various differentiation stages, granulocytes. diffHiC, TADbit, directionality index. Hi-C and RNA-seq data on GEO [GSE99163](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99163). <details>
    <summary>Paper</summary>
    Johanson, Timothy M. “Transcription-Factor-Mediated Supervision of Global Genome Architecture Maintains B Cell Identity.” Nature Immunology 19 (2018): 14. https://doi.org/10.1038/s41590-018-0234-8
</details>

- TADs in Drosophila, Hi-C and RNA-seq in four cell lines of various origin. dCTCF, SMC3, and Su(Hw) are weakly enriched at TAD boundaries. Transcription and active chromatin (H3K27ac, H3K4me1, H3K4me3, H3K36me3, H4K16ac) are associated with TAD boundaries. Also, BEAF-32 and CP190. Hierarchical TADs. Housekeeping genes tend to be near TAD boundaries and in inter-TAD regions. TAD boundary prediction using regression, modeling to associate TADs with bands, investigation of the hierarchy. Heavy use of the Armatus TAD caller. RNA-seq and replicate Hi-C data, high correlation, merged into 20kb resolution.  [GEO GSE69013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69013)
    - Ulianov, Sergey V., Ekaterina E. Khrameeva, Alexey A. Gavrilov, Ilya M. Flyamer, Pavel Kos, Elena A. Mikhaleva, Aleksey A. Penin, et al. “[Active Chromatin and Transcription Play a Key Role in Chromosome Partitioning into Topologically Associating Domains](https://doi.org/10.1101/gr.196006.115).” Genome Research 26, no. 1 (January 2016)

# Differential Hi-C

- Liquid-liquid phase separation (LLPS) in haematological cancers is associated with intrinsically disordered regions (IDRs) of NUP98-HOXA TF chimera and induces CTCF-independent chromatin loops enriched in proto-oncogenes. Many biochemical assays, imaging, mass-spec, ChIP-seq, RNA-seq. All data at [GEO GSE144643](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144643). In situ Hi-C (HEK293FT kidney cells, IDR wild type and mutated, biological and technical replicates) at [GEO GSE143465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143465). <details>
    <summary>Paper</summary>
    Ahn, Jeong Hyun, Eric S. Davis, Timothy A. Daugird, Shuai Zhao, Ivana Yoseli Quiroga, Hidetaka Uryu, Jie Li, et al. “Phase Separation Drives Aberrant Chromatin Looping and Cancer Development.” Nature, June 23, 2021. https://doi.org/10.1038/s41586-021-03662-5.
</details>

- WIZ (widely interspaced zinc finger-containing protein) - new loop-organizing protein, colocalizes with CTCF and cohesin across the genome. Loss of WIZ increases cohesin occupancy and DNA loops. WIZ maintains proper gene expression and stem cell identity. Arima, Juicer. [GEO GSE137285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137285) - RNA-seq, ChIP-seq, Hi-C replicates in WT and WIZdel mouse ESCs.
    - Justice, Megan, Zachary M. Carico, Holden C. Stefan, and Jill M. Dowen. “[A WIZ/Cohesin/CTCF Complex Anchors DNA Loops to Define Gene Expression and Cell Identity](https://doi.org/10.1016/j.celrep.2020.03.067).” Cell Reports 31, no. 2 (April 2020)

- 3D chromatin reorganization during different types of cellular senescence, replicative (RS) and oncogene-induced (OIS over time course). Senescence-associated heterochromatin loci (SAHFs), formed with the help of DNMT1 via regulation of MMGA2 expression. WI38 primary fibroblasts. OIS - gain in long-range contacts. diffHiC analysis, differential regions enriched in H3K9me3. TADkit for 3D modeling, [visualization](https://vre.multiscalegenomics.eu/data_repositories/data_senescence.php). Data (Hi-C replicates, different conditions, timecourse, H3K4me3/H3K9me3/H3K27ac ChIP-seq, RNA-seq) [GEO GSE130306](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130306)
    - Sati, Satish, Boyan Bonev, Quentin Szabo, Daniel Jost, Paul Bensadoun, Francois Serra, Vincent Loubiere, et al. “[4D Genome Rewiring during Oncogene-Induced and Replicative Senescence](https://doi.org/10.1016/j.molcel.2020.03.007).” Molecular Cell, March 2020

- X chromosome sex differences in Drosophila. Male X chromosome has two-fold upregulation of gene expression, more mid/long-range interactions, weaker boundaries marked by BEAF-32, CP190, Chromator, and CLAMP, a dosage compensation complex cofactor. Less negative slope in distance-dependent decay of interactions, less clustered top scoring interactions (more randomness), more open structure overall. Local score differentiator (LSD-score) to call differential TAD boundaries in CNV-independent manner - more non-matching boundaries than autosomes, ~20% appearing and ~35% disappearing boundaries. Enrichment in epigenomic marks identified stronger boundary association with MSL (male-specific lethal complex) and CLAMP binding. Many other experimental observations. hiclib, hicpipe processing. [R implementation of LSD differential TAD analysis](https://bitbucket.org/koustavpal1988/fly_dc_structuralchanges_2018/src/master/), Hi-C data in bedGraph format [GEO GSE94115](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94115), [Tweet](https://twitter.com/3DGenome_KPal/status/1198911861493837830?s=20)
    - Pal, Koustav, Mattia Forcato, Daniel Jost, Thomas Sexton, Cédric Vaillant, Elisa Salviato, Emilia Maria Cristina Mazza, Enrico Lugli, Giacomo Cavalli, and Francesco Ferrari. “[Global Chromatin Conformation Differences in the Drosophila Dosage Compensated Chromosome X](https://doi.org/10.1038/s41467-019-13350-8).” Nature Communications, (December 2019)

- Hi-C TAD comparison between normal prostate cells (RWPE1) and two prostate cancer cells (C42B, 22Rv1). TADs (TopDom-called) become smaller in cancer, switch epigenetic states. FOXA1 promoter has more loop anchors in cancer. Androgen receptor (AR) locus has chromatin structure changed around it (Figure 6). Loop investigation called with Fit-HiC, motifs (NOMe-seq) enriched in loop-associated enhancers different between normal and cancer. HiTC visualization. Figure 1a, Supplementary Figure 3, 5 - examples/coordinates of TAD boundary/length changes.
- Data For RWPE1, C42B, 22Rv1 cell lines: [GEO GSE118629](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629). In situ Hi-C, 4-cutter MboI,  replicated, text-based sparse matrices at 10kb and 40kb resolution, raw and ICE-normalized, hg19. H3K9me3, H3K27me3, H3K36me3, RNA-seq.
- [Supplementary data](https://www.nature.com/articles/s41467-019-12079-8#Sec22): Data 2 - TAD coordinates and annotations; Data 3 - differentially expressed genes in smaller TADs; Data 4 - gene expression changes in TADs switching epigenomic state; Data 5 - enhancer-promoter loops; Data 6 - coordinates of nucleosome-depleted regions; Data 7 - all differentially expressed genes; Data 8 - target genes of FOXA1-bound enhancers; Data 9 - overexpressed genes with more enhancer-promoter loops
    - Rhie, Suhn Kyong, Andrew A. Perez, Fides D. Lay, Shannon Schreiner, Jiani Shi, Jenevieve Polin, and Peggy J. Farnham. “[A High-Resolution 3D Epigenomic Map Reveals Insights into the Creation of the Prostate Cancer Transcriptome](https://doi.org/10.1038/s41467-019-12079-8).” Nature Communications, (December 2019)
    
- DNA methylation linked with 3D genomics. Methylation directs PRC-dependent 3D organization of mouse ESCs. Hypomethylation in mouse ESCs driven to naive pluripotency in two inhibitors (2i) is accopmanied by redistribution of polycomb H3K27me3 mark and decompaction of chromatin. Focus on HoxC, HoxD loci. Hi-C data processed with distiller and other cool-related tools. RNA-seq, H3K37me3 ChIPseq of Mouse ESCs grown in serum and 2i conditions. Hi-C data in replicates [GEO GSE124342](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124342)
    - McLaughlin, Katy, Ilya M. Flyamer, John P. Thomson, Heidi K. Mjoseng, Ruchi Shukla, Iain Williamson, Graeme R. Grimes, et al. “[DNA Methylation Directs Polycomb-Dependent 3D Genome Re-Organization in Naive Pluripotency](https://doi.org/10.1016/j.celrep.2019.10.031).” Cell Reports 29, no. 7 (November 2019)

- RNA transcription inhibition minimally affects TADs, weakens TAD boundaries. K562, RNAse inhibition before/after crosslinking (bXL/aXL), actinomycin D (complete transcriptional arrest) treatment. Processing using cword, 40kb resolution. Data with replicates of each condition, [GEO GSE114337](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114337)
    - Barutcu, A Rasim, Benjamin J Blencowe, and John L Rinn. “[Differential Contribution of Steady‐state RNA and Active Transcription in Chromatin Organization](https://doi.org/10.15252/embr.201948068).” EMBO Reports 20, no. 10 (October 4, 2019). 

- Comparison of the 3D structure of human and chimpanzee induced puripotent stem cells. Lower-order pairwise interactions are relatively conserved, but higher-order, such as TADs, differ. HiCUP and HOMER for Hi-C data processing to 10kb resolution. cyclic loess normalization, limma for significant interaction definition, Arrowhead on combined replicated wot detect TADs.  Association of differential chromatin interactions with gene expression. PyGenomeTracks for visualization. [Workflowr code](https://ittaieres.github.io/HiCiPSC/), Processed Hi-C data (4 human and 4 chimp iPSCs) [GEO GSE122520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122520)
    - Eres, Ittai E., Kaixuan Luo, Chiaowen Joyce Hsiao, Lauren E. Blake, and Yoav Gilad. “[Reorganization of 3D Genome Structure May Contribute to Gene Regulatory Evolution in Primates](https://doi.org/10.1371/journal.pgen.1008278).”  PLOS Genetics 15, no. 7 (July 19, 2019)

- In situ HiC libraries in biological replicates (n=2) for several hematopoietic celltypes (200mio reads per replicate) with 
a focus on the B cell lineage in mice. The authors investigate the role of the transcription factor Pax5 towards its supervisiory role of organizing the 3D genome architecture throughout B cell differentiation. The raw data are available via [GEO GSE99151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99151)
    - Timothy M. Johanson, Aaron T. L. Lun, Hannah D. Coughlan, Tania Tan, Gordon K. Smyth, Stephen L. Nutt & Rhys S. Allan. "[Transcription-factor-mediated supervision of global genome architecture maintains B cell identity](https://www.nature.com/articles/s41590-018-0234-8)." Nature Immunology, (2018)

- DNA loop changes during macrophage development (THP-1 monocyte to macrophage development under 72h PMA treatment). In situ Hi-C (pbn reads, 10kb resolution), RNA-seq, ATAC-seq, CTCF and H3K27ac ChIP-seq. Formation of multi-hubs at key macrophage genes. Differential (dynamic, DESeq2-detected) loops are enriched for AP-1, more enriched in H3K27ac, in contrast to static loops. Association between local H3K27ac and transcription level with distal DNA elements with elevated H3K27ac. Very few genes and lower H3K27ac signal in lost loops, more genes and H3K27ac signal in gained loops. Fold changes in H3K27ac signal positively correlate with DNA looping. Macrophage development-specific gene ontology enrichments. Network analysis for multi-loop multi-enhancer activation hubs identification. [GEO GSE96800 ChIP-seq, ATAC-seq, RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96800), [Two Hi-C samples, THP-1 PMA-treated and untreated, SRA PRJNA385337](https://www.ncbi.nlm.nih.gov/bioproject/385337). 
    - [Supplemental material](https://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30603-2#supplementaryMaterial): 
        - `Table S1`. DNA Loops in Untreated THP-1 Cells, 16067. Text, hg19 genomic coordinates, columns: anchor1_chrom anchor1_start anchor1_end anchor2_chrom anchor2_start anchor2_end sample -log10(P) anchor1_strand anchor2_strand
        - `Table S2`. DNA Loops in PMA-Treated THP-1 Cells, 16335.
        - `Table S3`. Differential Loops
    - Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. “[Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development](https://doi.org/10.1016/j.molcel.2017.08.006).” Molecular Cell, (September 2017)


# Timecourse Hi-C

- Vara, Covadonga, Andreu Paytuví-Gallart, Yasmina Cuartero, François Le Dily, Francisca Garcia, Judit Salvà-Castro, Laura Gómez-H, et al. “[Three-Dimensional Genomic Structure and Cohesin Occupancy Correlate with Transcriptional Activity during Spermatogenesis](https://doi.org/10.1016/j.celrep.2019.06.037).” Cell Reports, (July 2019) - 3D structure changes during spermatogenesis in mouse. Hi-C, RNA-seq, CTCF/REC8/RAD21L ChIP-seq. Description of biology of each stage (Fibroblasts, spermatogonia, leptonema/zygonema, pachynema/diplonema, round spermatids, sperm), and A/B compartment and TAD analysis (TADbit, insulation score), data normalized with ICE. Integration with differential expression. Changes in distribution of CTCF and cohesins (REC8 and RAD21L). Key tools: BBDuk (BBMap), TADbit, HiCExplorer, HiCRep, DeepTools. Data (no replicates) [GEO GSE132054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132054)

- Paulsen, Jonas, Tharvesh M. Liyakat Ali, Maxim Nekrasov, Erwan Delbarre, Marie-Odile Baudement, Sebastian Kurscheid, David Tremethick, and Philippe Collas. “[Long-Range Interactions between Topologically Associating Domains Shape the Four-Dimensional Genome during Differentiation](https://doi.org/10.1038/s41588-019-0392-0).” Nature Genetics, April 22, 2019 - Long-range TAD-TAD interactions form cliques (>3 TAD interacting) are enriched in B compartments and LADs, downregulated gene expression. Graph representation of TAD interactions. Quantifying statistical significance of between-TAD interactions. TAD boundaries are conserved. TAD cliques are dynamic. Permutation test preserving distances. Armatus for TAD detection. hiclib for data processing, Juicebox for visualization. Data: Time course differentiation or human adipose stem cells (day 0, 1, and 3). Hi-C (two replicates), Lamin B1 ChIP-seq, H3K9me3. [GEO GSE109924](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109924). Also used mouse ES differentiation (Bonev 2017), mouse B cell reprogramming (Stadhouders 2018), scHi-C (Nagano 2017)

- Du, Zhenhai, Hui Zheng, Bo Huang, Rui Ma, Jingyi Wu, Xianglin Zhang, Jing He, et al. “[Allelic Reprogramming of 3D Chromatin Architecture during Early Mammalian Development](https://doi.org/10.1038/nature23263).” Nature, (12 2017) - Developmental time course Hi-C. Mouse early development. low-input Hi-C technology (sisHi-C). TADs are initially absent, then gradually appeared. HiCPro mapping, Pearson correlation on low-resolution matrices, allele resolving. Data:  [GEO GSE82185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82185)

- Hug, Clemens B., Alexis G. Grimaldi, Kai Kruse, and Juan M. Vaquerizas. “[Chromatin Architecture Emerges during Zygotic Genome Activation Independent of Transcription](https://doi.org/10.1016/j.cell.2017.03.024).” Cell, (06 2017) - TADs appearing during zygotic genome activation, independent of transcription. TAD boundaries are enriched in housekeeping genes, colocalize in 3D. Drosophila. Insulation score for boundary detection. Overlap analysis of TAD boundaries. [Processed Hi-C matrices at 5kb resolution (replicates merged, .cool format) and TAD boundaries at nuclear cycle 12, 13, 14, and 3-4 hours post fertilization](https://github.com/vaquerizaslab/Hug-et-al-Cell-2017-Supp-Site)

- Ke, Yuwen, Yanan Xu, Xuepeng Chen, Songjie Feng, Zhenbo Liu, Yaoyu Sun, Xuelong Yao, et al. “[3D Chromatin Structures of Mature Gametes and Structural Reprogramming during Mammalian Embryogenesis](https://doi.org/10.1016/j.cell.2017.06.029).” Cell, (July 13, 2017) - 3D timecourse changes during embryo development, from zygotic (no TADs, many long-range interactions) to 2-, 4-, 8-cell, blastocyst and E7.5 mature embryos (TADs established after several rounds of DNA replication). A/B compartments associated with un/methylatied CpGs, respectively. PC1, directionality index, insulation score to define compartments and TADs, these metrics increase in magnitude/strength during maturation. Enrichment in CTCF, SMC1, H3K4me3, H3K27ac, H3K9ac, H3K4me1, depletion in H3K9me3, H3K36me3, H3K27me3. The compartment strength is weaker in maternal vs. paternal genomes. Covariance for each gene vs. boundary score across the timecourse. Relative TAD intensity changes. [Hi-C and RNA-seq data at different stages, some replicates](ttp://bigd.big.ac.cn/bioproject/browse/PRJCA000241)

# Promoter-capture Hi-C

- SIPs, super-interactive promoters in five hematopoietic cell types (Erythrocyte, Macrophage/monophage, megakaryocyte, naive CD4 T-cells, Neutrophils). Reanalysis of promoter-capture Hi-C data from Javierre et al., “Lineage-Specific Genome Architecture Links Enhancers and Non-Coding Disease Variants to Target Gene Promoters.” study. CHiCAGO pipeline. Promoter-interacting regions (PIRs) interacting with SIPs are more enriched in cell type-specific ATAC-seq peaks, GWAS variants for relevant cell types. SIP-associated genes are higher expressed in relevant cells. Some SIPs are shared across cell lines. Super-SIPs.
    - [Additional File 1](https://www.biorxiv.org/content/biorxiv/early/2021/03/16/2021.03.15.435494/DC2/embed/media-2.xlsx?download=true) - Cell type-specific SIPs and genes.
    - [Additional File 2](https://www.biorxiv.org/content/biorxiv/early/2021/03/16/2021.03.15.435494/DC3/embed/media-3.xlsx?download=true) - Cell type-specific SIPs and GWAS variants
    - Lagler, Taylor M., Yuchen Yang, Yuriko Harigaya, Vijay G. Sankaran, Ming Hu, Alexander P. Reiner, Laura M. Raffield, Jia Wen, and Yun Li. “[Super Interactive Promoters Provide Insight into Cell Type-Specific Regulatory Networks in Blood Lineage Cell Types](https://doi.org/10.1101/2021.03.15.435494).” Preprint. Genetics, March 16, 2021. 

- Genome-wide maps linking disease variants to genes. Activity-By-Contact (ABC) Model. 72 diseases and complex traits (non-specific, no psychiatric), linking 5046 fine-mapped GWAS signals to 2249 genes. 577 genes influence multiple phenotypes. Nearly half enhancers regulate multiple genes.[Table S7](https://www.biorxiv.org/content/biorxiv/early/2020/09/03/2020.09.01.278093/DC6/embed/media-6.tsv?download=true) - Summary of diseases and traits.[Table S9](https://www.biorxiv.org/content/biorxiv/early/2020/09/03/2020.09.01.278093/DC8/embed/media-8.tsv?download=true) - ABC-Max predictions for 72 diseases and complex traits.
    - Nasser, Joseph, Drew T Bergman, Charles P Fulco, Philine Guckelberger, Benjamin R Doughty, Tejal A Patwardhan, Thouis R Jones, et al. “[Genome-Wide Maps of Enhancer Regulation Connect Risk Variants to Disease Genes](https://doi.org/10.1101/2020.09.01.278093),” bioRxiv, September 03, 2020.

- Promoter-enhancer contacts occur in cohesin-dependent and cohesin-independent manner. Promoter Capture Hi-C on degradation of cohesin (SCC1 subunit) and CTCF (both targeted by auxin-inducible degron and mEGFP reporter) in G1-synchronized HeLa cells. The majority of promoter contacts are lost (associated with transcriptional changes, SLAM-seq) but some are retained and gained. Cohesin-independent promoter contacts interact with active enhancers. Cohesin-dependent interactions are typically longer and associated with CTCF, while cohesin-independent interactions are shorter and associated with active promoters and enhancers. HiCUP, CHiCAGO, Chicdiff. [Processed data](https://osf.io/brzuc/), replicates of promoter-capture Hi-C data [GEO GSE145735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145735), replicates of SLAM-seq data [GEO GSE145734](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145734)
    - Thiecke, Michiel J., Gordana Wutz, Matthias Muhar, Wen Tang, Stephen Bevan, Valeriya Malysheva, Roman Stocsits et al. "[Cohesin-dependent and-independent mechanisms mediate chromosomal contacts between promoters and enhancers](https://doi.org/10.1016/j.celrep.2020.107929)." Cell reports, July 21 (2020)

- Promoter-enhancer predictions in 131 cell types and tissues using the Activity-By-Contact (ABC) Model, based on chromatin state (ATAC-seq) and 3D folding (consensus Hi-C). ABC model assumes an element’s quantitative effect on a gene should depend on its strength as an enhancer (Activity) weighted by how often it comes into 3D contact with the promoter of the gene (Contact), and that the relative contribution of an element on a gene’s expression (as assayed by the proportional decrease in expression following CRISPR-inhibition) should depend on that element’s effect divided by the total effect of all elements. Outperforms distance-based methods, 3D-based only, machine learning approaches. [Enhancer-promoter predictions for GM12878, K562, liver, LNCAP, mESCs, NCCIT cells](https://osf.io/uhnb4/), more at [Engreitz Lab page](https://www.engreitzlab.org/resources/). GitHub repository [broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).
    - Fulco, Charles P., Joseph Nasser, Thouis R. Jones, Glen Munson, Drew T. Bergman, Vidya Subramanian, Sharon R. Grossman, et al. “[Activity-by-Contact Model of Enhancer–Promoter Regulation from Thousands of CRISPR Perturbations](https://doi.org/10.1038/s41588-019-0538-0).” Nature Genetics 51, no. 12 (December 2019)

- Promoter-enhancer interactions. Promoter-capture Hi-C, 27 human cell lines. Well-formatted data and hg19 genomic coordinates [Supplementary material](https://www.nature.com/articles/s41588-019-0494-8#Sec23) and http://www.3div.kr/capture_hic
    - Jung, Inkyung, Anthony Schmitt, Yarui Diao, Andrew J. Lee, Tristin Liu, Dongchan Yang, Catherine Tan, et al. “[A Compendium of Promoter-Centered Long-Range Chromatin Interactions in the Human Genome](https://doi.org/10.1038/s41588-019-0494-8).” Nature Genetics, September 9, 2019

- [Promoter capture Hi-C in 17 blood cell types](https://www.chicp.org/). Chromatin interactions are cell type-specific. >50% interactions are one-to-one. Enriched in H3K27ac and H3K4me1 (active enhancers). GWAS loci enriched in PIRs. Table S3 lists prioritized genes/SNPs, for autoimmune diseases. Used CHiCAGO to identify strongly interacting regions. Data has active promoter-enhancer links. More than 2,500 potential disease-associated genes are linked to GWAS SNPs. https://osf.io/u8tzp/
    - Javierre, Biola M., Oliver S. Burren, Steven P. Wilder, Roman Kreuzhuber, Steven M. Hill, Sven Sewitz, Jonathan Cairns, et al. “[Lineage-Specific Genome Architecture Links Enhancers and Non-Coding Disease Variants to Target Gene Promoters](https://doi.org/10.1016/j.cell.2016.09.037).” Cell, (November 17, 2016)

# Single-cell Hi-C

See [Notes on single-cell Hi-C technologies, tools, and data](https://github.com/mdozmorov/scHiC_notes) repository

# Micro-C

See the [Micro-C](https://github.com/mdozmorov/HiC_tools#micro-c) section in the [HiC_tools](https://github.com/mdozmorov/HiC_tools) repository

# GAM

Genome Architecture Mapping data

- [GEO GSE64881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881) - mouse ES cell GAM matrices at 1Mb and 30kb resolution
    - Beagrie, Robert A., Antonio Scialdone, Markus Schueler, Dorothee C. A. Kraemer, Mita Chotalia, Sheila Q. Xie, Mariano Barbieri, et al. “[Complex Multi-Enhancer Contacts Captured by Genome Architecture Mapping](https://doi.org/10.1038/nature21411).” Nature 543, no. 7646 (23 2017)

# Imaging

- [MERFISH](https://github.com/BogdanBintu/ChromatinImaging) - Super-resolution imaging technology, reconstruction 3D structure in single cells at 30kb resolution, 1.2Mb region of Chr21 in IMR90 cells. Distance maps obtained by microscopy show small distance for loci within, and larger between, TADs. TAD-like structures exist in single cells. 2.5Mb region of Chr21 in HCT116 cells, cohesin depletion does not abolish TADs, only alter their preferential positioning. Multi-point (triplet) interactions are prevalent. TAD boundaries are highly heterogeneous in single cells. , diffraction-limited and STORM (stochastic optical reconstruction microscopy) imaging. [GitHub](https://github.com/BogdanBintu/ChromatinImaging)
    - Bintu, Bogdan, Leslie J. Mateo, Jun-Han Su, Nicholas A. Sinnott-Armstrong, Mirae Parker, Seon Kinrot, Kei Yamaya, Alistair N. Boettiger, and Xiaowei Zhuang. “[Super-Resolution Chromatin Tracing Reveals Domains and Cooperative Interactions in Single Cells](https://doi.org/10.1126/science.aau1783).” Science, (October 26, 2018)

- Single-cell level massively multiplexed FISH (MERFISH, sequential genome imaging) to measure 3D genome structure in context of gene expression and nuclear structures. Approx. 650 loci, 50kb resolution, on chr21 10.4-46.7Mb from the hg38 genome assembly, IMR90 cells, population average from approx. 12K chr21 copies, multiple rounds of hybridization. Investigation of TADs, A/B compartments, 87% agreement with bulk Hi-C. Association with cell type markers, transcription. Genome-scale imaging using barcodes, 1041 30kb loci covering autosomes and chrX of IMR90, over 5K cells, 5 replicates. [Processed multiplexed FISH data and more, TXT format](https://zenodo.org/record/3928890), [GitHub](https://github.com/ZhuangLab/Chromatin_Analysis_2020_cell)
    - Su, Jun-Han, Pu Zheng, Seon S. Kinrot, Bogdan Bintu, and Xiaowei Zhuang. “[Genome-Scale Imaging of the 3D Organization and Transcriptional Activity of Chromatin](https://doi.org/10.1016/j.cell.2020.07.032).” Cell, August 2020

- [Parser of multiplexed single-cell imaging data from Bintu et al. 2018 and Su et al. 2020](https://github.com/agalitsyna/DPDchrom_input_parser) - Take 3D coordinates of the regions as input and write the distance and contact matrices for these datasets.

# CTCF

[Notes on CTCF motifs and data](CTCF/README.md)

# Integrative Hi-C

- 3D structure mediates the effect of genetic variants on gene expression. 317 lymphoblastoid (LCL) and 78 fibroblast (FIB) cell lines, Hi-C data from Rao et al. 2014 paper. Regulatory elements identified from H3K4me1, H3K4me3, H3K27ac ChIP-seq. The regulatory activity is structured in 12,583 well-delimited cis-regulatory domains (CRDs) that respect the local chromatin organization into topologically associating domains (TADs) but constitute finer organization. 30 trans-regulatory hubs (TRHs) formed by CDRs on distinct chromosomes, associated with AB compartments and allelic regulation. [Processed data](https://zenodo.org/record/2572871) - cQTLs - variants associated with chromatin peak activity; (cis/trans) eQTLs - variants associated with gene expression; aCRD-QTLs - variants associated with CRD activity; sCRD-QTLs - variants associated with CRD structure; chromatin peaks, and CRDs. For LCL and FIB cell lines, coordinates in hg19. <details>
    <summary>Paper</summary>
    Delaneau, O., M. Zazhytska, C. Borel, G. Giannuzzi, G. Rey, C. Howald, S. Kumar, et al. “Chromatin Three-Dimensional Interactions Mediate Genetic Effects on Gene Expression.” Science (New York, N.Y.) 364, no. 6439 (03 2019). https://doi.org/10.1126/science.aat8266.
</details>

- Du, Qian, Grady C. Smith, Phuc Loi Luu, James M. Ferguson, Nicola J. Armstrong, C. Elizabeth Caldon, Elyssa M. Campbell, et al. “[DNA Methylation Is Required to Maintain Both DNA Replication Timing Precision and 3D Genome Organization Integrity](https://doi.org/10.1016/j.celrep.2021.109722).” Cell Reports, (September 2021)
    - Link between DNA methylation, 3D genome regulation, replication timing. Hypomethylation, shift in partially methylated domain boundaries, are associated with the disruption of 3D genome compartmentalization, replication timing heterogeneity, loss of allele-specific replication, change in gene expression. Wide (non-canonical) H3K4me3-H3K8me3 domains protect late replication. [GEO GSE158011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158011) -  Repli-seq, in situ Hi-C in HCT116 and DKO1 cells, 10X Genomics single cell RNA-seq, Nanopore sequencing; [GEO GSE58638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58638) - Histone ChIP-seq, and more references in Methods. [Scripts](https://github.com/qianxidu/Replication_Timing_Du_et_al_2021) and many tools in Methods.

- Zhang, Ruochi, and Jian Ma. “[MATCHA: Probing Multi-Way Chromatin Interaction with Hypergraph Representation Learning](https://doi.org/10.1016/j.cels.2020.04.004).” Cell Systems 10, no. 5 (May 2020)
    - GM129878-specific [SPRITE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114242), [Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), [single-cell Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84920) data (also, [Repli-seq](https://data.4dnucleome.org/experiments-repliseq/4DNEXI55T28T/)).
    - Drosophila S2 cell line [ChIA-Drop](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109355), [Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99104)

# Misc

- [NucPosDB](https://generegulation.org/nucposdb/) - nucleosome positioning database, nucleomics of cell-free DNA (cfDNA) nucleomics (reflects apoptosis, necrosis, or NETosis processes). Experimental in vivo nucleosome maps, stable nucleosomes of the human genome, cfDNA datasets processed with [NucTools](https://homeveg.github.io/nuctools/), tools for nucleosome data analysis and prediction of nucleosome formation from DNA sequences. Different organisms (more than 16), cell types, diseases, technologies (MNase-seq,H3 ChIP-seq, MH-seq, MPE-seq, MiSeq, NOME-seq, RED-seq), raw/processed (BED-like) data. Average nucleosome profiles differ between healthy and diseas (inflammation, cancer) states.
    - Mariya Shtumpf, Kristan V Piroeva, Shivam P Agrawal, Divya R Jacob, Vladimir B. Teif "[NucPosDB: a database of nucleosome positioning in vivo and nucleosomics of cell-free DNA](https://doi.org/10.1101/2021.11.24.469884)"
bioRxiv 2021.11.24 

- Prioritization of COVID-19 candidate genes using 3D chromosomal topology. Applying COGS (Capture Hi-C Omnibus Gene Score), a statistical pipeline for linking GWAS variants with their target genes based on 3D chromatin interaction data. [COVID-19 GWAS data](https://www.covid19hg.org/results/r5/). Promoter-capture Hi-C data from Javierre et al., “Lineage-Specific Genome Architecture Links Enhancers and Non-Coding Disease Variants to Target Gene Promoters” and Ho et al. "TOP1 inhibition therapy protects against SARS-CoV-2-induced lethal inflammation" studies ([17 human primary cell types data](https://osf.io/u8tzp/) and [SARS-CoV-2-infected lung carcinoma cells data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164533)). Four prioritization approaches, summary in [Supplementary Table S4](https://www.frontiersin.org/articles/10.3389/fgene.2021.745672/full#supplementary-material). Biological analysis. 
    - Thiecke, Michiel J., Emma J. Yang, Oliver S. Burren, Helen Ray-Jones, and Mikhail Spivakov. “[Prioritisation of Candidate Genes Underpinning COVID-19 Host Genetic Traits Based on High-Resolution 3D Chromosomal Topology](https://doi.org/10.3389/fgene.2021.745672).” Frontiers in Genetics 12 (October 25, 2021)

- [Consensus Hi-C matrix, FTP](ftp://ftp.broadinstitute.org/outgoing/lincRNA/average_hic/average_hic.v2.191020.tar.gz), mean of 10 cell line-specific matrices, from Nasser, Joseph, Drew T Bergman, Charles P Fulco, Philine Guckelberger, Benjamin R Doughty, Tejal A Patwardhan, Thouis R Jones, et al. “[Genome-Wide Maps of Enhancer Regulation Connect Risk Variants to Disease Genes](https://doi.org/10.1101/2020.09.01.278093),” bioRxiv, September 03, 2020

- Akdemir, Kadir C. “[Somatic Mutation Distributions in Cancer Genomes Vary with Three-Dimensional Chromatin Structure](https://doi.org/10.1038/s41588-020-0708-0).” Nature Genetics, 05 October 2020 - Spatial distribution of mutations in cancer. 3000 tumor-normal pair whole-genome datasets from 42 human cancer types. Common TAD boundaries from 5 cell types. Sharp increase in somatic mutational load is co-localized with topologically associating domain boundaries, but not CNVs. Distinct mutational signatures affect active/inactive domains. Late-replicating genomic regions and inactive domains acquire a higher mutational load. Kataegis regions overlapped with TAD boundaries and domains with higher transcriptional activity. Various public and newly generated Hi-C data (DNA-repair proficient and deficient colon cancer cell lines) in [Supplementary Table 1](https://www.nature.com/articles/s41588-020-0708-0#Sec23). [GitHub](https://github.com/kcakdemir/MutationalDistribution). [Tweet](https://twitter.com/kcakdemir/status/1313133819038490632?s=20)

- Sauerwald, Natalie, and Carl Kingsford. “[Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types](https://doi.org/10.1093/bioinformatics/bty265).” Bioinformatics (Oxford, England), (July 1, 2018) - Analysis of TAD similarity using variation of information (VI) metric as a local distance measure. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. [Kingsford-Group/localtadsim](https://github.com/Kingsford-Group/localtadsim). 23 human Hi-C datasets, Hi-C Pro processed into 100kb matrices, Armatus to call TADs. 

- Sauerwald, Natalie, Akshat Singhal, and Carl Kingsford. “[Analysis of the Structural Variability of Topologically Associated Domains as Revealed by Hi-C](https://doi.org/10.1093/nargab/lqz008).” NAR Genomics and Bioinformatics, 30 September 2019 - TAD variability among 137 Hi-C samples (including replicates, 69 if not) from 9 studies. HiCrep, Jaccard, TADsim to measure similarity. Variability does not come from genetics. Introduction to TADs. 10-70% of TAD boundaries differ between replicates. 20-80% differ between biological conditions. Much less variation across individuals than across tissue types. Lab -specific source of variation - in situ vs. dilution ligation protocols, restriction enzymes not much. HiCpro to 100kb data, ICE-normalization, Armatus for TAD calling. Table 1 - all studies and accession numbers.

- McCole, Ruth B., Jelena Erceg, Wren Saylor, and Chao-Ting Wu. “[Ultraconserved Elements Occupy Specific Arenas of Three-Dimensional Mammalian Genome Organization](https://doi.org/10.1016/j.celrep.2018.06.031).” Cell Reports, (July 10, 2018) - Ultraconserved elements analysis in the context of 3D genomic structures (TADs, boundaries, loop anchors). Enriched (obseerved/expected overlaps) in domains, depleted in boundaries, no enrichment in loops. Separate analysis for exonic, intronic, intergenic UCEs. Human and mouse Hi-C data. Supplementary tables - coordinates of UCEs, more. [rmccole/UCEs_genome_organization](https://github.com/rmccole/UCEs_genome_organization)
    - `McCole_2018` - [Supplementary material](https://www.cell.com/action/showImagesData?pii=S2211-1247%2818%2930941-0)
    - `mmc2.xlsx` - Table S1. Hi-C Datasets, genomic coordinates of human/mouse pooled domains/boundaries, cell-specific domains/boundaries
    - `mmc3.xlsx` - Table S2. Depletion/Enrichment Analysis. "C" and "D" sheets have genomic coordinates of hg19/mm9 UCEs and their Intergenic/intronic/exonic subsets.

- Nagano, Takashi, Csilla Várnai, Stefan Schoenfelder, Biola-Maria Javierre, Steven W. Wingett, and Peter Fraser. “[Comparison of Hi-C Results Using in-Solution versus in-Nucleus Ligation](https://doi.org/10.1186/s13059-015-0753-7).” Genome Biology (August 26, 2015) - comparing _in situ_ and _in solution_ HiC ligation protocol. Mouse liver cells and human ES cells. Two biological and two technical replicates. [GEO GSE70181](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70181)

- Trussart, M., F. Serra, D. Bau, I. Junier, L. Serrano, and M. A. Marti-Renom. “[Assessing the Limits of Restraint-Based 3D Modeling of Genomes and Genomic Domains](https://doi.org/10.1093/nar/gkv221).” Nucleic Acids Research, (April 20, 2015) - TADbit - modeling 3D structures from Hi-C data. Hi-C matrix simulation methods. [The contact maps and the underlying structures](http://sgt.cnag.cat/3dg/datasets/) - simulated and real datasets, text files, square interaction matrices

- [Genomic coordinates of replication domains boundaries (mm9, hg19, multiple cell lines), TAD boundaries (hg19, IMR90, 40kb and 20kb resolution)](http://mouseencode.org/publications/mcp05/)









