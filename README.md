# Hi-C data

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

A (continuously updated) collection of references to Hi-C data. Predominantly human/mouse Hi-C data, with replicates. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Large collections](#large-collections)
  - [Lieberman-Aiden lab](#lieberman-aiden-lab)
  - [Leonid Mirny lab](#leonid-mirny-lab)
  - [Bing Ren lab](#bing-ren-lab)
  - [Feng Yue lab](#feng-yue-lab)
  - [4D Nucleome Data Portal](#4d-nucleome-data-portal)
- [Cancer](#cancer)
- [Tissue-specific](#tissue-specific)
  - [ENCODE](#encode)
  - [Brain](#brain)
  - [Cell lines](#cell-lines)
  - [Non-human data](#non-human-data)
- [Differential Hi-C](#differential-hi-c)
- [Timecourse Hi-C](#timecourse-hi-c)
- [Promoter-enhancer interactions](#promoter-enhancer-interactions)
- [Integrative Hi-C](#integrative-hi-c)
- [CTCF](#ctcf)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Large collections

- [3DIV](http://kobic.kr/3div/download) - database of uniformly processed 315 Hi-C datasets, 80 human cell/tissue types. Bait-centric (SNP rsID, gene name, hg19 coordinates) visualization of long-range interactions in context of epigenomic (histone, enhancers) signals, numerical results. Custom BWA-MEM pipeline, Bias, distance effect removed. Coordinates of significant interactions, with annotations, are available for (FTP) download, http://kobic.kr/3div/download
    - Yang, Dongchan, Insu Jang, Jinhyuk Choi, Min-Seo Kim, Andrew J Lee, Hyunwoong Kim, Junghyun Eom, Dongsup Kim, Inkyung Jung, and Byungwook Lee. “[3DIV: A 3D-Genome Interaction Viewer and Database](https://doi.org/10.1093/nar/gkx1017).” Nucleic Acids Research 46, no. D1 (January 4, 2018)

- [Chorogenome](http://chorogenome.ie-freiburg.mpg.de/) resource: Processed data (Hi-C, ChIP-seq) for Drosophila, Mouse, Human, http://chorogenome.ie-freiburg.mpg.de/
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

[3D Genome Browser](http://3dgenome.fsm.northwestern.edu) - Classical datasets for TAD/loop identification, provided as raw and normalized matrices, genomic coordinates of TADs/loops, tools for various 3C data analysis.


## 4D Nucleome Data Portal

- https://data.4dnucleome.org/ - downloadable data from key chromosome conformation capture papers. 
- https://www.4dnucleome.org/software.html - alphabetical list of Hi-C software.


# Cancer

- Changes in 3D genome are associated with CNVs in multiple myeloma cells (RPMI-8226 trt- and tetraploid, U266 nearly diploid). The number of TADs increases by \~25%, they become smaller, \~20% switch compartment. ICE normalization better accounts for CNVs than HiCNorm. CNV breakpoints overlap with TAD boundaries. 40kb resolution, replicates. [Code](https://github.com/ChengLiLab/myeloma), Hi-C, WGS, RNA-seq data [GSE87585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87585)
    - Wu, Pengze, Tingting Li, Ruifeng Li, Lumeng Jia, Ping Zhu, Yifang Liu, Qing Chen, Daiwei Tang, Yuezhou Yu, and Cheng Li. “[3D Genome of Multiple Myeloma Reveals Spatial Genome Disorganization Associated with Copy Number Variations](https://doi.org/10.1038/s41467-017-01793-w).” Nature Communications 8, no. 1 (December 2017)

- BRCA gene targets regulated by SNPs - Capture-C of chromatin interactions centered on causal variants and promoters of causal genes (Variant- and Promoter Capture Hi-C) in six human mammary epithelial (B80T5, MCF10A) and breast cancer (MCF7, T47D, MDAMB231, Hs578T) cell lines. HindIII fragments, CHiCAGO and Peaky for significant interaction calling. PCA on interactions separates cell types, significant interactions enriched in epigenomic elements. 651 target genes at 139 independent breast cancer risk signals. Table 1 - top priority target genes. [HiCUP-processed capture Hi-C data (hg19)](https://osf.io/2cnw7/), [code](https://github.com/jmbeesley/Beesley_GenomeBiol2019), [Supplementary tables](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1877-y#Sec35), Tables S11 - 651 target genes,
    - Beesley, Jonathan, Haran Sivakumaran, Mahdi Moradi Marjaneh, Luize G. Lima, Kristine M. Hillman, Susanne Kaufmann, Natasha Tuano, et al. “[Chromatin Interactome Mapping at 139 Independent Breast Cancer Risk Signals](https://doi.org/10.1186/s13059-019-1877-y).” Genome Biology 21, no. 1 (December 2020)

- Curtaxins drugs affect 3D genome by DNA intercalation but without inducing DNA damage, compromise enhancer-promoter interactions, suppress oncogene expression, including MYC family genes, downregulates survival genes, partially disrupt TAD borders, decreases short-range interactions, the level of spatial segregation of the A/B compartments, depletes CTCF but not other factors. Hi-C in HT1080 fibrosarcoma cells. Data: Hi-C and CTCF ChIP-seq in duplicates [GSE122463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122463), gene expression in MM1.S and HeLa S3 cells [GSE117611](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117611), H3K27ac [GSE117409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117409), nascent RNA transcription [GSE107633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107633)
    - Kantidze, Omar L., Artem V. Luzhin, Ekaterina V. Nizovtseva, Alfiya Safina, Maria E. Valieva, Arkadiy K. Golov, Artem K. Velichko, et al. “[The Anti-Cancer Drugs Curaxins Target Spatial Genome Organization](https://doi.org/10.1038/s41467-019-09500-7).” Nature Communications 10, no. 1 (December 2019). 

- 3D genomics of glioblastoma. Replicate samples from three patients. Sub-5kb-resolution Hi-C data, integration with ChIP- and RNA-seq. Data: Six Hi-C replicates, [EGAS00001003493](https://ega-archive.org/studies/EGAS00001003493), ChIP-seq [GSE121601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121601), RNA-seq data [EGAS00001003700](https://www.ebi.ac.uk/ega/studies/EGAS00001003700). [Processed data](https://wangftp.wustl.edu/hubs/johnston_gallo/)
    - Johnston, Michael J., Ana Nikolic, Nicoletta Ninkovic, Paul Guilhamon, Florence M.G. Cavalli, Steven Seaman, Franz J. Zemp, et al. “[High-Resolution Structural Genomics Reveals New Therapeutic Vulnerabilities in Glioblastoma](https://doi.org/10.1101/gr.246520.118).” Genome Research 29, no. 8 (August 2019) 

- Ten non-replicated Hi-C datasets. Two human lymphoblastoid cell lines with known chromosomal translocations (FY1199 and DD1618), transformed mouse cell line (EKLF), six human brain tumours: five glioblastomas ( GB176, GB180, GB182, GB183 and GB238) and one anaplastic astrocytoma (AA86), a normal human cell line control (GM07017). [GSE81879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81879)
- Harewood, Louise, Kamal Kishore, Matthew D. Eldridge, Steven Wingett, Danita Pearson, Stefan Schoenfelder, V. Peter Collins, and Peter Fraser. “[Hi-C as a Tool for Precise Detection and Characterisation of Chromosomal Rearrangements and Copy Number Variation in Human Tumours](https://doi.org/10.1186/s13059-017-1253-8).” Genome Biology 18, no. 1 (December 2017). 

- Prostate cancer, normal. RWPE1 prostate epithelial cells transfected with GFP or ERG oncogene. Two biological and up to four technical replicates. [GSE37752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37752)
    - Rickman, David S., T. David Soong, Benjamin Moss, Juan Miguel Mosquera, Jan Dlabal, Stéphane Terry, Theresa Y. MacDonald, et al. “[Oncogene-Mediated Alterations in Chromatin Conformation](https://doi.org/10.1073/pnas.1112570109).” Proceedings of the National Academy of Sciences of the United States of America 109, no. 23 (June 5, 2012)

- Taberlay, Phillippa C., Joanna Achinger-Kawecka, Aaron T. L. Lun, Fabian A. Buske, Kenneth Sabir, Cathryn M. Gould, Elena Zotenko, et al. “[Three-Dimensional Disorganization of the Cancer Genome Occurs Coincident with Long-Range Genetic and Epigenetic Alterations](https://doi.org/10.1101/gr.201517.115).” Genome Research 26, no. 6 (June 2016)
- Cancer, normal Hi-C. Prostate epithelial cells, PC3, LNCaP. Two-three replicates. [GSE73785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73785)

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

Search query for any type of Hi-C data, e.g., human brain, https://www.encodeproject.org/search/?type=Experiment&assay_slims=3D+chromatin+structure&assay_title=Hi-C&organ_slims=brain

- Brain microvascular endothelial cell, https://www.encodeproject.org/experiments/ENCSR507AHE/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105544
- Brain pericytes, https://www.encodeproject.org/experiments/ENCSR323QIP/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105513
- Astrocyte of the cerebellum, https://www.encodeproject.org/experiments/ENCSR011GNI/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105194
- Neuroblastoma, derived from a bone marrow metastasis, SK-N-DZ not treated and treated with dimethyl sulfoxide, https://www.encodeproject.org/experiments/ENCSR105KFX/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105275
- Neuroepithelioma, derived from metastatis, SK-N-MC, https://www.encodeproject.org/experiments/ENCSR834DXR/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105914

## Brain

- Won, Hyejung, Luis de la Torre-Ubieta, Jason L. Stein, Neelroop N. Parikshak, Jerry Huang, Carli K. Opland, Michael J. Gandal, et al. “Chromosome Conformation Elucidates Regulatory Relationships in Developing Human Brain.” Nature 538, no. 7626 (October 27, 2016): 523–27. doi:10.1038/nature19847. https://www.ncbi.nlm.nih.gov/pubmed/27760116. Two brain regions: the cortical and subcortical plate (CP), consisting primarily of post-mitotic neurons and the germinal zone (GZ), containing primarily mitotically active neural progenitors. Three replicates per condition. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77565. Controlled access.

- Bonev, Boyan, Netta Mendelson Cohen, Quentin Szabo, Lauriane Fritsch, Giorgio L. Papadopoulos, Yaniv Lubling, Xiaole Xu, et al. “Multiscale 3D Genome Rewiring during Mouse Neural Development.” Cell 171, no. 3 (October 2017): 557–572.e24. https://doi.org/10.1016/j.cell.2017.09.043.
    - Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107. Four HiC replicates in each condition. Mouse embryonic stem cells (ESs), neural progenitors (NPCs), and cortical neurons (CNs), purified NPC and CN populations from neocortex (ncx_NPC, ncx_CN). Replicated RNA-seq and ChIP-seq (H3K4me3, H4K9me3, H3K27ac, H3K36me3).
    - [Bonev-Cavalli_mmc1.xlsx] - Table S1. Summary Statistics for the Datasets, http://www.cell.com/cms/attachment/2111760282/2083800642/mmc1.xlsx

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation.” Molecular Systems Biology 11, no. 12 (December 23, 2015): 852–852. https://doi.org/10.15252/msb.20156492.
    - [GSE59027] - mouse embryonic stem cells (ESC), neuronal progenitor cells (NPC) and neurons. Two datasets per cell type, digested using HindIII and NcoI enzymes. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59027. Genomic coordinates for TADs identified from NcoI datasets are provided in http://msb.embopress.org/content/msb/11/12/852/DC5/embed/inline-supplementary-material-5.xls?download=true

- 5C libraries generated in Beagan et al. in pluripotent mouse ES cells and multipotent neural progenitor cells were downloaded from GEO accession numbers GSM1974095, GSM1974096, GSM1974099, and GSM1974100 (Beagan et al. 2016). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68582

## Cell lines

- Haarhuis, Judith H.I., Robin H. van der Weide, Vincent A. Blomen, J. Omar Yáñez-Cuna, Mario Amendola, Marjon S. van Ruiten, Peter H.L. Krijger, et al. “The Cohesin Release Factor WAPL Restricts Chromatin Loop Extension.” Cell 169, no. 4 (May 2017): 693-707.e14. https://doi.org/10.1016/j.cell.2017.04.013. - WAPL, cohesin's antagonist, DNA release factor, restricts loop length and prevents looping between incorrectly oriented CTCF sites. Together with SCC2/SCC4 complex, WAPL promotes correct assembly of chromosomal structures. WAPL WT and KO Hi-C, RNA-seq, ChIP-seq for CTCF and SMC1. Also, SCC4 KO and combined SCC4-WAPL KO Hi-C. Potential role of WAPL in mitosis chromosome condensation. Tools: HiC-Pro processing, HICCUPS, HiCseq, DI, SomaticSniper for variant calling. Data (Hi-C in custom paired BED format) : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95015

- Grubert, Fabian, Judith B. Zaugg, Maya Kasowski, Oana Ursu, Damek V. Spacek, Alicia R. Martin, Peyton Greenside, et al. “Genetic Control of Chromatin States in Humans Involves Local and Distal Chromosomal Interactions.” Cell 162, no. 5 (August 2015): 1051–65. https://doi.org/10.1016/j.cell.2015.07.048. - seven Hi-C replicates on Gm12878 cell line, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62742

- Naumova, Natalia, Maxim Imakaev, Geoffrey Fudenberg, Ye Zhan, Bryan R. Lajoie, Leonid A. Mirny, and Job Dekker. “Organization of the Mitotic Chromosome.” Science (New York, N.Y.) 342, no. 6161 (November 22, 2013): 948–53. https://doi.org/10.1126/science.1236083. - E-MTAB-1948 - 5C and Hi-C chromosome conformation capture study on metaphase chromosomes from human HeLa, HFF1 and K562 cell lines across the cell cycle. Two biological and two technical replicates. https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1948/samples/

- Jessica Zuin et al., “Cohesin and CTCF Differentially Affect Chromatin Architecture and Gene Expression in Human Cells,” Proceedings of the National Academy of Sciences of the United States of America 111, no. 3 (January 21, 2014): 996–1001, https://doi.org/10.1073/pnas.1317788111. - CTCF and cohesin (RAD21 protein) are enriched in TAD boundaries. Depletion experiments. Different effect on inter- and intradomain interactions. Loss of cohesin leads to loss of local interactions, but TADs remained. Loss of CTCF leads to both loss of local and increase in inter-domain interactions. Different gene expression changes. TAD structures remain largely intact. Data: Hi-C, RNA-seq, RAD21 ChIP-seq for control and depleted RAD21 and CTCF in HEK293 hepatocytes. Two replicates in each condition. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44267

## Non-human data

- TADs in Drosophila, Hi-C and RNA-seq in four cell lines of various origin. dCTCF, SMC3, and Su(Hw) are weakly enriched at TAD boundaries. Transcription and active chromatin (H3K27ac, H3K4me1, H3K4me3, H3K36me3, H4K16ac) are associated with TAD boundaries. Also, BEAF-32 and CP190. Hierarchical TADs. Housekeeping genes tend to be near TAD boundaries and in inter-TAD regions. TAD boundary prediction using regression, modeling to associate TADs with bands, investigation of the hierarchy. Heavy use of the Armatus TAD caller. RNA-seq and replicate Hi-C data, high correlation, merged into 20kb resolution.  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69013
    - Ulianov, Sergey V., Ekaterina E. Khrameeva, Alexey A. Gavrilov, Ilya M. Flyamer, Pavel Kos, Elena A. Mikhaleva, Aleksey A. Penin, et al. “Active Chromatin and Transcription Play a Key Role in Chromosome Partitioning into Topologically Associating Domains.” Genome Research 26, no. 1 (January 2016): 70–84. https://doi.org/10.1101/gr.196006.115.

# Differential Hi-C

- RNA transcription inhibition minimally affects TADs, weakens TAD boundaries. K562, RNAse inhibition before/after crosslinking (bXL/aXL), actinomycin D (complete transcriptional arrest) treatment. Processing using cword, 40kb resolution. Data with replicates of each condition, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114337
    - Barutcu, A Rasim, Benjamin J Blencowe, and John L Rinn. “[Differential Contribution of Steady‐state RNA and Active Transcription in Chromatin Organization](https://doi.org/10.15252/embr.201948068).” EMBO Reports 20, no. 10 (October 4, 2019). 

- Comparison of the 3D structure of human and chimpanzee induced puripotent stem cells. Lower-order pairwise interactions are relatively conserved, but higher-order, such as TADs, differ. HiCUP and HOMER for Hi-C data processing to 10kb resolution. cyclic loess normalization, limma for significant interaction definition, Arrowhead on combined replicated wot detect TADs.  Association of differential chromatin interactions with gene expression. PyGenomeTracks for visualization. Workflowr code https://ittaieres.github.io/HiCiPSC/, Processed Hi-C data (4 human and 4 chimp iPSCs) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122520
    - Eres, Ittai E., Kaixuan Luo, Chiaowen Joyce Hsiao, Lauren E. Blake, and Yoav Gilad. “[Reorganization of 3D Genome Structure May Contribute to Gene Regulatory Evolution in Primates](https://doi.org/10.1371/journal.pgen.1008278).” Edited by Harmit S. Malik. PLOS Genetics 15, no. 7 (July 19, 2019): e1008278. 

- 3D chromatin reorganization during different types of cellular senescence, replicative (RS) and oncogene-induced (OIS over time course). Senescence-associated heterochromatin loci (SAHFs), formed with the help of DNMT1 via regulation of MMGA2 expression. WI38 primary fibroblasts. OIS - gain in long-range contacts. diffHiC analysis, differential regions enriched in H3K9me3. TADkit for 3D modeling, visualization at https://vre.multiscalegenomics.eu/data_repositories/data_senescence.php. Data (Hi-C replicates, different conditions, timecourse, H3K4me3/H3K9me3/H3K27ac ChIP-seq, RNA-seq) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130306
    - Sati, Satish, Boyan Bonev, Quentin Szabo, Daniel Jost, Paul Bensadoun, Francois Serra, Vincent Loubiere, et al. “4D Genome Rewiring during Oncogene-Induced and Replicative Senescence.” Molecular Cell, March 2020, S1097276520301556. https://doi.org/10.1016/j.molcel.2020.03.007.

- X chromosome sex differences in Drosophila. Male X chromosome has two-fold upregulation of gene expression, more mid/long-range interactions, weaker boundaries marked by BEAF-32, CP190, Chromator, and CLAMP, a dosage compensation complex cofactor. Less negative slope in distance-dependent decay of interactions, less clustered top scoring interactions (more randomness), more open structure overall. Local score differentiator (LSD-score) to call differential TAD boundaries in CNV-independent manner - more non-matching boundaries than autosomes, ~20% appearing and ~35% disappearing boundaries. Enrichment in epigenomic marks identified stronger boundary association with MSL (male-specific lethal complex) and CLAMP binding. Many other experimental observations. hiclib, hicpipe processing. R implementation of LSD differential TAD analysis https://bitbucket.org/koustavpal1988/fly_dc_structuralchanges_2018/src/master/, Hi-C data in bedGraph format https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94115 [Tweet](https://twitter.com/3DGenome_KPal/status/1198911861493837830?s=20)
    - Pal, Koustav, Mattia Forcato, Daniel Jost, Thomas Sexton, Cédric Vaillant, Elisa Salviato, Emilia Maria Cristina Mazza, Enrico Lugli, Giacomo Cavalli, and Francesco Ferrari. “Global Chromatin Conformation Differences in the Drosophila Dosage Compensated Chromosome X.” Nature Communications 10, no. 1 (December 2019): 5355. https://doi.org/10.1038/s41467-019-13350-8.

- DNA methylation linked with 3D genomics. Methylation directs PRC-dependent 3D organization of mouse ESCs. Hypomethylation in mouse ESCs driven to naive pluripotency in two inhibitors (2i) is accopmanied by redistribution of polycomb H3K27me3 mark and decompaction of chromatin. Focus on HoxC, HoxD loci. Hi-C data processed with distiller and other cool-related tools. RNA-seq, H3K37me3 ChIPseq of Mouse ESCs grown in serum and 2i conditions. Hi-C data in replicates https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124342
    - McLaughlin, Katy, Ilya M. Flyamer, John P. Thomson, Heidi K. Mjoseng, Ruchi Shukla, Iain Williamson, Graeme R. Grimes, et al. “DNA Methylation Directs Polycomb-Dependent 3D Genome Re-Organization in Naive Pluripotency.” Cell Reports 29, no. 7 (November 2019): 1974-1985.e6. https://doi.org/10.1016/j.celrep.2019.10.031.

- Hi-C TAD comparison between normal prostate cells (RWPE1) and two prostate cancer cells (C42B, 22Rv1). TADs (TopDom-called) become smaller in cancer, switch epigenetic states. FOXA1 promoter has more loop anchors in cancer. Androgen receptor (AR) locus has chromatin structure changed around it (Figure 6). Loop investigation called with Fit-HiC, motifs (NOMe-seq) enriched in loop-associated enhancers different between normal and cancer. HiTC visualization. Figure 1a, Supplementary Figure 3, 5 - examples/coordinates of TAD boundary/length changes.
- Data For RWPE1, C42B, 22Rv1 cell lines: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629. In situ Hi-C, 4-cutter MboI,  replicated, text-based sparse matrices at 10kb and 40kb resolution, raw and ICE-normalized, hg19. H3K9me3, H3K27me3, H3K36me3, RNA-seq.
- [Supplementary data](https://www.nature.com/articles/s41467-019-12079-8#Sec22): Data 2 - TAD coordinates and annotations; Data 3 - differentially expressed genes in smaller TADs; Data 4 - gene expression changes in TADs switching epigenomic state; Data 5 - enhancer-promoter loops; Data 6 - coordinates of nucleosome-depleted regions; Data 7 - all differentially expressed genes; Data 8 - target genes of FOXA1-bound enhancers; Data 9 - overexpressed genes with more enhancer-promoter loops
    - Rhie, Suhn Kyong, Andrew A. Perez, Fides D. Lay, Shannon Schreiner, Jiani Shi, Jenevieve Polin, and Peggy J. Farnham. “A High-Resolution 3D Epigenomic Map Reveals Insights into the Creation of the Prostate Cancer Transcriptome.” Nature Communications 10, no. 1 (December 2019): 4154. https://doi.org/10.1038/s41467-019-12079-8.



# Timecourse Hi-C

- Du, Zhenhai, Hui Zheng, Bo Huang, Rui Ma, Jingyi Wu, Xianglin Zhang, Jing He, et al. “Allelic Reprogramming of 3D Chromatin Architecture during Early Mammalian Development.” Nature 547, no. 7662 (12 2017): 232–35. https://doi.org/10.1038/nature23263. - Developmental time course Hi-C. Mouse early development. low-input Hi-C technology (sisHi-C). TADs are initially absent, then gradually appeared. HiCPro mapping, Pearson correlation on low-resolution matrices, allele resolving. Data:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82185

- Hug, Clemens B., Alexis G. Grimaldi, Kai Kruse, and Juan M. Vaquerizas. “Chromatin Architecture Emerges during Zygotic Genome Activation Independent of Transcription.” Cell 169, no. 2 (06 2017): 216-228.e19. https://doi.org/10.1016/j.cell.2017.03.024. - TADs appearing during zygotic genome activation, independent of transcription. TAD boundaries are enriched in housekeeping genes, colocalize in 3D. Drosophila. Insulation score for boundary detection. Overlap analysis of TAD boundaries. Processed Hi-C matrices at 5kb resolution (replicates merged, .cool format) and TAD boundaries at nuclear cycle 12, 13, 14, and 3-4 hours post fertilization are at  https://github.com/vaquerizaslab/Hug-et-al-Cell-2017-Supp-Site

- Ke, Yuwen, Yanan Xu, Xuepeng Chen, Songjie Feng, Zhenbo Liu, Yaoyu Sun, Xuelong Yao, et al. “3D Chromatin Structures of Mature Gametes and Structural Reprogramming during Mammalian Embryogenesis.” Cell 170, no. 2 (July 13, 2017): 367-381.e20. https://doi.org/10.1016/j.cell.2017.06.029. - 3D timecourse changes during embryo development, from zygotic (no TADs, many long-range interactions) to 2-, 4-, 8-cell, blastocyst and E7.5 mature embryos (TADs established after several rounds of DNA replication). A/B compartments associated with un/methylatied CpGs, respectively. PC1, directionality index, insulation score to define compartments and TADs, these metrics increase in magnitude/strength during maturation. Enrichment in CTCF, SMC1, H3K4me3, H3K27ac, H3K9ac, H3K4me1, depletion in H3K9me3, H3K36me3, H3K27me3. The compartment strength is weaker in maternal vs. paternal genomes. Covariance for each gene vs. boundary score across the timecourse. Relative TAD intensity changes. Hi-C and RNA-seq data at different stages, some replicates, http://bigd.big.ac.cn/bioproject/browse/PRJCA000241 

- Paulsen, Jonas, Tharvesh M. Liyakat Ali, Maxim Nekrasov, Erwan Delbarre, Marie-Odile Baudement, Sebastian Kurscheid, David Tremethick, and Philippe Collas. “Long-Range Interactions between Topologically Associating Domains Shape the Four-Dimensional Genome during Differentiation.” Nature Genetics, April 22, 2019. https://doi.org/10.1038/s41588-019-0392-0. - Long-range TAD-TAD interactions form cliques (>3 TAD interacting) are enriched in B compartments and LADs, downregulated gene expression. Graph representation of TAD interactions. Quantifying statistical significance of between-TAD interactions. TAD boundaries are conserved. TAD cliques are dynamic. Permutation test preserving distances. Armatus for TAD detection. hiclib for data processing, Juicebox for visualization. Data: Time course differentiation or human adipose stem cells (day 0, 1, and 3). Hi-C (two replicates), Lamin B1 ChIP-seq, H3K9me3. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109924. Also used mouse ES differentiation (Bonev 2017), mouse B cell reprogramming (Stadhouders 2018), scHi-C (Nagano 2017)

- Vara, Covadonga, Andreu Paytuví-Gallart, Yasmina Cuartero, François Le Dily, Francisca Garcia, Judit Salvà-Castro, Laura Gómez-H, et al. “Three-Dimensional Genomic Structure and Cohesin Occupancy Correlate with Transcriptional Activity during Spermatogenesis.” Cell Reports 28, no. 2 (July 2019): 352-367.e9. https://doi.org/10.1016/j.celrep.2019.06.037. - 3D structure changes during spermatogenesis in mouse. Hi-C, RNA-seq, CTCF/REC8/RAD21L ChIP-seq. Description of biology of each stage (Fibroblasts, spermatogonia, leptonema/zygonema, pachynema/diplonema, round spermatids, sperm), and A/B compartment and TAD analysis (TADbit, insulation score), data normalized with ICE. Integration with differential expression. Changes in distribution of CTCF and cohesins (REC8 and RAD21L). Key tools: BBDuk (BBMap), TADbit, HiCExplorer, HiCRep, DeepTools. Data (no replicates) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132054

# Promoter-enhancer interactions

- Promoter-enhancer interactions. Promoter-capture Hi-C, 27 human cell lines. Well-formatted data and hg19 genomic coordinates [Supplementary material](https://www.nature.com/articles/s41588-019-0494-8#Sec23) and http://www.3div.kr/capture_hic
    - Jung, Inkyung, Anthony Schmitt, Yarui Diao, Andrew J. Lee, Tristin Liu, Dongchan Yang, Catherine Tan, et al. “A Compendium of Promoter-Centered Long-Range Chromatin Interactions in the Human Genome.” Nature Genetics, September 9, 2019. https://doi.org/10.1038/s41588-019-0494-8.

- Hi-C promoter capture in 17 blood cell types, sorted. Chromatin interactions are cell type-specific. >50% interactions are one-to-one. Enriched in H3K27ac and H3K4me1 (active enhancers). GWAS loci enriched in PIRs. Table S3 lists prioritized genes/SNPs, for autoimmune diseases. Used CHiCAGO to identify strongly interacting regions. Data has active promoter-enhancer links. More than 2,500 potential disease-associated genes are linked to GWAS SNPs. https://www.chicp.org/
    - Javierre, Biola M., Oliver S. Burren, Steven P. Wilder, Roman Kreuzhuber, Steven M. Hill, Sven Sewitz, Jonathan Cairns, et al. “Lineage-Specific Genome Architecture Links Enhancers and Non-Coding Disease Variants to Target Gene Promoters.” Cell 167, no. 5 (November 17, 2016): 1369-1384.e19. https://doi.org/10.1016/j.cell.2016.09.037.

# CTCF

[Notes on CTCF motifs and data](CTCF/README.md)

# Integrative Hi-C

- Zhang, Ruochi, and Jian Ma. “[MATCHA: Probing Multi-Way Chromatin Interaction with Hypergraph Representation Learning](https://doi.org/10.1016/j.cels.2020.04.004).” Cell Systems 10, no. 5 (May 2020)
    - GM129878-specific [SPRITE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114242), [Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), [single-cell Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84920) data (also, [Repli-seq](https://data.4dnucleome.org/experiments-repliseq/4DNEXI55T28T/)).
    - Drosophila S2 cell line [ChIA-Drop](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109355), [Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99104)

# Misc

- Sauerwald, Natalie, and Carl Kingsford. “Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types.” Bioinformatics (Oxford, England) 34, no. 13 (July 1, 2018): i475–83. https://doi.org/10.1093/bioinformatics/bty265. - Analysis of TAD similarity using variation of information (VI) metric as a local distance measure. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. https://github.com/Kingsford-Group/localtadsim. 23 human Hi-C datasets, Hi-C Pro processed into 100kb matrices, Armatus to call TADs. 

- Sauerwald, Natalie, Akshat Singhal, and Carl Kingsford. “Analysis of the Structural Variability of Topologically Associated Domains as Revealed by Hi-C.” BioRxiv, January 1, 2018, 498972. https://doi.org/10.1101/498972. - TAD variability among 137 Hi-C samples (including replicates, 69 if not) from 9 studies. HiCrep, Jaccard, TADsim to measure similarity. Variability does not come from genetics. Introduction to TADs. 10-70% of TAD boundaries differ between replicates. 20-80% differ between biological conditions. Much less variation across individuals than across tissue types. Lab -specific source of variation - in situ vs. dilution ligation protocols, restriction enzymes not much. HiCpro to 100kb data, ICE-normalization, Armatus for TAD calling. Table 1 - all studies and accession numbers.

- McCole, Ruth B., Jelena Erceg, Wren Saylor, and Chao-Ting Wu. “Ultraconserved Elements Occupy Specific Arenas of Three-Dimensional Mammalian Genome Organization.” Cell Reports 24, no. 2 (July 10, 2018): 479–88. https://doi.org/10.1016/j.celrep.2018.06.031. - Ultraconserved elements analysis in the context of 3D genomic structures (TADs, boundaries, loop anchors). Enriched (obseerved/expected overlaps) in domains, depleted in boundaries, no enrichment in loops. Separate analysis for exonic, intronic, intergenic UCEs. Human and mouse Hi-C data. Supplementary tables - coordinates of UCEs, more. https://github.com/rmccole/UCEs_genome_organization
    - [McCole_2018] - Supplementary material, https://www.cell.com/action/showImagesData?pii=S2211-1247%2818%2930941-0
    - [mmc2.xlsx] - Table S1. Hi-C Datasets, genomic coordinates of human/mouse pooled domains/boundaries, cell-specific domains/boundaries
    - [mmc3.xlsx] - Table S2. Depletion/Enrichment Analysis. "C" and "D" sheets have genomic coordinates of hg19/mm9 UCEs and their Intergenic/intronic/exonic subsets.

- Nagano, Takashi, Csilla Várnai, Stefan Schoenfelder, Biola-Maria Javierre, Steven W. Wingett, and Peter Fraser. “Comparison of Hi-C Results Using in-Solution versus in-Nucleus Ligation.” Genome Biology 16 (August 26, 2015): 175. https://doi.org/10.1186/s13059-015-0753-7. - comparing _in situ_ and _in solution_ HiC ligation protocol. Mouse liver cells and human ES cells. Two biological and two technical replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70181

- Trussart, M., F. Serra, D. Bau, I. Junier, L. Serrano, and M. A. Marti-Renom. “Assessing the Limits of Restraint-Based 3D Modeling of Genomes and Genomic Domains.” Nucleic Acids Research 43, no. 7 (April 20, 2015): 3465–77. https://doi.org/10.1093/nar/gkv221. - TADbit - modeling 3D structures from Hi-C data. Hi-C matrix simulation methods. The contact maps and the underlying structures are at http://sgt.cnag.cat/3dg/datasets/ - simulated and real datasets, text files, square interaction matrices

- Genomic coordinates of replication domains boundaries (mm9, hg19, multiple cell lines), TAD boundaries (hg19, IMR90, 40kb and 20kb resolution) http://mouseencode.org/publications/mcp05/

- List of 80 studies (315 Hi-C experiments) from different tissues. Plus 30 extra datasets (Supplementary Table 1). ChIP-seq experiments of histone modification marks, and their QC statistics (Supplementary Table 2). https://academic.oup.com/nar/article/46/D1/D52/4584622#107180936 and http://kobic.kr/3div/statistics. From Yang, Dongchan, Insu Jang, Jinhyuk Choi, Min-Seo Kim, Andrew J. Lee, Hyunwoong Kim, Junghyun Eom, Dongsup Kim, Inkyung Jung, and Byungwook Lee. “3DIV: A 3D-Genome Interaction Viewer and Database.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D52–57. https://doi.org/10.1093/nar/gkx1017.









