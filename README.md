# Hi-C data

A (continuously updated) collection of references to Hi-C data. Predominantly human/mouse Hi-C data, with replicates. Search for file name to find its description. Issues and/or Pull requests to add other data are welcome!


* [Large collections](#large-collections)
  * [Lieberman-Aiden lab](#lieberman-aiden-lab)
  * [Leonid Mirny lab](#leonid-mirny-lab)
  * [Bing Ren lab](#bing-ren-lab)
  * [Feng Yue lab](#feng-yue-lab)
  * [4D Nucleome Data Portal](#4d-nucleome-data-portal)
* [Cancer](#cancer)
* [Tissue-specific](#tissue-specific)
  * [ENCODE](#encode)
  * [Brain](#brain)
  * [Cell lines](#cell-lines)
* [Single cell Hi-C](#single-cell-hi-c)
* [Misc](#misc)

# Large collections

## Lieberman-Aiden lab

All HiC data released by Lieberman-Aiden group. Links to Amazon storage and GEO studies. http://aidenlab.org/data.html

- Lieberman-Aiden, Erez, Nynke L. van Berkum, Louise Williams, Maxim Imakaev, Tobias Ragoczy, Agnes Telling, Ido Amit, et al. “Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome.” Science (New York, N.Y.) 326, no. 5950 (October 9, 2009): 289–93. https://doi.org/10.1126/science.1181369. Gm12878, K562 cells. HindIII, NcoI enzymes. Two-three replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199

- Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. “A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping.” Cell 159, no. 7 (December 18, 2014): 1665–80. https://doi.org/10.1016/j.cell.2014.11.021. - Human Gm12878, K562, IMR90, NHEC, HeLa cells, Mouse CH12 cells. Different digestion enzymes (HindIII, NcoI, Mbol, DpnII), different dilutions. Up to 35 biological replicates for Gm12878. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
    - `Rao_Aiden_2014_mmc2.xlsx` - Hi-C meta-data. [Source](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/1bcd7dea-7dbe-45af-8664-4ddd7fb54bc6/mmc2.xlsx)

- Sanborn, Adrian L., Suhas S. P. Rao, Su-Chen Huang, Neva C. Durand, Miriam H. Huntley, Andrew I. Jewett, Ivan D. Bochkov, et al. “Chromatin Extrusion Explains Key Features of Loop and Domain Formation in Wild-Type and Engineered Genomes.” Proceedings of the National Academy of Sciences of the United States of America 112, no. 47 (November 24, 2015): E6456-6465. https://doi.org/10.1073/pnas.1518552112. HAP1, derived from chronic myelogenous leukemia cell line. Replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74072
    - `Sanborn_Aiden_2015_st01.xlsx` - Sheet "S6 - Hi-C experiments" contains information about Hi-C experiments. [Source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4664323/bin/pnas.1518552112.st01.xlsx)

- Rao, Suhas S.P., Su-Chen Huang, Brian Glenn St Hilaire, Jesse M. Engreitz, Elizabeth M. Perez, Kyong-Rim Kieffer-Kwon, Adrian L. Sanborn, et al. “Cohesin Loss Eliminates All Loop Domains.” Cell 171, no. 2 (2017): 305–320.e24. https://doi.org/10.1016/j.cell.2017.09.026. - HCT-116 human colorectal carcinoma cells. Timecourse, replicates under different conditions. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104334
    - `Rao_Aiden_2017_mmc1.xlsx` - Table S1. Summary Statistics for the Datasets. [Source](http://www.cell.com/cms/attachment/2112693605/2084046365/mmc1.xlsx)


## Leonid Mirny lab

http://mirnylab.mit.edu/

Data from multiple studies, in one place, in `.cool` format: ftp://cooler.csail.mit.edu/coolers. Convert to any other format with `cooler` https://cooler.readthedocs.io/en/latest/index.html


## Bing Ren lab

http://chromosome.sdsc.edu/mouse/hi-c/download.html

Raw and normalized chromatin interaction matrices and TADs defined with DomainCaller. Mouse ES, cortex, Human ES, IMR90 fibroblasts. Two replicates per condition. GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35156, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43070

- Dixon, Jesse R., Siddarth Selvaraj, Feng Yue, Audrey Kim, Yan Li, Yin Shen, Ming Hu, Jun S. Liu, and Bing Ren. “Topological Domains in Mammalian Genomes Identified by Analysis of Chromatin Interactions.” Nature 485, no. 7398 (April 11, 2012): 376–80. https://doi.org/10.1038/nature11082.

- Jin, Fulai, Yan Li, Jesse R. Dixon, Siddarth Selvaraj, Zhen Ye, Ah Young Lee, Chia-An Yen, Anthony D. Schmitt, Celso A. Espinoza, and Bing Ren. “A High-Resolution Map of the Three-Dimensional Chromatin Interactome in Human Cells.” Nature 503, no. 7475 (November 14, 2013): 290–94. https://doi.org/10.1038/nature12644.

- Schmitt, Anthony D., Ming Hu, Inkyung Jung, Zheng Xu, Yunjiang Qiu, Catherine L. Tan, Yun Li, et al. “A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome.” Cell Reports 17, no. 8 (November 2016): 2042–59. https://doi.org/10.1016/j.celrep.2016.10.061. Normal human cells, brain (dorsolateral prefrontal cortex, hippocampus), adrenal, bladder, lung, ovary, pancreas, etc. 21 human cell lines and primary tissues. Some replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112. Used in `HiCDB` paper https://doi.org/10.1093/nar/gky789


## Feng Yue lab

Classical datasets for TAD idenrification, provided as raw and normalized matrices. http://promoter.bx.psu.edu/hi-c/download.html


## 4D Nucleome Data Portal

https://data.4dnucleome.org/ - downloadable data from key chromosome conformation capture papers. https://www.4dnucleome.org/software.html - alphabetical list of Hi-C software.


# Cancer

- Harewood, Louise, Kamal Kishore, Matthew D. Eldridge, Steven Wingett, Danita Pearson, Stefan Schoenfelder, V. Peter Collins, and Peter Fraser. “Hi-C as a Tool for Precise Detection and Characterisation of Chromosomal Rearrangements and Copy Number Variation in Human Tumours.” Genome Biology 18, no. 1 (December 2017). https://doi.org/10.1186/s13059-017-1253-8. - ten non-replicated Hi-C datasets. Two human lymphoblastoid cell lines with known chromosomal translocations (FY1199 and DD1618), transformed mouse cell line (EKLF), six human brain tumours: five glioblastomas ( GB176, GB180, GB182, GB183 and GB238) and one anaplastic astrocytoma (AA86), a normal human cell line control (GM07017). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81879

- Rickman, David S., T. David Soong, Benjamin Moss, Juan Miguel Mosquera, Jan Dlabal, Stéphane Terry, Theresa Y. MacDonald, et al. “Oncogene-Mediated Alterations in Chromatin Conformation.” Proceedings of the National Academy of Sciences of the United States of America 109, no. 23 (June 5, 2012): 9083–88. https://doi.org/10.1073/pnas.1112570109. Prostate cancer, normal. RWPE1 prostate epithelial cells transfected with GFP or ERG oncogene. Two biological and up to four technical replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37752

- Taberlay, Phillippa C., Joanna Achinger-Kawecka, Aaron T. L. Lun, Fabian A. Buske, Kenneth Sabir, Cathryn M. Gould, Elena Zotenko, et al. “Three-Dimensional Disorganization of the Cancer Genome Occurs Coincident with Long-Range Genetic and Epigenetic Alterations.” Genome Research 26, no. 6 (June 2016): 719–31. https://doi.org/10.1101/gr.201517.115. Cancer, normal Hi-C. Prostate epithelial cells, PC3, LNCaP. Two-three replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73785

- Barutcu AR, Lajoie BR, McCord RP, Tye CE et al. Chromatin interaction analysis reveals changes in small chromosome and telomere clustering between epithelial and breast cancer cells. Genome Biol 2015 Sep 28;16:214. PMID: 26415882. Breast cancer. Epithelial (MCF-10A) and breast cancer (MCF-7) cells. Tumor vs. normal comparison, replicate comparison. Two replicates for each. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733
    - The data was reanalyzed in Fritz, Andrew J., Prachi N. Ghule, Joseph R. Boyd, Coralee E. Tye, Natalie A. Page, Deli Hong, David J. Shirley, et al. “Intranuclear and Higher-Order Chromatin Organization of the Major Histone Gene Cluster in Breast Cancer.” Journal of Cellular Physiology 233, no. 2 (February 2018): 1278–90. https://doi.org/10.1002/jcp.25996. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98552


- Le Dily F, Baù D, Pohl A, Vicent GP et al. Distinct structural transitions of chromatin topological domains correlate with coordinated hormone-induced gene regulation. Genes Dev 2014 Oct 1;28(19):2151-62. PMID: 25274727. Breast cancer. T47D-MTLV cell line. 3D response to progesterone, integrative analysis, effect of cutting enzymes. Hi-C at 0h and 1h time points, with different enzymes. RNA-seq and ChIP-Seq available. No replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53463

- Tordini, Fabio, Marco Aldinucci, Luciano Milanesi, Pietro Liò, and Ivan Merelli. “The Genome Conformation As an Integrator of Multi-Omic Data: The Example of Damage Spreading in Cancer.” Frontiers in Genetics 7 (November 15, 2016). https://doi.org/10.3389/fgene.2016.00194. Breast cancer. MCF-7 cell line. 3D response to estrogen, time course (0, 0.5h, 1h, 4h, 24h), replicate comparison. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51687


# Tissue-specific

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
    - `Bonev-Cavalli_mmc1.xlsx` - Table S1. Summary Statistics for the Datasets, http://www.cell.com/cms/attachment/2111760282/2083800642/mmc1.xlsx

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation.” Molecular Systems Biology 11, no. 12 (December 23, 2015): 852–852. https://doi.org/10.15252/msb.20156492.
    - `GSE59027` - mouse embryonic stem cells (ESC), neuronal progenitor cells (NPC) and neurons. Two datasets per cell type, digested using HindIII and NcoI enzymes. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59027. Genomic coordinates for TADs identified from NcoI datasets are provided in http://msb.embopress.org/content/msb/11/12/852/DC5/embed/inline-supplementary-material-5.xls?download=true

- 5C libraries generated in Beagan et al. in pluripotent mouse ES cells and multipotent neural progenitor cells were downloaded from GEO accession numbers GSM1974095, GSM1974096, GSM1974099, and GSM1974100 (Beagan et al. 2016). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68582

## Cell lines

- Grubert, Fabian, Judith B. Zaugg, Maya Kasowski, Oana Ursu, Damek V. Spacek, Alicia R. Martin, Peyton Greenside, et al. “Genetic Control of Chromatin States in Humans Involves Local and Distal Chromosomal Interactions.” Cell 162, no. 5 (August 2015): 1051–65. https://doi.org/10.1016/j.cell.2015.07.048. - seven Hi-C replicates on Gm12878 cell line, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62742

- Naumova, Natalia, Maxim Imakaev, Geoffrey Fudenberg, Ye Zhan, Bryan R. Lajoie, Leonid A. Mirny, and Job Dekker. “Organization of the Mitotic Chromosome.” Science (New York, N.Y.) 342, no. 6161 (November 22, 2013): 948–53. https://doi.org/10.1126/science.1236083. - E-MTAB-1948 - 5C and Hi-C chromosome conformation capture study on metaphase chromosomes from human HeLa, HFF1 and K562 cell lines across the cell cycle. Two biological and two technical replicates. https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1948/samples/

- Jessica Zuin et al., “Cohesin and CTCF Differentially Affect Chromatin Architecture and Gene Expression in Human Cells,” Proceedings of the National Academy of Sciences of the United States of America 111, no. 3 (January 21, 2014): 996–1001, https://doi.org/10.1073/pnas.1317788111. - CTCF and cohesin (RAD21 protein) are enriched in TAD boundaries. Depletion experiments. Different effect on inter- and intradomain interactions. Loss of cohesin leads to loss of local interactions, but TADs remained. Loss of CTCF leads to both loss of local and increase in inter-domain interactions. Different gene expression changes. TAD structures remain largely intact. Data: Hi-C, RNA-seq, RAD21 ChIP-seq for control and depleted RAD21 and CTCF in HEK293 hepatocytes. Two replicates in each condition. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44267


# Single cell Hi-C

- Nagano, Takashi, Yaniv Lubling, Tim J. Stevens, Stefan Schoenfelder, Eitan Yaffe, Wendy Dean, Ernest D. Laue, Amos Tanay, and Peter Fraser. “Single-Cell Hi-C Reveals Cell-to-Cell Variability in Chromosome Structure.” Nature 502, no. 7469 (October 3, 2013): 59–64. https://doi.org/10.1038/nature12593. - Single-cell Hi-C. Mouse Th1 cells, 11 samples. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48262

- Takashi Nagano et al., “Cell-Cycle Dynamics of Chromosomal Organization at Single-Cell Resolution,” Nature 547, no. 7661 (July 5, 2017): 61–67, https://doi.org/10.1038/nature23001. Single-cell Hi-C, mouse embryonic stem cells, diploid and haploid, over cell cycle, 45 samples. GEO link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489 and data on the pipeline page https://bitbucket.org/tanaylab/schic2/overview

- Ramani, Vijay, Xinxian Deng, Ruolan Qiu, Kevin L Gunderson, Frank J Steemers, Christine M Disteche, William S Noble, Zhijun Duan, and Jay Shendure. “Massively Multiplex Single-Cell Hi-C.” Nature Methods 14, no. 3 (January 30, 2017): 263–66. https://doi.org/10.1038/nmeth.4155. - Single-cell Hi-C. Cell lines derived from mouse (primary mouse embryonic fibroblasts (MEFs), and the ‘Patski’ embryonic fibroblast line) and human cells (HeLa S3, the HAP1 cell line, K562, and GM12878. Human data - controlled access. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84920

- Stevens, Tim J., David Lando, Srinjan Basu, Liam P. Atkinson, Yang Cao, Steven F. Lee, Martin Leeb, et al. “3D Structures of Individual Mammalian Genomes Studied by Single-Cell Hi-C.” Nature, March 13, 2017. https://doi.org/10.1038/nature21429. - Single-cell Hi-C. Mouse embryonic cell lines. Eight samples. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80280. 100kb five single-cell HiC. TADs are dynamic, A/B compartments, LADs, enhancers/promoters are consistent. 3D clustering of active histone marks, highly expressed genes. Co-expression of genes within TAD boundaries. Videos at http://www.nature.com/nature/journal/v544/n7648/full/nature21429.html#supplementary-information. Supplementary material has processing pipeline description, https://github.com/TheLaueLab/nuc_processing

 


# Misc

- Sauerwald, Natalie, and Carl Kingsford. “Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types.” Bioinformatics (Oxford, England) 34, no. 13 (July 1, 2018): i475–83. https://doi.org/10.1093/bioinformatics/bty265. - Analysis of TAD similarity using variation of information (VI) metric as a local distance measure. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. https://github.com/Kingsford-Group/localtadsim. 23 human Hi-C datasets, Hi-C Pro processed into 100kb matrices, Armatus to call TADs. 

- Sauerwald, Natalie, Akshat Singhal, and Carl Kingsford. “Analysis of the Structural Variability of Topologically Associated Domains as Revealed by Hi-C.” BioRxiv, January 1, 2018, 498972. https://doi.org/10.1101/498972. - TAD variability among 137 Hi-C samples (including replicates, 69 if not) from 9 studies. HiCrep, Jaccard, TADsim to measure similarity. Variability does not come from genetics. Introduction to TADs. 10-70% of TAD boundaries differ between replicates. 20-80% differ between biological conditions. Much less variation across individuals than across tissue types. Lab -specific source of variation - in situ vs. dilution ligation protocols, restriction enzymes not much. HiCpro to 100kb data, ICE-normalization, Armatus for TAD calling. Table 1 - all studies and accession numbers.

- McCole, Ruth B., Jelena Erceg, Wren Saylor, and Chao-Ting Wu. “Ultraconserved Elements Occupy Specific Arenas of Three-Dimensional Mammalian Genome Organization.” Cell Reports 24, no. 2 (July 10, 2018): 479–88. https://doi.org/10.1016/j.celrep.2018.06.031. - Ultraconserved elements analysis in the context of 3D genomic structures (TADs, boundaries, loop anchors). Enriched (obseerved/expected overlaps) in domains, depleted in boundaries, no enrichment in loops. Separate analysis for exonic, intronic, intergenic UCEs. Human and mouse Hi-C data. Supplementary tables - coordinates of UCEs, more. https://github.com/rmccole/UCEs_genome_organization
    - `McCole_2018` - Supplementary material, https://www.cell.com/action/showImagesData?pii=S2211-1247%2818%2930941-0
    - `mmc2.xlsx` - Table S1. Hi-C Datasets, genomic coordinates of human/mouse pooled domains/boundaries, cell-specific domains/boundaries
    - `mmc3.xlsx` Table S2. Depletion/Enrichment Analysis. "C" and "D" sheets have genomic coordinates of hg19/mm9 UCEs and their Intergenic/intronic/exonic subsets.

- Nagano, Takashi, Csilla Várnai, Stefan Schoenfelder, Biola-Maria Javierre, Steven W. Wingett, and Peter Fraser. “Comparison of Hi-C Results Using in-Solution versus in-Nucleus Ligation.” Genome Biology 16 (August 26, 2015): 175. https://doi.org/10.1186/s13059-015-0753-7. - comparing _in situ_ and _in solution_ HiC ligation protocol. Mouse liver cells and human ES cells. Two biological and two technical replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70181

- Trussart, M., F. Serra, D. Bau, I. Junier, L. Serrano, and M. A. Marti-Renom. “Assessing the Limits of Restraint-Based 3D Modeling of Genomes and Genomic Domains.” Nucleic Acids Research 43, no. 7 (April 20, 2015): 3465–77. https://doi.org/10.1093/nar/gkv221. - TADbit - modeling 3D structures from Hi-C data. Hi-C matrix simulation methods. The contact maps and the underlying structures are at http://sgt.cnag.cat/3dg/datasets/ - simulated and real datasets, text files, square interaction matrices

- Genomic coordinates of replication domains boundaries (mm9, hg19, multiple cell lines), TAD boundaries (hg19, IMR90, 40kb and 20kb resolution) http://mouseencode.org/publications/mcp05/

- List of 80 studies (315 Hi-C experiments) from different tissues. Plus 30 extra datasets (Supplementary Table 1). ChIP-seq experiments of histone modification marks, and their QC statistics (Supplementary Table 2). https://academic.oup.com/nar/article/46/D1/D52/4584622#107180936 and http://kobic.kr/3div/statistics. From Yang, Dongchan, Insu Jang, Jinhyuk Choi, Min-Seo Kim, Andrew J. Lee, Hyunwoong Kim, Junghyun Eom, Dongsup Kim, Inkyung Jung, and Byungwook Lee. “3DIV: A 3D-Genome Interaction Viewer and Database.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D52–57. https://doi.org/10.1093/nar/gkx1017.









