# CTCF notes

- CTCF is an insulator protein strongly enriched in TAD boundaries.

- CTCF is enriched near TAD boundaries, "near" is defined as +/- 20kb around the boundary (Dixon et.al. 2012)

- Considering linear genome, two convergent CTCF binding sites are very likely to outline a chromatin loop. 

- Cell type specificity of CTCF binding - "Widespread plasticity in CTCF occupancy linked to DNA methylation" [PMID: 22955980]



# How to find genomic coordinates and orientation of CTCF binding?

## Scanning genomic sequence for CTCF motifs
    
- https://www.biostars.org/p/278267/
- CTCFBSDB 2.0: A database for CTCF binding sites and genome organization, predicted CTCF binding sites download from http://insulatordb.uthsc.edu/
- `CTCFBSDB_PWM.mat` - variants of CTCF motifs, described http://insulatordb.uthsc.edu/help_new.php
- CTCF binds to the consensus sequence CCGCGNGGNGGCAG, from https://en.wikipedia.org/wiki/CTCF. Can scan FASTA file with FIMO, http://meme-suite.org/meme_4.11.2/tools/fimo
- The position frequency matrix of CTCF for human was downloaded from Jaspar 2016 (http://jaspar.genereg.net). CTCF motif occurrences were identified by the FIMO package (V4.11.1) with the P-value < 1e-5. In total, 110879 motif occurrences were identified.
- Motif finder-identified CTCF sites: Potential CTCF motifs across provided genomes are available at http://hicfiles.s3.amazonaws.com/internal/motifs/GENOME_ID.motifs.txt (e.g. http://hicfiles.s3.amazonaws.com/internal/motifs/hg19.motifs.txt). hg19, hg38, mm9, and mm10 supported

## Predicting CTCF motifs

- From Oti, Martin, Jonas Falck, Martijn A. Huynen, and Huiqing Zhou. “CTCF-Mediated Chromatin Loops Enclose Inducible Gene Regulatory Domains.” BMC Genomics 17 (March 22, 2016): 252. https://doi.org/10.1186/s12864-016-2516-6. - CTCF loops investigation in multiple tissues. Max size - 200kb. Enclose regulatory domains of enhancer-regulated genes. Within loops - enrichment in enhancer-related marks. on the boundaries - histone marks and housekeeping genes from Eisenberg E, Levanon EY. Human housekeeping genes, revisited. Predict CTCF loops from ChIP-seq peaks. CTCF orientation method - should be oriented into the loop. Predicted CTCF sites: https://zenodo.org/record/29423
    - `fimo_ctcfmotifs_MA0139_hg19_2.5e-4.bed` - genome-wide CTCF motifs in human genome (hg19) detected by FIMO tool. From https://zenodo.org/record/29423. 1310708 CTCF motifs. Columns: chromosome, start, end, name, score, strand, p-value, q-value, sequence.
    - `ctcf_predictedloops_ENCODE_chipseq_datasets.tar.gz` - Predicted CTCF loops for 100 ENCODE ChIP-seq datasets. 100 files named like `predictedloops_wgEncodeAwgTfbsBroadGm12878CtcfUniPk_prop04.bed`. Columns: chromosome, start, end, paired coordinates and score, score, strand as dot, start coordinate of the first in pair, end coordinate of the second in pair, 16711680, 2, comma-separated width of CTCF sites, comma-separated something.

## Experimentally obtained ChIP-seq data  

- CTCF footprinting with MNase HiChIP in K562. Short (<80bp) CTCF-protected fragments and longer (>120bp) nucleosome-protected fragments, alignment with upstream (16bp) and core (19bp) motif parts. Region Capture Micro-C in mESCs to identify CTCF-cohesin occupancy, cohesin depletion with auxin. Fully extruded loops are rare. Active chromatin impedes extrusion. Integration with other ChIP-seq and Hi-C datasets. CTCF Analyzer (with) Multinomial Estimation (CAMEL), a tool to detect significant CTCF footprints at near base-pair resolution. [GitHub](https://github.com/aryeelab/cohesin_extrusion_reproducibility) with detailed scripts for all prepeocessing and analyses. GEO [GSE285087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285087) PAIRS data. <details>
    <summary>Paper</summary>
    Sept, Corriene E., Y. Esther Tak, Viraat Goel, Mital S. Bhakta, Christian G. Cerda-Smith, Haley M. Hutchinson, Marco Blanchette, et al. “High-Resolution CTCF Footprinting Reveals Impact of Chromatin State on Cohesin Extrusion.” Nature Communications 16, no. 1 (May 15, 2025): 4506. https://doi.org/10.1038/s41467-025-57775-w.
</details>

- CTCF clusters at TAD boundaries, over extended genomic intervals, CTCF clusters correlate with insulation score. Nano-C and 4C-seq on mESC cells, detailed dissection of CTCF clustering and contribution to domain boundary formation. [Processed ChIP-seq, Nano-C, 4C-seq data](https://data.mendeley.com/datasets/g7b4z8957z/1). [Supplementary material](https://www.biorxiv.org/content/10.1101/2021.04.15.440007v1.supplementary-material): CTCF ChIP-seq peaks in mESCs, over 83K peaks with at least one significant CTCF binding motif, mm10, Extended Data Table 1. Coordinates of TADs, genome-wide insulation scores and genome-wide derivative insulation scores are provided in Extended Data Tables 3-5.
    - Chang, Li-Hsin, Sourav Ghosh, Andrea Papale, Mélanie Miranda, Vincent Piras, Jéril Degrouard, Mallory Poncelet, et al. “[A Complex CTCF Binding Code Defines TAD Boundary Structure and Function](https://doi.org/10.1101/2021.04.15.440007).” Preprint. Genetics, April 15, 2021.

- Previously published ChIP-seq data for CTCF from mouse, macaque, and dog livers (Schmidt et al., 2012, http://www.cell.com/cell/abstract/S0092-8674(11)01507-8)

- Oriented CTCFs for GM12878, from Aidenlab https://aidenlab.org/data.html, https://hicfiles.s3.amazonaws.com/external/GM12878_CTCF_orientation.bed

- mESC CTCF, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36027

- mESC epigenomic marks, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31039

- RenLab Hi-C and CTCF data http://chromosome.sdsc.edu/mouse/download.html

- Find all available CTCF datasets at http://cistrome.org/db/#/. Need to select relevant tissue-specific marks.

- Cheng Y, Ma Z, Kim B-H, Wu W, Cayting P, Boyle AP, Sundaram V, Xing X, Dogan N, Li J, et al. 2014. Principles of regulatory information conservation between mouse and human. Nature 515: 371–375. - CTCF in human (ENCODE Gm12878 and K562 cells) and mouse (MEL and CH12 cells), data in Table S2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4343047/bin/NIHMS664131-supplement-Table_S2.xlsx

# Other transcription factors

- CEBPB, CMYC, CTCF, JUND, MAFK, P300, POL2, POLR2A, RAD21, SMC3, TAF1, and TBP for hESCs, and CEBPB, CTCF, MAFK, POLR2A, and RAD21 - From Arboretum-Hi-C paper

- In addition to CTCF, ZNF143, YY1, DNAse, H3K36me3, TSSs, RNA Pol II, SP1, ZNF274, SIX5. From Hong, Seungpyo, and Dongsup Kim. “Computational Characterization of Chromatin Domain Boundary-Associated Genomic Elements.” Nucleic Acids Research 45, no. 18 (October 13, 2017): 10403–14. https://doi.org/10.1093/nar/gkx738.

- Four members of the cohesin complex (STAG2, SMC3, SMC1A, and RAD21) - from M. Ryan Corces and Victor G. Corces, “The Three-Dimensional Cancer Genome,” Current Opinion in Genetics & Development 36 (2016)

- Histone mark overlap at boundaries - expected enrichment in H3K4me3, depletion in H3K4me1, Figure 3B in [@wang:2017aa]. Dixon's observations of H3K4me3, H3K27ac enrichment at TAD boundaries, depletion of H3K9me3, H3K27me3 enriched across resolutions. Genes, especially, housekeeping, are enriched near TAD boundaries. High- and extreme occupancy target regions (TF binding) are also enriched. CTCF enrichment is general. See also [@Malik:2015aa]

# Housekeeping genes

- Eisenberg, Eli, and Erez Y. Levanon. “Human Housekeeping Genes, Revisited.” Trends in Genetics: TIG 29, no. 10 (October 2013): 569–74. https://doi.org/10.1016/j.tig.2013.05.010.

- Downloadable list: http://www.tau.ac.il/~elieis/HKG/, http://www.tau.ac.il/~elieis/HKG/HK_genes.txt

- Housekeeping genes from single-cell RNA-seq, human, mouse, http://shiny.maths.usyd.edu.au/scHK/, paper https://www.biorxiv.org/content/early/2017/12/06/229815

# References

- Wit, Elzo de, Erica S. M. Vos, Sjoerd J. B. Holwerda, Christian Valdes-Quezada, Marjon J. A. M. Verstegen, Hans Teunissen, Erik Splinter, Patrick J. Wijchers, Peter H. L. Krijger, and Wouter de Laat. “CTCF Binding Polarity Determines Chromatin Looping.” Molecular Cell 60, no. 4 (November 19, 2015): 676–84. https://doi.org/10.1016/j.molcel.2015.09.023. http://www.cell.com/molecular-cell/abstract/S1097-2765(15)00762-5 - CTCF forward-reverse (convergent) orientation is needed to form loops. Good intro about TADs, boundaries, gene coexpression

- Vietri Rudan, Matteo, Christopher Barrington, Stephen Henderson, Christina Ernst, Duncan T. Odom, Amos Tanay, and Suzana Hadjur. “Comparative Hi-C Reveals That CTCF Underlies Evolution of Chromosomal Domain Architecture.” Cell Reports 10, no. 8 (March 3, 2015): 1297–1309. https://doi.org/10.1016/j.celrep.2015.02.004. http://www.cell.com/cell-reports/abstract/S2211-1247(15)00112-6. - CTCF role in chromatin structure. Convergent orientation is important. CTCF sites and TAD boundaries are conserved across organisms. Short mentioning of CTCF orientation detection with MEME
