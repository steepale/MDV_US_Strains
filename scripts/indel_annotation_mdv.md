Represent MDV Indels in Table
================
Alec Steep
May 10, 2019

We called Indels from different strains of Marek’s Disease Virus with
VarScan2. This script organizes those calls, annotates which genes they
occurred in, and present them in a table for publication.

## Resources, Dependencies, and Environment

Set options and working directory:

``` r
knitr::opts_chunk$set(echo = TRUE)
# Set the working directory
wd = '/Users/Alec/Documents/Bioinformatics/MDV_Project/snp_phylo'
knitr::opts_knit$set(root.dir = wd)
# Comprable to:
# setwd(wd)
```

In case you need to install dependencies:

``` r
#BiocManager::install("EDASeq")
```

Load dependencies:

``` r
library(purrr)
pacs...man <- c("dplyr","tibble","stringr", "data.table","GenomicRanges","GenomicFeatures",
                "Biostrings","BSgenome","AnnotationHub","BSgenome.Ggallus.UCSC.galGal5",
                "org.Gg.eg.db","stringdist",'tidyr','BSgenome.Hsapiens.UCSC.hg19','cluster',
                'biomaRt','ggplot2','GenomicRanges')
purrr::walk(pacs...man, .f = function(X) {
        do.call("library", list(X)) 
})
```

#### Functions

Make the ‘not in’ operator:

``` r
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
```

## Load Indels and Annotate

#### Create a BSgenome object

``` r
mdv_fa <- './data/Md5B40BAC.fasta'
mdv_seq = readDNAStringSet(mdv_fa)[[1]]
mdv_string <- readDNAStringSet(mdv_fa)
```

#### Load indel calls from each sample into a dataframe

``` r
samples = c('686_wt_TAAGGCGA-TAGATCGC_L001_R','AC722_3_GCTACGCT-GTAAGGAG_L001_R','AC730_S6_3A_CGAGGCTG-ACTGCATA_L001_R','AC738_B1_AAGAGGCA-AAGGAGTA_L001_R','AC739_CON1_GTAGAGGA-CTAAGCCT_L001_R')
indel_df <- data.frame()
for (sample in samples) {
        # Load in the indels for each virus
        infile = paste0('./data/',sample,'_aln_sorted_mpileup_varscan_indel.txt')
        i_df <- read.table(infile, header = TRUE, sep = '\t')
        head(i_df)
        i_df$Sample <- sample
        indel_df <- rbind(indel_df, i_df)
}
head(indel_df)
```

    ##                         Chrom Position Ref    Cons Reads1 Reads2 VarFreq
    ## 1 gi|306447494|gb|HQ149525.1|      543   T   -A/-A     96    530  84.13%
    ## 2 gi|306447494|gb|HQ149525.1|      770   T    */+G    170     21  10.94%
    ## 3 gi|306447494|gb|HQ149525.1|     4302   T    */+A   1002    897   38.6%
    ## 4 gi|306447494|gb|HQ149525.1|     5280   T -TA/-TA    162   1138   87.4%
    ## 5 gi|306447494|gb|HQ149525.1|     7007   T   */-CC    188    331  62.45%
    ## 6 gi|306447494|gb|HQ149525.1|     7839   C    */-G    563     84  12.12%
    ##   Strands1 Strands2 Qual1 Qual2        Pvalue MapQual1 MapQual2 Reads1Plus
    ## 1        2        2    35    34 2.032596e-254        1        1         86
    ## 2        2        2    39    37  2.666539e-07        1        1        166
    ## 3        2        2    32    32  0.000000e+00        1        1        691
    ## 4        2        2    34    34  0.000000e+00        1        1        133
    ## 5        2        2    33    29 3.117920e-135        1        1        138
    ## 6        2        2    34    34  2.896042e-27        1        1        517
    ##   Reads1Minus Reads2Plus Reads2Minus VarAllele
    ## 1          10        156         374        -A
    ## 2           4         19           2        +G
    ## 3         311        190         707        +A
    ## 4          29        470         668       -TA
    ## 5          50        107         224       -CC
    ## 6          46         49          35        -G
    ##                            Sample
    ## 1 686_wt_TAAGGCGA-TAGATCGC_L001_R
    ## 2 686_wt_TAAGGCGA-TAGATCGC_L001_R
    ## 3 686_wt_TAAGGCGA-TAGATCGC_L001_R
    ## 4 686_wt_TAAGGCGA-TAGATCGC_L001_R
    ## 5 686_wt_TAAGGCGA-TAGATCGC_L001_R
    ## 6 686_wt_TAAGGCGA-TAGATCGC_L001_R

#### Annotate indels in GRanges

``` r
indel_gr <- makeGRangesFromDataFrame(indel_df,
                                     keep.extra.columns = TRUE,
                                     ignore.strand=TRUE,
                                     seqnames.field='Chrom',
                                     start.field = 'Position',
                                     end.field = 'Position',
                                     starts.in.df.are.0based=FALSE)
```

#### Load the gene annotation

``` r
infile = "./data/GallidHerpesvirus2_md5_hq149525_annotation.txt"
ann_df <- read.table(infile, header = TRUE, sep = '\t', comment = '#')
ann_df <- ann_df %>% dplyr::select('START','STOP','STRAND','GENE')
ann_df$CHROM <- "gi|306447494|gb|HQ149525.1|"
head(ann_df)
```

    ##   START STOP STRAND    GENE                       CHROM
    ## 1  1762 2477      - R-LORF1 gi|306447494|gb|HQ149525.1|
    ## 2  2858 3081      - R-LORF2 gi|306447494|gb|HQ149525.1|
    ## 3  3139 3446      + R-LORF3 gi|306447494|gb|HQ149525.1|
    ## 4  3946 4178      - R-LORF4 gi|306447494|gb|HQ149525.1|
    ## 5  4154 4493      + R-LORF5 gi|306447494|gb|HQ149525.1|
    ## 6  5664 6280      - R-LORF6 gi|306447494|gb|HQ149525.1|

#### Annotate the annotation in the Granges objects

``` r
ann_gr <- makeGRangesFromDataFrame(ann_df,
                                   keep.extra.columns = TRUE,
                                   ignore.strand=FALSE,
                                   seqnames.field='CHROM',
                                   start.field = 'START',
                                   end.field = 'STOP',
                                   starts.in.df.are.0based=FALSE,
                                   strand.field = 'STRAND')
```

#### Find the overlaps between indels and gene anootation

``` r
olaps <- findOverlaps(ann_gr, indel_gr)
f1 <- factor(subjectHits(olaps),
             levels=seq_len(subjectLength(olaps)))
genes <- splitAsList(mcols(ann_gr)[['GENE']][queryHits(olaps)], f1)
```

#### Data munging

``` r
# Annotate the variants
elementMetadata(indel_gr)[[ 'GENE' ]] <- genes
# Add the column to metadata
mcols(indel_gr)$GENE <- as.character(mcols(indel_gr)$GENE)

# Convert the GRanges object into a dataframe
var_df <- as.data.frame(indel_gr)
head(var_df)
```

    ##                      seqnames start  end width strand Ref    Cons Reads1
    ## 1 gi|306447494|gb|HQ149525.1|   543  543     1      *   T   -A/-A     96
    ## 2 gi|306447494|gb|HQ149525.1|   770  770     1      *   T    */+G    170
    ## 3 gi|306447494|gb|HQ149525.1|  4302 4302     1      *   T    */+A   1002
    ## 4 gi|306447494|gb|HQ149525.1|  5280 5280     1      *   T -TA/-TA    162
    ## 5 gi|306447494|gb|HQ149525.1|  7007 7007     1      *   T   */-CC    188
    ## 6 gi|306447494|gb|HQ149525.1|  7839 7839     1      *   C    */-G    563
    ##   Reads2 VarFreq Strands1 Strands2 Qual1 Qual2        Pvalue MapQual1
    ## 1    530  84.13%        2        2    35    34 2.032596e-254        1
    ## 2     21  10.94%        2        2    39    37  2.666539e-07        1
    ## 3    897   38.6%        2        2    32    32  0.000000e+00        1
    ## 4   1138   87.4%        2        2    34    34  0.000000e+00        1
    ## 5    331  62.45%        2        2    33    29 3.117920e-135        1
    ## 6     84  12.12%        2        2    34    34  2.896042e-27        1
    ##   MapQual2 Reads1Plus Reads1Minus Reads2Plus Reads2Minus VarAllele
    ## 1        1         86          10        156         374        -A
    ## 2        1        166           4         19           2        +G
    ## 3        1        691         311        190         707        +A
    ## 4        1        133          29        470         668       -TA
    ## 5        1        138          50        107         224       -CC
    ## 6        1        517          46         49          35        -G
    ##                            Sample    GENE
    ## 1 686_wt_TAAGGCGA-TAGATCGC_L001_R    <NA>
    ## 2 686_wt_TAAGGCGA-TAGATCGC_L001_R    <NA>
    ## 3 686_wt_TAAGGCGA-TAGATCGC_L001_R R-LORF5
    ## 4 686_wt_TAAGGCGA-TAGATCGC_L001_R    <NA>
    ## 5 686_wt_TAAGGCGA-TAGATCGC_L001_R    <NA>
    ## 6 686_wt_TAAGGCGA-TAGATCGC_L001_R    <NA>

``` r
# Adjust sample names 
x <- str_split(var_df$Sample,'-')
# Grab first element of nested list
x <- sapply(x, `[[`, 1)
var_df$Sample <- substr(x,1,nchar(x)-9)

# Compress columns
indel_df <- var_df
indel_df <- unite(indel_df, 
      col = "Reads1:Reads2:VarFreq:Strands1:Strands2:Qual1:Qual2,Pvalue", 
      c("Reads1","Reads2","VarFreq","Strands1","Strands2","Qual1","Qual2","Pvalue"), 
      sep =":", 
      remove = TRUE)
indel_df <- unite(indel_df, 
                  col = "MapQual1:MapQual2:Reads1Plus:Reads1Minus:Reads2Plus:Reads2Minus:VarAllele", 
                  c("MapQual1","MapQual2","Reads1Plus","Reads1Minus","Reads2Plus","Reads2Minus","VarAllele"), 
                  sep =":", 
                  remove = TRUE)

# Adjust the GENE column
indel_df <- indel_df %>% mutate(Gene = ifelse(!is.na(GENE), GENE, '.'))
indel_df <- indel_df %>% dplyr::select('seqnames','start','Ref','Cons',
                          "Reads1:Reads2:VarFreq:Strands1:Strands2:Qual1:Qual2,Pvalue",
                          "MapQual1:MapQual2:Reads1Plus:Reads1Minus:Reads2Plus:Reads2Minus:VarAllele",
                          "Sample","Gene")
names(indel_df) <- c('chr','pos','ref','alt',
                     "reads1:reads2:varfreq:strands1:strands2:qual1:qual2,pvalue",
                     "mapqual1:mapqual2:reads1plus:reads1minus:reads2plus:reads2minus:varallele",
                     "sample","gene")
indel_df <- indel_df %>% dplyr::select('chr','pos','ref','alt',"sample","gene",
                                               "reads1:reads2:varfreq:strands1:strands2:qual1:qual2,pvalue",
                                               "mapqual1:mapqual2:reads1plus:reads1minus:reads2plus:reads2minus:varallele")

# Sort the column by position
indel_df <- indel_df %>% arrange(pos) 

# Replace sample names
indel_df$sample <- str_replace(indel_df$sample, 'AC722_3','AC722')
indel_df$sample <- str_replace(indel_df$sample, 'AC730_S6_3A','AC730')
indel_df$sample <- str_replace(indel_df$sample, 'AC738_B1','AC738')
indel_df$sample <- str_replace(indel_df$sample, 'AC739_CON1','AC739')
indel_df$sample <- str_replace(indel_df$sample, '686_wt','686')

# Adjust the indel ref and alt fields
head(indel_df)
```

    ##                           chr pos ref   alt sample gene
    ## 1 gi|306447494|gb|HQ149525.1| 543   T -A/-A    686    .
    ## 2 gi|306447494|gb|HQ149525.1| 543   T -A/-A  AC722    .
    ## 3 gi|306447494|gb|HQ149525.1| 543   T -A/-A  AC730    .
    ## 4 gi|306447494|gb|HQ149525.1| 543   T -A/-A  AC738    .
    ## 5 gi|306447494|gb|HQ149525.1| 543   T -A/-A  AC739    .
    ## 6 gi|306447494|gb|HQ149525.1| 770   T  */+G    686    .
    ##   reads1:reads2:varfreq:strands1:strands2:qual1:qual2,pvalue
    ## 1              96:530:84.13%:2:2:35:34:2.03259634144297e-254
    ## 2              42:252:85.42%:2:2:35:34:1.87980150404253e-122
    ## 3               14:127:88.81%:2:2:35:34:7.84568865961562e-65
    ## 4                 3:28:90.32%:2:2:33:33:1.28569734915105e-14
    ## 5               59:197:76.65%:2:2:34:34:1.15758610319791e-88
    ## 6                170:21:10.94%:2:2:39:37:2.6665394471759e-07
    ##   mapqual1:mapqual2:reads1plus:reads1minus:reads2plus:reads2minus:varallele
    ## 1                                                      1:1:86:10:156:374:-A
    ## 2                                                        1:1:39:3:60:192:-A
    ## 3                                                         1:1:13:1:34:93:-A
    ## 4                                                           1:1:1:2:9:19:-A
    ## 5                                                        1:1:56:3:56:141:-A
    ## 6                                                         1:1:166:4:19:2:+G

``` r
# Adjust alt names
x <- str_split(indel_df$alt,'/')
# Grab first element of nested list
indel_df$alt <- sapply(x, `[[`, 2)

# Select the final columns before adjustment for convenience with printing
indel_df <- indel_df %>% dplyr::select('chr','pos','ref','alt',"sample","gene")

# Adjust the ref and alt fields to be characters
indel_df$ref <- as.character(indel_df$ref)
indel_df$alt <- as.character(indel_df$alt)

# Adjust the ref field
indel_df <- indel_df %>% 
        mutate(ref = case_when(grepl('-', alt) ~ paste0(ref,str_remove(alt,'-')),
                               grepl('+',alt) ~ as.character(ref)))
# Adjust the alt field
indel_df <- indel_df %>% mutate(alt = case_when(grepl('-', alt) ~ paste0(substr(ref,1,1)),
                                         grepl('+', alt) ~ paste0(ref,str_remove(alt,'\\+'))))

head(indel_df)
```

    ##                           chr pos ref alt sample gene
    ## 1 gi|306447494|gb|HQ149525.1| 543  TA   T    686    .
    ## 2 gi|306447494|gb|HQ149525.1| 543  TA   T  AC722    .
    ## 3 gi|306447494|gb|HQ149525.1| 543  TA   T  AC730    .
    ## 4 gi|306447494|gb|HQ149525.1| 543  TA   T  AC738    .
    ## 5 gi|306447494|gb|HQ149525.1| 543  TA   T  AC739    .
    ## 6 gi|306447494|gb|HQ149525.1| 770   T  TG    686    .

#### Create a final table for publication

``` r
# Annotate each sample for mutation, transform data
# Group the samples by all but samples column, cat the samples column
indel_df <- indel_df %>%
        group_by(chr,pos,ref, alt, gene) %>%
        summarize(sample = toString(sample))

# Add additional columns to generate table grouped by variants across sample cohort
indel_df <- indel_df %>% mutate(wt_686 = ifelse(grepl('686',sample), 1,0))
indel_df <- indel_df %>% mutate(AC722 = ifelse(grepl('AC722',sample), 1,0))
indel_df <- indel_df %>% mutate(AC730 = ifelse(grepl('AC730',sample), 1,0))
indel_df <- indel_df %>% mutate(AC738 = ifelse(grepl('AC738',sample), 1,0))
indel_df <- indel_df %>% mutate(AC739 = ifelse(grepl('AC739',sample), 1,0))

# Remove the sample column
indel_df <- indel_df %>% dplyr::select(-sample)
head(indel_df)
```

    ## # A tibble: 6 x 10
    ## # Groups:   chr, pos, ref, alt [6]
    ##   chr                 pos ref   alt   gene   wt_686 AC722 AC730 AC738 AC739
    ##   <fct>             <int> <chr> <chr> <chr>   <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 gi|306447494|gb|…   543 TA    T     .           1     1     1     1     1
    ## 2 gi|306447494|gb|…   770 T     TG    .           1     1     0     0     0
    ## 3 gi|306447494|gb|…   892 A     AC    .           0     0     1     0     0
    ## 4 gi|306447494|gb|…   896 TG    T     .           0     0     0     0     1
    ## 5 gi|306447494|gb|…  1259 A     AC    .           0     0     0     1     0
    ## 6 gi|306447494|gb|…  1870 C     CA    R-LOR…      0     0     0     0     1

#### Save the dataframe to file

``` r
write.table(indel_df, file = './data/mdv_indels_genes.txt', sep ='\t',
            quote = FALSE, row.names = FALSE)
```

#### Session Info

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggplot2_3.1.1                      
    ##  [2] biomaRt_2.38.0                     
    ##  [3] cluster_2.0.8                      
    ##  [4] BSgenome.Hsapiens.UCSC.hg19_1.4.0  
    ##  [5] tidyr_0.8.3                        
    ##  [6] stringdist_0.9.5.1                 
    ##  [7] org.Gg.eg.db_3.7.0                 
    ##  [8] BSgenome.Ggallus.UCSC.galGal5_1.4.2
    ##  [9] AnnotationHub_2.14.5               
    ## [10] BSgenome_1.50.0                    
    ## [11] rtracklayer_1.42.2                 
    ## [12] Biostrings_2.50.2                  
    ## [13] XVector_0.22.0                     
    ## [14] GenomicFeatures_1.34.8             
    ## [15] AnnotationDbi_1.44.0               
    ## [16] Biobase_2.42.0                     
    ## [17] GenomicRanges_1.34.0               
    ## [18] GenomeInfoDb_1.18.2                
    ## [19] IRanges_2.16.0                     
    ## [20] S4Vectors_0.20.1                   
    ## [21] BiocGenerics_0.28.0                
    ## [22] data.table_1.12.2                  
    ## [23] stringr_1.4.0                      
    ## [24] tibble_2.1.1                       
    ## [25] dplyr_0.8.0.1                      
    ## [26] purrr_0.3.2                        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.0                    bit64_0.9-7                  
    ##  [3] shiny_1.3.1                   assertthat_0.2.1             
    ##  [5] interactiveDisplayBase_1.20.0 BiocManager_1.30.4           
    ##  [7] blob_1.1.1                    GenomeInfoDbData_1.2.0       
    ##  [9] Rsamtools_1.34.1              yaml_2.2.0                   
    ## [11] progress_1.2.0                pillar_1.3.1                 
    ## [13] RSQLite_2.1.1                 lattice_0.20-38              
    ## [15] glue_1.3.1                    digest_0.6.18                
    ## [17] promises_1.0.1                colorspace_1.4-1             
    ## [19] plyr_1.8.4                    htmltools_0.3.6              
    ## [21] httpuv_1.5.1                  Matrix_1.2-17                
    ## [23] XML_3.98-1.19                 pkgconfig_2.0.2              
    ## [25] zlibbioc_1.28.0               xtable_1.8-4                 
    ## [27] scales_1.0.0                  later_0.8.0                  
    ## [29] BiocParallel_1.16.6           withr_2.1.2                  
    ## [31] SummarizedExperiment_1.12.0   lazyeval_0.2.2               
    ## [33] cli_1.1.0                     magrittr_1.5                 
    ## [35] crayon_1.3.4                  mime_0.6                     
    ## [37] memoise_1.1.0                 evaluate_0.13                
    ## [39] fansi_0.4.0                   tools_3.5.0                  
    ## [41] prettyunits_1.0.2             hms_0.4.2                    
    ## [43] matrixStats_0.54.0            munsell_0.5.0                
    ## [45] DelayedArray_0.8.0            compiler_3.5.0               
    ## [47] rlang_0.3.4                   grid_3.5.0                   
    ## [49] RCurl_1.95-4.12               bitops_1.0-6                 
    ## [51] rmarkdown_1.12                gtable_0.3.0                 
    ## [53] DBI_1.0.0                     R6_2.4.0                     
    ## [55] GenomicAlignments_1.18.1      knitr_1.22                   
    ## [57] utf8_1.1.4                    bit_1.1-14                   
    ## [59] stringi_1.4.3                 Rcpp_1.0.1                   
    ## [61] tidyselect_0.2.5              xfun_0.6
