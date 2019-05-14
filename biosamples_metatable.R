############################################################
##### Resources and Dependencies ###########################
############################################################

# Set the working directory
setwd('/Users/Alec/Documents/Bioinformatics/MDV_Project/MDV_comp_Hans_John/data')

# Load the dependencies
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("genoCN", version = "3.8")

# Load dependencies
pacs...man <- c("dplyr","tidyr", "stringr","ggplot2",'tibble') 
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

###########################################################################################
################ Functions ################################################################
###########################################################################################

# Make the 'not in' operator
####################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
####################################################

###########################################################################################
################ BioSamples Metatable #####################################################
###########################################################################################

# Load the sample labels into a dataframe
df <- read.csv(file = 'samples.txt', header = FALSE) 
df <- as_tibble(df)
names(df) <- 'sample_name'

## Add organism (broadcast)
df$organism <- 'Gallid alphaherpesvirus 2'


# Extract isolates, date, loc
sam_split <- read.table(text = as.character(df$sample_name), sep = "_", colClasses = "character")
isolates <- sam_split[[1]]
date <- sam_split[[3]]
geo <- sam_split[[4]]

## Add isolates column
df$isolate <- as.character(isolates)

## Add host column
df$host <- 'Gallus gallus'

## Add lab host column
df$lab_host <- 'Anas platyrhynchos: Duck embryo fibroblast cells'

## Add collection date
df$collection_date <- date

## Adjust the geo info

# Replace states
geo <- replace(geo, geo=='MI', 'USA: Michigan')
geo <- replace(geo, geo=='AL', 'USA: Alabama')
geo <- replace(geo, geo=='CO', 'USA: Colorado')
geo <- replace(geo, geo=='DE', 'USA: Delaware')
geo <- replace(geo, geo=='CO', 'USA: Colorado')
geo <- replace(geo, geo=='NC', 'USA: North Carolina')
geo <- replace(geo, geo=='CA', 'USA: California')
geo <- replace(geo, geo=='IA', 'USA: Iowa')
geo <- replace(geo, geo=='WI', 'USA: Wisconsin')
geo <- replace(geo, geo=='AR', 'USA: Arkansas')
geo <- replace(geo, geo=='MD', 'USA: Maryland')
geo <- replace(geo, geo=='PA', 'USA: Pennsylvania')
geo <- replace(geo, geo=='ME', 'USA: Maine')
geo <- replace(geo, geo=='OH', 'USA: Ohio')
geo <- replace(geo, geo=='NE', 'USA: Nebraska')
geo <- replace(geo, geo=='NY', 'USA: New York')
geo <- replace(geo, geo=='VA', 'USA: Virginia')
geo <- replace(geo, geo=='GA', 'USA: Georgia')
geo <- replace(geo, geo=='MA', 'USA: Massachusetts')

# Add the geo column
df$geo_loc_name <- geo

# Add isolation source column
df$isolation_source <- 'Chicken'

# Save the file
write.table(df, file = 'metadata_seq_sub_mdvR.txt', 
            row.names = FALSE, sep = '\t', quote = FALSE)

###########################################################################################
################ SRA Metatable ###########################################################
###########################################################################################

# Read in file names
file_names <- read.csv(file = 'sra_filenames.txt', header = FALSE)
file_names <- file_names$V1

# Seperate the first files
file1 <- as.character(file_names[grepl("R1", file_names)])
file1 <- as.factor(file1)
# Seperate the second files
file2 <- as.character(file_names[grepl("R2", file_names)])
file2 <- as.factor(file2)

# Generate the SRA dataframe
sra_df <- as.tibble(file1)
names(sra_df) <- 'filename'

# Create a list of samples that underwent WGS
wgs <- c('686_wt_TAAGGCGA-TAGATCGC_L001_R1_001.fastq.gz', 
         '686_wt_TAAGGCGA-TAGATCGC_L001_R2_001.fastq.gz',
         'AC722_3_GCTACGCT-GTAAGGAG_L001_R1_001.fastq.gz',
         'AC722_3_GCTACGCT-GTAAGGAG_L001_R2_001.fastq.gz',
         'AC730_S6_3A_CGAGGCTG-ACTGCATA_L001_R1_001.fastq.gz',
         'AC730_S6_3A_CGAGGCTG-ACTGCATA_L001_R2_001.fastq.gz',
         'AC738_B1_AAGAGGCA-AAGGAGTA_L001_R1_001.fastq.gz',
         'AC738_B1_AAGAGGCA-AAGGAGTA_L001_R2_001.fastq.gz',
         'AC739_CON1_GTAGAGGA-CTAAGCCT_L001_R1_001.fastq.gz',
         'AC739_CON1_GTAGAGGA-CTAAGCCT_L001_R2_001.fastq.gz')

### Add the filename2
sra_df$filename2 <- file2
### filetype
sra_df$filetype <- 'fastq'
### design_description
sra_df <- mutate(sra_df, design_description = ifelse(filename %in% wgs, 
                                                       "Illumina MiSeq 150bp paired-end whole genome sequencing", 
                                                       "Illumina MiSeq 150bp paired-end DNA targeted sequencing"))

### Instrument_model
sra_df$Instrument_model <- 'Illumina MiSeq'
### platform
sra_df$platform <- 'ILLUMINA'
### library_layout 
sra_df$library_layout <- 'paired'
### library_selection
sra_df <- mutate(sra_df, library_selection = ifelse(filename %in% wgs, 
                                                     "RANDOM", 
                                                     "PCR"))
### Add the bioproject_accession
sra_df$bioproject_accession <- 'PRJNA526753'
### library_ID
# Remove end of strings
interum <- substr(sra_df$filename,1,nchar(as.character(sra_df$filename))-21)                        

# Take the last characters
##########################################
substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
}
#########################################
# Collect the lib ids
lib_id <- substrRight(interum, 17)
# Add lib ids
sra_df$library_ID <- lib_id

### library_source
sra_df$library_source <- 'GENOMIC'

### library_strategy
sra_df <- mutate(sra_df, library_strategy = ifelse(filename %in% wgs, 
                                                    "WGS", 
                                                    "AMPLICON"))

### raw_id
raw_id <- substr(sra_df$filename,1,nchar(as.character(sra_df$filename))-39)
sra_df$raw_id <- raw_id

# Load in the hash table
hash_table <- read.csv(file = 'sample-id_file-id.txt', header = FALSE, sep = '\t', comment.char = '#')
names(hash_table) <- c('sample_id', 'raw_id')

### sample_id
sra_df <- dplyr::left_join(sra_df, hash_table, by = "raw_id")

### title
sra_df <- mutate(sra_df, title = ifelse(filename %in% wgs, 
                                                   paste0("WGS of MDV: ", sra_df$sample_id), 
                                                   paste0("AMPLICON DNA-Seq: ", sra_df$sample_id)))


# Load the BioSample labels into a dataframe
df <- read.csv(file = 'BioSampleObjects_mdv.txt', header = TRUE, sep = '\t') 
df <- as_tibble(df)
df <- df[1:6]
names(df) <- c('biosample_accession', 'sample_id', 'spuid', 'organism', 'tax_id','isolate')
names(df)
df <- df %>% dplyr::select(biosample_accession, sample_id)

### biosample_accession
sra_df <- dplyr::left_join(sra_df, df, by = "sample_id")

### Reorder the columns accordingly
names(sra_df)

sra_df <- sra_df %>% dplyr::select(bioproject_accession, biosample_accession, library_ID, title, library_strategy,
                         library_source, library_selection, library_layout, platform, Instrument_model,
                         design_description, filetype, filename, filename2)

### Save the file
write.table(sra_df, file = 'SRA_metadata_mdv.txt', 
            row.names = FALSE, sep = '\t', quote = FALSE)


