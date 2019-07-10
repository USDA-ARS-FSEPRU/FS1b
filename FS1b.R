#R scripts for analyzing and visualizing data from the research paper "Shifts in the nasal microbiota of swine in response to different routes of oxytetracycline administration"

#################################################################################################################################################################################################

#SECOND SECTION (after processing sequences in mothur): Creating OTU Table

#Purpose: Create OTU table with R using specific files generated from mothur
  
#Files needed:
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#FS1bfinal.outsingletons.abund.taxonomy.csv
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared
#FS1bfinal.singleton.abund.subsample.shared.csv

#The output files from mothur needed in R for this section and subsequent sections, aside from fasta files, are text files that can be saved as csv for ease of use in R.

#Load library package
library(tidyverse)

#Set your working directory with the command 
setwd("<directory path>")

#If you have to stop midway, you can save R objects using the function
save.image(file="<yourfilename>.RData")

#When you return, you can use the following function to reload datasets written with the save.image or save functions
load("<yourfilename>.RData")

#To start creating OTU table, edit taxonomy file
#Save "FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy" file generated from mothur as a csv file in a spreadsheet editor. Named file as "FS1bfinal.outsingletons.abund.taxonomy.csv"
taxonomy <- read.csv("FS1bfinal.outsingletons.abund.taxonomy.csv") #Import this csv file from working directory using "read.csv" function
taxonomy$Taxonomy <- gsub('*\\(.*?\\) *', '', taxonomy$Taxonomy) #Substitute (100) and variations of that with nothing ('') from the Taxonomy column
taxonomy[1:6,3] #Show column number 3, rows 1 through 6 in 'taxonomy' dataframe
write.csv(taxonomy, file = "FS1babundsingleton2000taxonomy.csv") #Write 'taxonomy' into a csv vile and open the csv file in a spreadsheet editor to remove the size and numbered rows

#Edit subsample.shared file
#Save "FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared" file generated from mothur as a csv file in a spreadsheet editor. Named file as "FS1bfinal.singleton.abund.subsample.shared.csv"
shared <- read.csv("FS1bfinal.singleton.abund.subsample.shared.csv", stringsAsFactors = FALSE)
head(shared) #Check on the first part of 'shared' dataframe
shared[1:6,1]
shared <- t(shared) #Transpose 'shared'
head(shared)
shared[1:6,1]
write.csv(shared, file = 'FS1babundsingletonshared2000.csv') #Open this csv file in a spreadsheet editor and remove the "V*", "label", and "numOtus" rows

#Read the revised FS1babundsingletonshared2000.csv file
shared <- read.csv("FS1babundsingletonshared2000.csv")
colnames(shared) [1] <- "OTU" #Rename first column of 'shared' to "OTU"
taxonomy <- read.csv("FS1babundsingleton2000taxonomy.csv")
OTUtable <- merge(shared, taxonomy, by.x ="OTU", by.y = "OTU") #Merge 'shared' and 'taxonomy' objects by OTU
head(OTUtable)
nrow(OTUtable) #Count number of rows in 'OTUtable'
ncol(OTUtable) #Count number of columns in 'OTUtable'
write.csv(OTUtable, file= "FS1babundsingleton2000OTUtable.csv") 
#Open this csv file in a spreadsheet editor and remove the first row (numbered rows) and "Size" column, and rename "Taxonomy" header to "taxonomy"

#Check OTU table
OTUtable <-read.csv("FS1babundsingleton2000OTUtable.csv", stringsAsFactors = FALSE)
head(OTUtable) #Check the first part of 'OTUtable' to make sure it looks ok

#################################################################################################################################################################################################

#THIRD SECTION: Creating phyloseq objects for each tissue

#Purpose: Create phyloseq objects to be used to calculate alpha and beta diversity measures for nasal and tonsil tissue samples. This section will also use the adonis function to determine the effect of time and treatment on the community structure of nasal and tonsil microbiota.

#Files needed:
#FS1babundsingleton2000OTUtable.csv
#FS1babundsingleton2000metadata.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(vegdist)

#Read files
otu <- read.csv("FS1babundsingleton2000OTUtable.csv", row.names=1) #Set column 1 as row names
meta <-read.csv("FS1babundsingleton2000metadata.csv")
dim(otu) #Check dimensions of 'otu'
head(otu[,1:5]) #Check the first part of 'otu' table
head(otu[,240:243]) #Check last part of 'otu' table
dim(meta) #Determine the dimensions of 'meta' dataframe. It has 242 rows and 5 columns
head(meta[,1:5]) #Check first part of 'meta' table

#Remove taxonomy from 'otu'
tax <- otu[,(242:243)] #Remove column 243 taxonomy and 242 to copy the row names from 'otu' to 'tax'
head(tax)
colnames(tax)[1] <- "delete" #Rename column 1 of 'tax' (formerly column 242) as "delete" which will be deleted later
head(tax)

#Modify 'otu' with only OTU count data
otu <- otu[,-243] #Remove column 243 taxonomy in 'otu' to have only OTU data
head(otu[,240:242]) 
dim(otu) #Dimensions of 'otu' show 1637 rows 242 columns

#Transpose 'otu' to match format of 'meta'
otu.trans <- t(otu) 
#Now rownames in 'otu.trans' are sample names, columns are OTUs
head(otu.trans[,1:5])
head(meta) #Row names are numbered, but we want sample names as row names
class(meta) #The type of class that 'meta' is is dataframe
class(otu) #dataframe

#Merge 'otu' and 'meta' data frames
otu.meta <- merge(meta, otu.trans, by.x=1, by.y=0) #Merge by names of the columns that are common to both x and y (columns with common names between the two data sets)
#by.x=1 means match by 1st column from 'meta'; y=0 means match by rownames in 'otu.trans'
head(otu.meta[,1:10])
class(otu.meta) #Check class type of 'otu.meta'. It should be a dataframe.

#Added an "All" column (combines 'Treatment' and 'Day' values) in new 'otu.meta2' dataframe
otu.meta2<- cbind(otu.meta) #Make second copy of 'otu.meta' to use to include an "All" column
otu.meta2$All <- with(otu.meta2, paste0(Day, sep=" ", Treatment)) #Combine "Day" and "Treatment" columns into an "All" column
head(otu.meta2) #Check first part of otu.meta2
dim(otu.meta2) #Check dimensions of 'otu.meta2' dataframe
head(otu.meta2[,1640:1643]) #Check the first part of the end of 'otu.meta2' dataframe
otu.meta2<- otu.meta2[,c(1:5,1643,6:1642)] #Reorder columns to have "All" column after "Treatment" column
head(otu.meta2[,1:10]) #Check the first part of the beginning of 'otu.meta2' dataframe
head(otu.meta2[,1640:1643])
write.csv(otu.meta2, file="FS1babundsingleton2000.otu.meta_all.csv")

#Subset nasal and tonsil samples
nw2 <- subset(otu.meta2, Tissue=="NW") #Return subset of 'otu.meta2' dataframe as 'nw2' where row values within 'otu.meta2' are "NW" in "Tissue" column
head(nw2[,1:10])
dim(nw2)
write.csv(nw2, file="FS1babundsingleton2000.otu.meta_all.nasalonly.csv")

tt2 <- subset(otu.meta2, Tissue=="TT")
head(tt2[,1:10])
tail(tt2[,1:10])
dim(tt2)
write.csv(tt2, file="FS1babundsingleton2000.otu.meta_all.tonsilonly.csv")


#Creating phyloseq objects for tonsil samples

#Pull out metadata from 'tt2' dataframe
head(tt2[,1:10])
meta.tt2 <- tt2[,1:6] #Take columns 1-6 ("Sample" to "All") to make 'meta.tt2'
head(meta.tt2)
row.names(meta.tt2) <- meta.tt2[,1] #Make column 1 be rownames for 'meta.tt2'
head(meta.tt2)
meta.tt2 <- meta.tt2[,-1] #Remove the extra "Sample" column
head(meta.tt2)
dim(meta.tt2)

#Create SAM metadata table phyloseq object
SAMtt2 = sample_data(meta.tt2, errorIfNULL = TRUE)
head(SAMtt2)
dim(SAMtt2)

#Pull out OTU data from 'tt2' dataframe
head(tt2[,1:10])
dim(tt2)
row.names(tt2) <- tt2[,1] #Make column 1 be rownames for 'tt2'
head(tt2[,1:10])
tt2 <- tt2[,-1] #Remove the extra "Sample" column
head(tt2[,1:10])
dim(tt2)
head(tt2[,1640:1642])
otu.tt2 <- tt2[,c(6:1642)] #Select "Sample" and "OTU" columns to create 'otu.tt2' dataframe
head(otu.tt2[,1:10])
dim(otu.tt2)
otu.tt2.trans <- t(otu.tt2) #Transpose 'otu.tt2' to have OTUs as rownames, sample names as column names
head(otu.tt2.trans[,1:10])
dim(otu.tt2.trans)

#Merge 'tax' back into 'otu.tt2.trans' for correct format and taxons
head(tax)
otu.tax.tt2 <- merge(otu.tt2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu.tax.tt2) #1637 38
head(otu.tax.tt2[,35:38])
head(otu.tax.tt2[,1:5])
row.names(otu.tax.tt2) <- otu.tax.tt2[,1] #Set first row of 'otu.tax.tt2' as rownames
head(otu.tax.tt2[,1:5])
otu.tax.tt2 <- otu.tax.tt2[,-1] #Remove first row aka extraneous "OTU" column from 'otu.tax.tt2'
head(otu.tax.tt2[,1:5])

#Split 'otu.tax.tt2' again
dim(otu.tax.tt2) 
head(otu.tax.tt2[,35:37])
otu.notax.tt2 <- otu.tax.tt2[,1:35] #Take rows 1-35 to make new dataframe 'otu.notax.tt2' (36 is "delete" column, 37 is "taxonomy" column)
head(otu.notax.tt2[,1:5])
head(otu.notax.tt2[,30:35])
dim(otu.notax.tt2)
class(otu.notax.tt2)
otu.notax.tt2 <- as.matrix(otu.notax.tt2) #Turn 'otu.notax.tt2' into a matrix class
class(otu.notax.tt2) #'otu.notax.tt2' is indeed a matrix

#Create OTU table phyloseq object
OTUtt2 = otu_table(otu.notax.tt2, taxa_are_rows = TRUE)
head(OTUtt2)
dim(OTUtt2)
head(otu.tax.tt2[,30:37])
tax.levels.tt2 <- separate(data = otu.tax.tt2, 
                       col = Taxonomy, 
                       into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#"separate" function separates "Taxonomy" column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax.levels.tt2) #Notice that "Species" column is blank
dim(tax.levels.tt2) #1637 43
tax.only.tt2 <- tax.levels.tt2[,37:42] #Keep only taxonomy columns "Kingdom" up to "Genus"
head(tax.only.tt2)
dim(tax.only.tt2) #1637 6
class(tax.only.tt2) #dataframe
tax.m.tt2 <- as.matrix(tax.only.tt2) #Convert class of 'tax.only.tt2' from dataframe to a matrix
class(tax.m.tt2) #matrix

#Create TAX taxonomy table phyloseq object
TAXtt2 = tax_table(tax.m.tt2)
head(TAXtt2)
dim(TAXtt2) #1637 6

#Create phyloseq object containing taxonomy, metadata, and OTU table
phyloseq.tt2 <- phyloseq(OTUtt2, SAMtt2, TAXtt2)
phyloseq.tt2 #View phyloseq object
#Output:
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1637 taxa and 35 samples ]
#sample_data() Sample Data:       [ 35 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1637 taxa by 6 taxonomic ranks ]

save(phyloseq.tt2, file="phyloseq.tt2.RData")


#Creating phyloseq objects for nasal samples

#Pull out metadata from 'nw2' dataframe
head(nw2[,1:10])
dim(nw2) #207 1643
meta.nw2 <- nw2[,1:6] #Take columns 1-5 to make 'meta.nw2'
head(meta.nw2)
row.names(meta.nw2) <- meta.nw2[,1] #Make column 1 be rownames for 'meta.nw2'
head(meta.nw2)
meta.nw2 <- meta.nw2[,-1] #Remove the extra "Sample" column
head(meta.nw2)
dim(meta.nw2) #207 5

#Create SAM metadata table phyloseq object
SAMnw2 = sample_data(meta.nw2, errorIfNULL = TRUE)
head(SAMnw2)
dim(SAMnw2) #207 5

#Pull out OTU data from 'nw2' dataframe
head(nw2[,1:10])
dim(nw2) #207 1643
head(nw2[,1640:1643])
otu.nw2 <- nw2[,c(1,7:1643)] #Select "Sample" and "OTU" columns to create 'otu.nw2' dataframe
head(otu.nw2[,1:10])
dim(otu.nw2) #207 1638
row.names(otu.nw2) <- otu.nw2[,1] #Make first column be rownames for 'otu.nw2'
head(otu.nw2[,1:10])
otu.nw2 <- otu.nw2[,-1] #Remove the first column
head(otu.nw2[,1:10])
otu.nw2.trans <- t(otu.nw2) #Transpose 'otu.nw2' to have OTUs as rownames, sample names as column names
head(otu.nw2.trans[,1:10])
dim(otu.nw2.trans) #1637 207
head(otu.nw2.trans[,200:207])

#Merge 'tax' back into 'otu.nw2.trans' for correct format and taxons
head(tax)
otu.tax.nw2 <- merge(otu.nw2.trans, tax, by=0) #Merge by rownames aka OTU rownames
dim(otu.tax.nw2) #1637 210
head(otu.tax.nw2[,200:210])
head(otu.tax.nw2[,1:5])
row.names(otu.tax.nw2) <- otu.tax.nw2[,1] #Set first row as rownames
head(otu.tax.nw2[,1:5])
otu.tax.nw2 <- otu.tax.nw2[,-1] #Remove first row, extraneous OTU column
head(otu.tax.nw2[,1:5])
dim(otu.tax.nw2) #1637 209

#Split 'otu.tax.nw2' again
dim(otu.tax.nw2) #1637 209
head(otu.tax.nw2[,205:209])
otu.notax.nw2 <- otu.tax.nw2[,1:207] #Take rows 1-207 to make new dataframe otu.notax.nw2 (208 is "delete" column, 209 is "taxonomy" column)
head(otu.notax.nw2[,1:5])
head(otu.notax.nw2[,200:207])
class(otu.notax.nw2) #dataframe
otu.notax.nw2 <- as.matrix(otu.notax.nw2) #Turn 'otu.notax.nw2' into a matrix class
class(otu.notax.nw2) #matrix

#Create OTU table phyloseq object
OTUnw2 = otu_table(otu.notax.nw2, taxa_are_rows = TRUE)
head(OTUnw2)
dim(OTUnw2) #1637 207
head(otu.tax.nw2[,200:209])
tax.levels.nw2 <- separate(data = otu.tax.nw2, 
                        col = Taxonomy, 
                        into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#Separate Taxonomy column into 7 separate columns labeled "Kingdom", "Phylum", "Class", etc.
head(tax.levels.nw2) #Notice that Species column is blank
dim(tax.levels.nw2) #1637 215
tax.only.nw2 <- tax.levels.nw2[,209:214] #Keep only taxonomy columns from "Kingdom" up to "Genus"
head(tax.only.nw2)
dim(tax.only.nw2) #1637 6
class(tax.only.nw2) #data.frame
tax.m.nw2 <- as.matrix(tax.only.nw2)
class(tax.m.nw2) #matrix

#Create TAX taxonomy table phyloseq object
TAXnw2 = tax_table(tax.m.nw2)
head(TAXnw2)
dim(TAXnw2) #1637 6

#Create phyloseq object containing taxonomy, metadata, and OTU table
phyloseq.nw2 <- phyloseq(OTUnw2, SAMnw2, TAXnw2)
phyloseq.nw2
#Output:
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1637 taxa and 207 samples ]
#sample_data() Sample Data:       [ 207 samples by 5 sample variables ]
#tax_table()   Taxonomy Table:    [ 1637 taxa by 6 taxonomic ranks ]

save(phyloseq.nw2, file="phyloseq.nw2.RData")


#Run adonis function to determine effect of time and treatment on structure of nasal, tonsil microbiota

#Use vegdist function from vegan package to run distance calculations instead of the distance function (original "distance" function that was used below is no longer available) and use those calculations to run through adonis test
#vegdist requires that phyloseq object's OTU table has OTUs listed in the columns and sample names listed in rows. Also, remove any OTUs with taxa_sums = 0 or non-numeric values. For example, this command can help remove OTUs with taxa_sums = 0: 
#OTU <- prune_taxa(taxa_sums(<yourOTUtable>) > 0, <yourOTUtable>)
#If you create a separate phyloseq object with this specific OTU table setup, you should be able to run the vegdist function without any errors and use the output to run through adonis function
      
#Nasal
dist.nw2 <- distance(phyloseq.nw2, method="bray") #Distance calculation using Bray-Curtis
adonis.nw2 <- as(sample_data(phyloseq.nw2), "data.frame")
set.seed(1) #Use set.seed function when running simulations to ensure all results are reproducible
full.nw2 <- adonis(dist.nw2~Day*Treatment, data=adonis.nw2, permutations=9999)
full.nw2 #Display results
#Output:
#Call:
#adonis(formula = dist.nw2 ~ Day * Treatment, data = adonis.nw2,      permutations = 9999) 
    
#Permutation: free
#Number of permutations: 9999
    
#Terms added sequentially (first to last)
    
#               Df  SumsOfSqs   MeanSqs F.Model   R2        Pr(>F)    
# Day             4 14.499      3.6249  22.0060   0.26863   1e-04 ***
# Treatment       2 1.918       0.9591  5.8226    0.03554   1e-04 ***
# Day:Treatment   8 5.932       0.7415  4.5015    0.10990   1e-04 ***
# Residuals     192 31.626      0.1647            0.58594           
# Total         206 53.976                        1.00000           
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    
###Day had largest effect on variation
    
#switched Treatment and Day order, obtained same conclusions: Day had largest effect on variation
full.nw2_2 <- adonis(dist.nw2~Treatment*Day, data=adonis.nw2, permutations=9999)
full.nw2_2
#Output:
#Call:
#  adonis(formula = dist.nw2 ~ Treatment * Day, data = adonis.nw2,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999

#Terms added sequentially (first to last)
  
#                 Df  SumsOfSqs MeanSqs F.Model   R2        Pr(>F)    
#  Treatment       2  1.959     0.9793  5.9449    0.03628   1e-04 ***
#  Day             4  14.459    3.6148  21.9449   0.26788   1e-04 ***
#  Treatment:Day   8  5.932     0.7415  4.5015    0.10990   1e-04 ***
#  Residuals     192  31.626    0.1647            0.58594           
#  Total         206  53.976                      1.00000           
#---

    
#Tonsil
dist.tt2 <- distance(phyloseq.tt2, method="bray")
adonis.tt2 <- as(sample_data(phyloseq.tt2), "data.frame")
set.seed(1)
full.tt2 <- adonis(dist.tt2~Day*Treatment, data=adonis.tt2, permutations=9999)
full.tt2
#Output:
#Call:
#adonis(formula = dist.tt2 ~ Day * Treatment, data = adonis.tt2,      permutations = 9999) 
  
#Permutation: free
#Number of permutations: 9999
  
#Terms added sequentially (first to last)

#               Df  SumsOfSqs MeanSqs   F.Model R2      Pr(>F)
#Day            2   0.2729    0.13644   1.1432  0.06593 0.2711
#Treatment      2   0.2616    0.13082   1.0961  0.06321 0.3268
#Day:Treatment  4   0.5012    0.12530   1.0498  0.12110 0.3719
#Residuals     26   3.1032    0.11935           0.74976       
#Total         34   4.1389                      1.00000 
  
###None have a more larger effect on variation than the other.
    
#switched Treatment and Day order, obtained same conclusions: none had a more larger effect on variation than the others
full.tt2_2 <- adonis(dist.tt2~Treatment*Day, data=adonis.tt2, permutations=9999)
full.tt2_2
#Output:
#Call:
#adonis(formula = dist.tt2 ~ Treatment * Day, data = adonis.tt2,      permutations = 9999) 

#Permutation: free
#Number of permutations: 9999
  
#Terms added sequentially (first to last)
  
#               Df  SumsOfSqs MeanSqs   F.Model R2      Pr(>F)
#Treatment      2   0.2662    0.13310   1.1152  0.06432 0.3049
#Day            2   0.2683    0.13416   1.1240  0.06483 0.2896
#Treatment:Day  4   0.5012    0.12530   1.0498  0.12110 0.3709
#Residuals     26   3.1032    0.11935           0.74976       
#Total         34   4.1389                      1.00000

#################################################################################################################################################################################################

#FOURTH SECTION: Beta diversity
    
#Purpose: This code generates non-metric multidimensional scaling ordination based on Bray-Curtis dissimilarities to create NMDS plots, and runs pairwise.adonis function to identify any significant differences in bacterial composition between treatment groups on a given day for each tissue.
    
#Files needed:
#phyloseq.nw2
#phyloseq.tt2

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
    
#Load these two functions from the funfuns R package (https://github.com/Jtrachsel/funfuns)
NMDS_ellipse
veganCovEllipse
    
#Nasal
    
#Setting up 'phyloseq.nw2' into dataframes for NMDS calculation
nw2.meta <- data.frame(phyloseq.nw2@sam_data) #Make 'phyloseq.nw2 sam_data' into dataframe
nw2.otu <- data.frame(t(phyloseq.nw2@otu_table)) #Make 'phyloseq.nw2 otu_table' into dataframe
class(nw2.meta) #data.frame
rownames(nw2.meta) == row.names(nw2.otu) #Make sure rownames between 'nw2.meta' and 'nw2.otu' match exactly
nw2.meta$numOTUS <- rowSums(nw2.otu > 1) #For rows with sums greater than 1 in 'nw2.otu', move rows and their respective sum values into "numOTUs" column in 'nw2.meta'
head(nw2.meta)
    
#NMDS calculation
nw2.otu[1:10,1:10]
nw_NMDS <- NMDS_ellipse(nw2.meta, nw2.otu, grouping_set = 'All')
#Output:
#Result: [1] "Ordination stress: 0.195744611098944"
    
#Separate meta data and ellipse data to two lists to make NMDS plot
head(nw_NMDS)
nw.metanmds <- nw_NMDS[[1]] #'nw.metanmds' has meta data + MDS calculations. Select this 1st list of 'nw_NMDS' using double brackets
nw.df_ell <- nw_NMDS[[2]] #'nw.df_ell' is accessing 2nd list from 'nw_NMDS' that has ellipse calculations
#Need two separate lists for making NMDS plot
nw.df_ell$group
head(nw.df_ell)
  
#Create "Day" and "Treatment" columns within 'nw.df_ell' for faceting purposes
nw.df_ell$Day <- sub('(D[0-9]+) ([A-Za-z]+)', '\\1', nw.df_ell$group) #Created "Day" column, '\\1' returns the first part of the regular expression (D[0-9]+) from 'nw.df_ell$group'
nw.df_ell$Treatment <- sub('(D[0-9]+) ([A-Za-z]+)', '\\2', nw.df_ell$group) #Created "Treatment" column, '\\2' returns the second part of the sub expression ([A-Za-z]+) from 'nw.df_ell$group'
head(nw.df_ell)
  
#Restructure level order for 'nw.metanmds' and 'nw.df_ell'
nw.metanmds$Day = factor(nw.metanmds$Day, levels = c("D0", "D4", "D7", "D11", "D14"))
nw.df_ell$Day = factor(nw.df_ell$Day, levels = c("D0", "D4", "D7", "D11", "D14"))
levels(nw.df_ell$Day)
    
#Renaming treatment groups Control, Injected and Infeed to NON, IM and IF, respectively, in 'nw.metanmds' and 'nw.df_ell' dataframes
nw.metanmds$Treatment2 = nw.metanmds$Treatment
nw.metanmds$Treatment2 <- as.character(nw.metanmds$Treatment2)
nw.metanmds$Treatment2[nw.metanmds$Treatment2 == 'Control'] <- "NON"
nw.metanmds$Treatment2[nw.metanmds$Treatment2 == 'Injected'] <- "IM"
nw.metanmds$Treatment2[nw.metanmds$Treatment2 == 'Infeed'] <- "IF"
nw.df_ell$Treatment2 = nw.df_ell$Treatment
nw.df_ell$Treatment2 <- as.character(nw.df_ell$Treatment2)
nw.df_ell$Treatment2[nw.df_ell$Treatment2 == 'Control'] <- "NON"
nw.df_ell$Treatment2[nw.df_ell$Treatment2 == 'Injected'] <- "IM"
nw.df_ell$Treatment2[nw.df_ell$Treatment2 == 'Infeed'] <- "IF"
    
#Creating plot from NMDS calculations
nasal.nmdsplot <- ggplot(data=nw.metanmds, aes(x=MDS1, y=MDS2, color=Treatment2)) + geom_point() + 
      geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
      geom_path(data=nw.df_ell, aes(x=NMDS1, y=NMDS2, color=Treatment2, group=group)) + 
      facet_wrap(~Day, scales = 'free') +
      scale_color_brewer(palette="Dark2") +
      theme_gray(base_size = 10) +
      theme(strip.text.x = element_text(size=15), axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
      labs(color="Treatment group")
nasal.nmdsplot2 <- nasal.nmdsplot + scale_colour_manual(values=c("#E69F00", "#56B4E9", "#999999")) + theme(legend.position = "right")
nasal.nmdsplot2
    
#Save 'nasal.nmdsplot2' as a .tiff for publication, 500dpi
ggsave("Nasal_NMDS_DayAndTreatment.tiff", plot=nasal.nmdsplot2, width = 10, height = 6, dpi = 500, units =c("in"))
    
#Using pairwise.adonis function
nw.adon <- pairwise.adonis(nw2.otu, nw2.meta$All) #Run pairwise.adonis on 'nw2.otu' OTU table and "All" column of 'nw2.meta' dataframe
#nw.adon contains all the pairwise comparisons
nw.adon$pairs #List all comparisons in the "pairs" column of 'nw.adon'
goodcomps.nw <- c(grep('D0 [A-Za-z]+ vs D0 [A-Za-z]+', nw.adon$pairs),
                  grep('D4 [A-Za-z]+ vs D4 [A-Za-z]+', nw.adon$pairs),
                  grep('D7 [A-Za-z]+ vs D7 [A-Za-z]+', nw.adon$pairs),
                  grep('D11 [A-Za-z]+ vs D11 [A-Za-z]+', nw.adon$pairs),
                  grep('D14 [A-Za-z]+ vs D14 [A-Za-z]+', nw.adon$pairs))
# "[A-Za-z]" matches all capital and lowercase letters
# "+" matches a whole word and not just one letter (if you didn't have "+", then it would match by one letter)
# "c" creates the vector, lumps all pairs of specific groups of interest together
# You want to make a vector to combine all the pairwise comparisons you're interested in (same day, different treatment group)
nw.adon.good <- nw.adon[goodcomps.nw,] #Rename 'goodcomps.nw' vector to 'nw.adon.good'
nw.adon.good
nw.adon.good$p.adjusted <- p.adjust(nw.adon.good$p.value, method = 'fdr') #"p.adjust" function returns a set of p-values adjusted with "fdr" method
nw.adon.good$p.adjusted2 <- round(nw.adon.good$p.adjusted, 3) #Round p-values to 3 decimal points and list in new "p.adjusted2" column
nw.adon.good$p.adjusted2[nw.adon.good$p.adjusted2 > 0.05] <- NA #For all p-values greater than 0.05, replace with "NA"
write.csv(nw.adon.good, file='nw.adon.good.txt', row.names=TRUE)
    
    
#Tonsil
    
#Setting up 'phyloseq.tt2' into dataframes for NMDS calculation
tt2.meta <- data.frame(phyloseq.tt2@sam_data) #Make phyloseq.tt2 sam_data into dataframe
tt2.otu <- data.frame(t(phyloseq.tt2@otu_table)) #Make phyloseq.tt2 otu_table into dataframe
class(tt2.meta) #data.frame
rownames(tt2.meta) == row.names(tt2.otu) #Make sure rownames between tt2.meta and tt2.otu match exactly. Yes they do.
tt2.meta$numOTUS <- rowSums(tt2.otu > 1)
head(tt2.meta)
    
#NMDS calculation
tt2.otu[1:10,1:10]
tt_NMDS <- NMDS_ellipse(tt2.meta, tt2.otu, grouping_set = 'All')
#Output:
#Result: [1] "Ordination stress: 0.180074413238642"
    
#Separate meta data and ellipse data to two lists to make NMDS plot
head(tt_NMDS)
tt.metanmds <- tt_NMDS[[1]] #'tt.metanmds' has meta data + MDS calculations. Select this 1st list of 'tt_NMDS' using double brackets
class(tt.metanmds) #data.frame

#Count number of samples in each "All" group in 'tt.metanmds' and summarize
tt.metanmds %>% group_by(All) %>% summarise(count=n())
#Output:
## A tibble: 9 x 2
#All          count
#<fct>        <int>
#1 D14 Control      2
#2 D14 Infeed       3
#3 D14 Injected     5
#4 D4 Control       5
#5 D4 Infeed        5
#6 D4 Injected      4
#7 D7 Control       3
#8 D7 Infeed        3
#9 D7 Injected      5

#Remove day 14 samples as there are too few sample points to plot
tt.metanmds2 <- tt.metanmds %>% filter(All != 'D14 Control')
head(tt.metanmds2)
tt.df_ell <- tt_NMDS[[2]] #'tt.df_ell' is accessing 2nd list from 'tt_NMDS' that has ellipse calculations
tt.df_ell2 <- tt.df_ell %>% filter(group != 'D14 Control') 
#Need two separate lists for making NMDS plot
tt.df_ell2$group
#15 levels
head(tt.df_ell2)

#Create "Day" and "Treatment" columns within 'tt.df_ell2' for faceting
tt.df_ell2$Day <- sub('(D[0-9]+) ([A-Za-z]+)', '\\1', tt.df_ell2$group) #'\\1' returns the first part of the regular expression (D[0-9]+) from tt.df_ell$group
tt.df_ell2$Treatment <- sub('(D[0-9]+) ([A-Za-z]+)', '\\2', tt.df_ell2$group) #'\\2' returns the second part of the sub expression ([A-Za-z]+) from tt.df_ell$group
head(tt.df_ell2)
    
#Restructure level order for 'tt.metanmds2' and 'tt.df_ell2'
tt.metanmds2$Day = factor(tt.metanmds2$Day, levels = c("D4", "D7", "D14"))
tt.df_ell2$Day = factor(tt.df_ell2$Day, levels = c("D4", "D7", "D14"))
levels(tt.df_ell2$Day)
    
#Renaming treatment group names to IM and IF in 'tt.metanmds' and 'tt.df_ell' dataframes
tt.metanmds2$Treatment2 = tt.metanmds2$Treatment
tt.metanmds2$Treatment2 <- as.character(tt.metanmds2$Treatment2)
tt.metanmds2$Treatment2[tt.metanmds2$Treatment2 == 'Injected'] <- "IM"
tt.metanmds2$Treatment2[tt.metanmds2$Treatment2 == 'Infeed'] <- "IF"
tt.df_ell2$Treatment2 = tt.df_ell2$Treatment
tt.df_ell2$Treatment2 <- as.character(tt.df_ell2$Treatment2)
tt.df_ell2$Treatment2[tt.df_ell2$Treatment2 == 'Injected'] <- "IM"
tt.df_ell2$Treatment2[tt.df_ell2$Treatment2 == 'Infeed'] <- "IF"
    
#Creating plot from NMDS calculations
ggplot(data=tt.metanmds2, aes(x=MDS1, y=MDS2, color=Treatment2)) + geom_point() + 
      geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY), alpha = 0.5) + 
      geom_path(data=tt.df_ell2, aes(x=NMDS1, y=NMDS2, color=Treatment2, group=group)) + 
      facet_wrap(~Day, scales = 'free') +
      scale_color_brewer(palette="Dark2") +
      theme_gray(base_size = 10) +
      theme(strip.text.x = element_text(size=15)) +
      labs(caption = 'Ordination stress = 0.18',color="Treatment group") 
    
#Using pairwise.adonis function
tt.adon <- pairwise.adonis(tt2.otu, tt2.meta$All) #Run pairwise.adonis on 'tt2.otu' OTU table and "All" column of 'tt2.meta' dataframe
#tt.adon contains all the pairwise comparisons
tt.adon$pairs #List all comparisons in the "pairs" column of 'tt.adon'
goodcomps.tt <- c(grep('D4 [A-Za-z]+ vs D4 [A-Za-z]+', tt.adon$pairs),
                  grep('D7 [A-Za-z]+ vs D7 [A-Za-z]+', tt.adon$pairs),
                  grep('D14 [A-Za-z]+ vs D14 [A-Za-z]+', tt.adon$pairs))
tt.adon.good <- tt.adon[goodcomps.tt,] #Rename 'this'goodcomps.tt' vector to 'tt.adon.good'
tt.adon.good
tt.adon.good$p.adjusted <- p.adjust(tt.adon.good$p.value, method = 'fdr')
tt.adon.good$p.adjusted2 <- round(tt.adon.good$p.adjusted, 3)
tt.adon.good$p.adjusted2[tt.adon.good$p.adjusted2 > 0.05] <- NA
write.csv(tt.adon.good, file='tt.adon.good.txt', row.names=TRUE)

#################################################################################################################################################################################################

#FIFTH SECTION: Alpha diversity
    
#Purpose: This code calculates alpha diversity metrics (Shannon, Inverse Simpson) that are plotted as box and whisker plots, and uses wilcoxon rank sum test to assess any significant differences in diversity between treatment groups on a given day for each tissue.

#Files needed:
#phyloseq.nw2
#phyloseq.tt2

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(scales)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library('wesanderson')
    
#Nasal
    
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
nw2.meta$shannon <- diversity(nw2.otu) #"diversity" is a vegan function. The default index is set at "shannon". I added a shannon index column in 'nw2.meta'
nw2.meta$invsimpson <- diversity(nw2.otu,index = 'invsimpson') #We used 'invsimpson' since it is easier to interpret than Simpson values and won't need to "inverse" the Simpson values to understand (With Simpson values, the lower the number, the higher the diversity)
levels(sample_data(nw2.meta)$Day) # Set the level order of values in "Day" column to "D0"  "D4"  "D7"  "D11" "D14"
    
#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within nw2.meta
nw.average.shannon.invsimpson.numOTUs <- aggregate(nw2.meta[, 5:7], list(nw2.meta$All), mean)
print(nw.average.shannon.invsimpson.numOTUs)
#Output:
#        Group.1  shannon invsimpson   numOTUS
#1    D0 Control 2.487904   6.992742  80.83333
#2     D0 Infeed 2.623665   5.260101  93.21053
#3   D0 Injected 2.493518   4.865172  74.84211
#4   D11 Control 3.567874  12.635943 130.62500
#5    D11 Infeed 1.796345   3.615875  27.71429
#6  D11 Injected 2.833157   8.368972  87.85714
#7   D14 Control 1.974733   3.377601  44.37500
#8    D14 Infeed 2.293317   7.677236  54.14286
#9  D14 Injected 1.890397   2.921545  49.83333
#10   D4 Control 3.574778  15.639098 123.27273
#11    D4 Infeed 3.065280   8.812949 101.40000
#12  D4 Injected 3.455888  12.502168 116.86364
#13   D7 Control 3.520855  13.793185 114.00000
#14    D7 Infeed 2.819091   6.328819  93.07143
#15  D7 Injected 2.822564   6.958273  93.80000
write.csv(nw.average.shannon.invsimpson.numOTUs, file="nw.average.shannon.invsimpson.num.OTUs.txt", row.names=TRUE)
nw2.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(nw2.meta$shannon, nw2.meta$All, p.adjust.method = 'none') #Calculate pairwise comparisons by "All" column of the shannon indices in "Shannon" column
print(nw2.pairwise.wilcox.shannon.test) #Look at the results of 'nw2.pairwise.wilcox.shannon.test'
nw2.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(nw2.meta$invsimpson, nw2.meta$All, p.adjust.method = 'none')
print(nw2.pairwise.wilcox.invsimpson.test)
    
#Change level names of 'nw2.meta' "Day" column from "D#" to "#"
levels(nw2.meta$Day)
levels(nw2.meta$Day) <- c("0", "4", "7", "11", "14")
    
#Renaming treatment group names from Control, Injected, and Infeed to NON, IM, and IF, respectively
nw2.meta$Treatment2 = nw2.meta$Treatment #Duplicate "Treatment" column and name as "Treatment2" column
nw2.meta$Treatment2 <- as.character(nw2.meta$Treatment2) #Covert values within "Treatment2" column to characters
nw2.meta$Treatment2[nw2.meta$Treatment2 == 'Control'] <- "NON" #Rename all "NON" characters in "Treatment2" column as "Control"
nw2.meta$Treatment2[nw2.meta$Treatment2 == 'Injected'] <- "IM"
nw2.meta$Treatment2[nw2.meta$Treatment2 == 'Infeed'] <- "IF"
    
#Generate a box and whisker plot of shannon (both shannon and inverse simpson diversity indices showed same trends; I chose to make a plot using shannon indices)
(nw2.shan <- ggplot(data = nw2.meta, aes(x=Treatment2, y=shannon, group=All, fill=Treatment2)) +
        geom_boxplot(position = position_dodge2(preserve = 'total')) +
        facet_wrap(~Day, scales = 'free') +
        scale_y_continuous(name = "Shannon diversity") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              strip.text.x = element_text(size=14),
              axis.title.y = element_text(size=15)) +
        theme(legend.position = "none"))
# "free" within "facet_wrap" allows each plot to customize the scale to the specific data set (no forced scaling applied to all plots)
# "position = position_dodge2(preserve = 'total')" fixes the ggplot box width, making them wider, prevents narrow boxes from forming in the plot
nw2.shan <- nw2.shan + scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999"))
nw2.shan
    
#Save 'nw2.shan' as a .tiff for publication, 500dpi
ggsave("Nasal_Shannon.tiff", plot=nw2.shan, width = 10, height = 5, dpi = 500, units =c("in"))
    
    
#Tonsil
    
#Calculating alpha diversity metrics: Shannon, Inverse Simpson
tt2.meta$shannon <- diversity(tt2.otu)
tt2.meta$invsimpson <- diversity(tt2.otu,index = 'invsimpson')
levels(sample_data(tt2.meta)$Day) #Levels are: "D4" "D7" "D14"
    
#Calculate the average shannon, invsimpson, numOTUs for each "All" subtype within 'tt2.meta'
tt.average.shannon.invsimpson.numOTUs <- aggregate(tt2.meta[, 6:8], list(tt2.meta$All), mean)
print(tt.average.shannon.invsimpson.numOTUs)
#Output:
#       Group.1  numOTUS  shannon invsimpson
#1  D14 Control 53.50000 2.900513  10.421222
#2   D14 Infeed 47.66667 2.840736  10.343016
#3 D14 Injected 51.00000 2.881732  10.383069
#4   D4 Control 60.40000 2.934585  10.827341
#5    D4 Infeed 52.80000 2.745635   9.454290
#6  D4 Injected 47.00000 2.432809   6.287067
#7   D7 Control 54.66667 2.776998  11.085440
#8    D7 Infeed 54.33333 3.031819  11.932736
#9  D7 Injected 45.60000 2.684307   7.741800
write.csv(tt.average.shannon.invsimpson.numOTUs, file="tt.average.shannon.invsimpson.num.OTUs.txt", row.names=TRUE)
tt2.pairwise.wilcox.shannon.test <- pairwise.wilcox.test(tt2.meta$shannon, tt2.meta$All, p.adjust.method = 'none')
print(tt2.pairwise.wilcox.shannon.test)
tt2.pairwise.wilcox.invsimpson.test <- pairwise.wilcox.test(tt2.meta$invsimpson, tt2.meta$All, p.adjust.method = 'none')
print(tt2.pairwise.wilcox.invsimpson.test)

#################################################################################################################################################################################################

#SIXTH SECTION: Magnitude of Change in Nasal Microbiota
    
#Purpose: This code plots the F-statistic from PERMANOVA pairwise comparisons of NON and either IM or IF groups, over time (displays the magnitude of change in the nasal bacterial community structure of the two antibiotic treatment groups relative to control)
    
#Files needed:
#Nasal_MagnitudeOfChange.csv

#Load library packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(philentropy)
library(cowplot)
library(ggrepel)
    
#The F-values were obtained from nw.adon.good.txt data generated in the "Beta diversity" section
#Nasal_MagnitudeOfChange.csv was created from nw.adon.good.txt data that was rearranged
#To make the Nasal_MagnitudeOfChange.csv file, open nw.adon.good.txt file (created from FS1b_beta_diversity_each_tissue_phyloseq.R) in excel, copy columns "F. Model" through "p.adjusted2" and paste in a separate spreadsheet. Add "Day" and "Treatment" columns and save as "Nasal_MagnitudeOfChange.csv".
    
nasal <- read.csv("Nasal_MagnitudeOfChange.csv")
class(nasal)
nasal$Day <- factor(nasal$Day) #Encode "Day" as a factor
nasal2 <- ggplot(data=nasal, aes(x=Day, y=F.Model, group=Treatment)) +
      geom_line(aes(color=Treatment)) + geom_point(aes(color=Treatment)) +
      ylab("PERMANOVA F vs control \n(difference relative to control)") +
      scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
      scale_color_manual(values=c("#E69F00", "#56B4E9")) +  
      geom_label_repel(aes(label=p.adjusted2), box.padding = 0.35, point.padding=0.5,segment.color = 'grey50') +
      theme_classic(base_size = 12) +
      theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.text=element_text(size=14), legend.title=element_text(size=14)) +
      labs(color="Treatment group")
nasal2
    
#Save 'nasal2' as a .tiff for publication, 500dpi
ggsave("Nasal_Magnitude.tiff", plot=nasal2, width = 10, height = 5, dpi = 500, units =c("in"))
    
#################################################################################################################################################################################################

#SEVENTH SECTION: Differential Abundance of Genera in Nasal Microbiota using DESeq2
    
#Purpose: This code uses DESeq2 package to identify nasal microbial genera that were differentially abundant between treatment groups and control group
    
#Files needed:
#Mothur shared file: FS1bfinal.outsingletons.abund.opti_mcc.shared
#Mothur constaxonomy file: FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#Metadata: FS1babundsingleton2000metadata.csv
#FS1b_FinalDiffAbundNasalGenus_IC.csv
#FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv

#Load library packages
library(DESeq2)
library(phyloseq)
library(ggplot2)
library(tidyr)
library("wesanderson")
library(plotly)
library(gapminder)
library(cowplot)
    
#Annotations
#fc = infeed, control
#ic = injected, control
#fi = infeed, injected
#if = injected, infeed
    
#Preparing objects for DESeq2: load files
otu2 <- import_mothur(mothur_shared_file = 'FS1bfinal.outsingletons.abund.opti_mcc.shared')
taxo2 <- import_mothur(mothur_constaxonomy_file = 'FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
taxo2
meta2 <- read.table(file = 'FS1babundsingleton2000metadata.csv', sep = ',', header = TRUE)
    
#Organize 'meta2' meta file
rownames(meta2) <- meta2$Sample #Make names in "Sample" become rownames for 'meta2'
meta2$Day <- gsub("D", "", meta2$Day) # remove "D"
meta2$Treatment2 = meta2$Treatment
meta2$Treatment2 <- as.character(meta2$Treatment2)
meta2$Treatment2[meta2$Treatment2 == 'Control'] <- "NON"
meta2$Treatment2[meta2$Treatment2 == 'Injected'] <- "IM"
meta2$Treatment2[meta2$Treatment2 == 'Infeed'] <- "IF"
meta2$Set <- paste(meta2$Day, meta2$Tissue, meta2$Treatment2, sep = '_')
    
#Make phyloseq object 'FS1b' (combines taxonomy, OTU, and metadata)
phy_meta2 <- sample_data(meta2) 
FS1b <- phyloseq(otu2, taxo2)
FS1b <- merge_phyloseq(FS1b, phy_meta2)   #Combines the 'phy_meta2' metadata with 'FS1b' phyloseq object
colnames(tax_table(FS1b)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS1b
    
#Prune samples from 'FS1b'
FS1b <- prune_samples(sample_sums(FS1b) > 2000, FS1b)  #This removes samples that have fewer than 2000 sequences associated with them.
FS1b <- prune_taxa(taxa_sums(FS1b) > 10, FS1b)        #This removes OTUs that occur less than 10 times globally
tax_table(FS1b) [1:5, 1:6] #Checking what's in 'tax_table' first 5 rows, first 6 columns
    
#Grouping OTUs by Genus using the tax_glom function 
FS1b.genus <- tax_glom(FS1b, taxrank = "Genus")
#This method merges species that have the same taxonomy at a certain taxanomic rank (in this case, by "Genus"). 
#Its approach is analogous to tip_glom, but uses categorical data instead of a tree. 
    
######################################## PRIMARY COMPARISONS TO MAKE ####################################################
    
# NMDS plot showed that disperion is different between days, so I subsetted by day and tissue
    
# Important comparisons to make (significant changes in beta diversity between treatments)
    
# Compare Days 4, 7, 11 NON compared to each of the two treatments 
# Compare Days 4, 7, 11 IM and IF
# Compare Day 14 NON and IF
    
# Other comparisons to make (no significant changes in beta diversity between treatments)
    
# Compare Days 0, 14 IM and IF
# Compare Day 0 NON compared to each of the two treatments
# Compare Day 14 NON and IM
    
##################################################### Day 0 ############################################################
    
############## Day 0 Nasal #########################
    
sample_data(FS1b.genus)
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D0.nw.De'
FS1b.D0.nw <- subset_samples(FS1b.genus, Day == 0 & Tissue == 'NW')
sample_sums(FS1b.D0.nw) #Returns the total number of individuals observed from each sample
FS1b.D0.nw <- prune_taxa(taxa_sums(FS1b.D0.nw) > 1, FS1b.D0.nw) #If taxa_sums is >1, then it will print that out in 'FS1b.D0.nw' object and not include any samples with taxa of sums <1.
rowSums(FS1b.D0.nw@otu_table)
FS1b.D0.nw.De <- phyloseq_to_deseq2(FS1b.D0.nw, ~ Set) # "~Set" : define the major sample covariate as the study design factor. This will be whatever you want to group data by, whatever column you used to designate ellipses with
FS1b.D0.nw.De <- DESeq(FS1b.D0.nw.De, test = "Wald", fitType = "parametric") #Differential expression analysis based on the negative binomial (aka Gamma-Poisson) distribution
    
######### Day 0 Nasal IF vs NON ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "0_NW_IF")
#Output:
#IF = 19
sum(meta2$Set == "0_NW_NON")
#Output:
#NON = 18
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.fc' table
res.D0.fc = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IF","0_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D0.fc
sigtab.D0.fc = res.D0.fc[which(res.D0.fc$padj < .05), ]
sigtab.D0.fc = cbind(as(sigtab.D0.fc, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.fc), ], "matrix")) #"cbind" combines all columns together, regardless of rownames (if you want to match by rownames, use merge function)
format(sigtab.D0.fc$padj, scientific = TRUE)
sigtab.D0.fc$newp <- format(round(sigtab.D0.fc$padj, digits = 3), scientific = TRUE)
sigtab.D0.fc$Treatment <- ifelse(sigtab.D0.fc$log2FoldChange >=0, "IF", "NON") #Assigning "IF" = yes, "NON" = no; it's important to make sure you have the correct group names in the "yes" and "no" position for "ifelse" function
sigtab.D0.fc
    
#Summarize 'sigtab.D0.fc' table
sum.sigtab.D0.fc <- summary(sigtab.D0.fc)
sum.sigtab.D0.fc
    
#plot 'sigtab.D0.fc'
deseq.D0.fc <- ggplot(sigtab.D0.fc, aes(x=reorder(rownames(sigtab.D0.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 0')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
# "reorder" makes the logfold changes appear in numerical order (largest logfold value at the top and at the bottom of the plot), making it easier to see
deseq.D0.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.fc'
sigtab.D0.fc$OTU <- rownames(sigtab.D0.fc)
sigtab.D0.fc$comp <- 'D0_nasal_IFvsNON'
    
#Create a final table ('final.nonsigtab') that lists all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D0.fc'. 
#Other within-day comparisons that had no significant changes in beta diversity between two treatment groups will be added to 'final.nonsigtab'.
final.nonsigtab <- sigtab.D0.fc
    
########## Day 0 Nasal IM vs NON  ####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "0_NW_IM")
#Output:
#IM = 19
sum(meta2$Set == "0_NW_NON")
#Output:
#NON = 18
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.ic' table
FS1b.D0.nw.De$Set
res.D0.ic = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IM","0_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D0.ic = res.D0.ic[which(res.D0.ic$padj < .05), ]
sigtab.D0.ic = cbind(as(sigtab.D0.ic, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.ic), ], "matrix"))
format(sigtab.D0.ic$padj, scientific = TRUE)
sigtab.D0.ic$newp <- format(round(sigtab.D0.ic$padj, digits = 3), scientific = TRUE)
sigtab.D0.ic$Treatment <- ifelse(sigtab.D0.ic$log2FoldChange >=0, "IM", "NON")
    
#Summarize 'sigtab.D0.ic' table
sum.sigtab.D0.ic <- summary(sigtab.D0.ic)
sum.sigtab.D0.ic
    
#Plot 'sigtab.D0.ic'
deseq.D0.ic <- ggplot(sigtab.D0.ic, aes(x=reorder(rownames(sigtab.D0.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 0')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4E9', '#999999'))
deseq.D0.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.ic'
sigtab.D0.ic$OTU <- rownames(sigtab.D0.ic)
sigtab.D0.ic
sigtab.D0.ic$comp <- 'D0_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D0.ic' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.ic)
    
######### Day 0 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "0_NW_IF")
#Output:
#IF = 19
sum(meta2$Set == "0_NW_IM")
#Output:
#IM = 19
    
#Extract results from 'FS1b.D0.nw.De' DESeq object and organize into 'sigtab.D0.if' table
res.D0.if = results(FS1b.D0.nw.De, contrast=c("Set","0_NW_IM","0_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D0.if
sigtab.D0.if = res.D0.if[which(res.D0.if$padj < .05), ]
sigtab.D0.if = cbind(as(sigtab.D0.if, "data.frame"), as(tax_table(FS1b.D0.nw)[rownames(sigtab.D0.if), ], "matrix"))
format(sigtab.D0.if$padj, scientific = TRUE)
sigtab.D0.if$newp <- format(round(sigtab.D0.if$padj, digits = 3), scientific = TRUE)
sigtab.D0.if$Treatment <- ifelse(sigtab.D0.if$log2FoldChange >=0, "IM", "IF")
sigtab.D0.if
    
#Summarize 'sigtab.D0.if' table
sum.sigtab.D0.if <- summary(sigtab.D0.if)
sum.sigtab.D0.if
    
#Plot 'sigtab.D0.if'
deseq.D0.if <- ggplot(sigtab.D0.if, aes(x=reorder(rownames(sigtab.D0.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 0')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust = 0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "IM"), values = c('#E69F00', '#56B4E9'))
deseq.D0.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D0.if'
sigtab.D0.if$OTU <- rownames(sigtab.D0.if)
sigtab.D0.if$comp <- 'D0_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D0.if' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D0.if)
    
    
    
##################################################### Day 4 ############################################################
    
############## Day 4 Nasal #########################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D4.nw.De'
FS1b.D4.nw <- subset_samples(FS1b.genus, Day == 4 & Tissue == 'NW')
sample_sums(FS1b.D4.nw)
FS1b.D4.nw <- prune_taxa(taxa_sums(FS1b.D4.nw) > 1, FS1b.D4.nw)
rowSums(FS1b.D4.nw@otu_table)
FS1b.D4.nw.De <- phyloseq_to_deseq2(FS1b.D4.nw, ~ Set)
FS1b.D4.nw.De <- DESeq(FS1b.D4.nw.De, test = "Wald", fitType = "parametric")
    
######### Day 4 Nasal IF vs NON ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "4_NW_IF")
#Output:
#IF = 20
sum(meta2$Set == "4_NW_NON")
#Output:
#NON = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.fc' table
res.D4.fc = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IF","4_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D4.fc
sigtab.D4.fc = res.D4.fc[which(res.D4.fc$padj < .05), ]
sigtab.D4.fc = cbind(as(sigtab.D4.fc, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.fc), ], "matrix"))
format(sigtab.D4.fc$padj, scientific = TRUE)
sigtab.D4.fc$newp <- format(round(sigtab.D4.fc$padj, digits = 3), scientific = TRUE)
sigtab.D4.fc$Treatment <- ifelse(sigtab.D4.fc$log2FoldChange >=0, "IF", "NON")
sigtab.D4.fc
    
#Plot 'sigtab.D4.fc'
deseq.D4.fc <- ggplot(sigtab.D4.fc, aes(x=reorder(rownames(sigtab.D4.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 4')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', "#999999"))
deseq.D4.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.fc'
sigtab.D4.fc$OTU <- rownames(sigtab.D4.fc)
sigtab.D4.fc$comp <- 'D4_nasal_IFvsNON'
#If there are duplicate OTU rownames, it'll add an extra "1, 2, 3" etc. in numerical order at the end of the rowname
    
#Create a final table ('final.sigtab') that lists all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D4.fc'. 
#Other within-day comparisons that had significant changes in beta diversity between two treatment groups will be added to 'final.sigtab'.
final.sigtab <- sigtab.D4.fc
    
########## Day 4 Nasal IM vs NON  ####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "4_NW_IM")
#Output:
#IM = 22
sum(meta2$Set == "4_NW_NON")
#Output:
#NON = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.ic' table
FS1b.D4.nw.De
res.D4.ic = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IM","4_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D4.ic = res.D4.ic[which(res.D4.ic$padj < .05), ]
sigtab.D4.ic = cbind(as(sigtab.D4.ic, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.ic), ], "matrix"))
format(sigtab.D4.ic$padj, scientific = TRUE)
sigtab.D4.ic$newp <- format(round(sigtab.D4.ic$padj, digits = 3), scientific = TRUE)
sigtab.D4.ic$Treatment <- ifelse(sigtab.D4.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D4.ic'
deseq.D4.ic <- ggplot(sigtab.D4.ic, aes(x=reorder(rownames(sigtab.D4.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 4')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM","NON"), values = c('#56B4E9', '#999999'))
deseq.D4.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.ic'
sigtab.D4.ic
sigtab.D4.ic$OTU <- rownames(sigtab.D4.ic)
sigtab.D4.ic
sigtab.D4.ic$comp <- 'D4_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D4.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D4.ic)
    
######### Day 4 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "4_NW_IF")
#Output:
#IF = 20
sum(meta2$Set == "4_NW_IM")
#Output:
#IM = 22
    
#Extract results from 'FS1b.D4.nw.De' DESeq object and organize into 'sigtab.D4.if' table
res.D4.if = results(FS1b.D4.nw.De, contrast=c("Set","4_NW_IM","4_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D4.if
sigtab.D4.if = res.D4.if[which(res.D4.if$padj < .05), ]
sigtab.D4.if = cbind(as(sigtab.D4.if, "data.frame"), as(tax_table(FS1b.D4.nw)[rownames(sigtab.D4.if), ], "matrix"))
format(sigtab.D4.if$padj, scientific = TRUE)
sigtab.D4.if$newp <- format(round(sigtab.D4.if$padj, digits = 3), scientific = TRUE)
sigtab.D4.if$Treatment <- ifelse(sigtab.D4.if$log2FoldChange >=0, "IM", "IF")
sigtab.D4.if
    
#Plot 'sigtab.D4.if'
deseq.D4.if <- ggplot(sigtab.D4.if, aes(x=reorder(rownames(sigtab.D4.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 4')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF' , 'IM'), values = c('#E69F00', '#56B4EF'))
deseq.D4.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D4.if'
sigtab.D4.if$OTU <- rownames(sigtab.D4.if)
sigtab.D4.if$comp <- 'D4_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D4.if' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D4.if)
    
################################################## Day 7 ############################################################
    
########## Day 7 Nasal #####################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D7.nw.De'
FS1b.D7.nw <- subset_samples(FS1b.genus, Day == 7 & Tissue == 'NW')
sample_sums(FS1b.D7.nw)
FS1b.D7.nw <- prune_taxa(taxa_sums(FS1b.D7.nw) > 1, FS1b.D7.nw)
rowSums(FS1b.D7.nw@otu_table)
FS1b.D7.nw.De <- phyloseq_to_deseq2(FS1b.D7.nw, ~ Set)
FS1b.D7.nw.De <- DESeq(FS1b.D7.nw.De, test = "Wald", fitType = "parametric")
FS1b.D7.nw.De$Set
    
########## Day 7 Nasal IF vs NON #####################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "7_NW_IF")
#Output:
#IF = 14
sum(meta2$Set == "7_NW_NON")
#Output:
#NON = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.fc' table
res.D7.fc = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IF","7_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D7.fc = res.D7.fc[which(res.D7.fc$padj < .05), ]
sigtab.D7.fc = cbind(as(sigtab.D7.fc, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.fc), ], "matrix"))
format(sigtab.D7.fc$padj, scientific = TRUE)
sigtab.D7.fc$newp <- format(round(sigtab.D7.fc$padj, digits = 3), scientific = TRUE)
sigtab.D7.fc$Treatment <- ifelse(sigtab.D7.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D7.fc'
deseq.D7.fc <- ggplot(sigtab.D7.fc, aes(x=reorder(rownames(sigtab.D7.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 7')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF","NON"), values = c('#E69F00','#999999'))
deseq.D7.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.fc'
sigtab.D7.fc$OTU <- rownames(sigtab.D7.fc)
sigtab.D7.fc$comp <- 'D7_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D7.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.fc)
    
############# Day 7 Nasal IM vs NON  ######################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "7_NW_IM")
#Output:
#IM = 15
sum(meta2$Set == "7_NW_NON")
#Output:
#NON = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.ic' table
FS1b.D7.nw.De
res.D7.ic = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_NON")) 
sigtab.D7.ic = res.D7.ic[which(res.D7.ic$padj < .05), ]
sigtab.D7.ic = cbind(as(sigtab.D7.ic, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.ic), ], "matrix"))
format(sigtab.D7.ic$padj, scientific = TRUE)
sigtab.D7.ic$newp <- format(round(sigtab.D7.ic$padj, digits = 3), scientific = TRUE)
sigtab.D7.ic$Treatment <- ifelse(sigtab.D7.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D7.ic'
deseq.D7.ic <- ggplot(sigtab.D7.ic, aes(x=reorder(rownames(sigtab.D7.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 7')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM","NON"), values = c('#56B4EF', '#999999'))
deseq.D7.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.ic'
sigtab.D7.ic$OTU <- rownames(sigtab.D7.ic)
sigtab.D7.ic$comp <- 'D7_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D7.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.ic)
    
######### Day 7 Nasal Wash IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "7_NW_IF")
#Output:
#IF = 14
sum(meta2$Set == "7_NW_IM")
#Output:
#IM = 15
    
#Extract results from 'FS1b.D7.nw.De' DESeq object and organize into 'sigtab.D7.if' table
res.D7.if = results(FS1b.D7.nw.De, contrast=c("Set","7_NW_IM","7_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D7.if
sigtab.D7.if = res.D7.if[which(res.D7.if$padj < .05), ]
sigtab.D7.if = cbind(as(sigtab.D7.if, "data.frame"), as(tax_table(FS1b.D7.nw)[rownames(sigtab.D7.if), ], "matrix"))
format(sigtab.D7.if$padj, scientific = TRUE)
sigtab.D7.if$newp <- format(round(sigtab.D7.if$padj, digits = 3), scientific = TRUE)
sigtab.D7.if$Treatment <- ifelse(sigtab.D7.if$log2FoldChange >=0, "IM", "IF")
sigtab.D7.if
    
#Plot 'sigtab.D7.if'
deseq.D7.if <- ggplot(sigtab.D7.if, aes(x=reorder(rownames(sigtab.D7.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 7')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4EF'))
deseq.D7.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D7.if'
sigtab.D7.if$OTU <- rownames(sigtab.D7.if)
sigtab.D7.if$comp <- 'D7_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D7.if' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D7.if)
    
################################################## Day 11 ############################################################
    
############## Day 11 Nasal ###############
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D11.nw.De'
FS1b.D11.nw <- subset_samples(FS1b.genus, Day == 11 & Tissue == 'NW')
sample_sums(FS1b.D11.nw)
FS1b.D11.nw <- prune_taxa(taxa_sums(FS1b.D11.nw) > 1, FS1b.D11.nw)
rowSums(FS1b.D11.nw@otu_table)
FS1b.D11.nw.De <- phyloseq_to_deseq2(FS1b.D11.nw, ~ Set)
FS1b.D11.nw.De <- DESeq(FS1b.D11.nw.De, test = "Wald", fitType = "parametric")
FS1b.D11.nw.De$Set
    
############## Day 11 Nasal IF vs NON ###############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "11_NW_IF")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.fc' table
res.D11.fc = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IF","11_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D11.fc = res.D11.fc[which(res.D11.fc$padj < .05), ]
sigtab.D11.fc = cbind(as(sigtab.D11.fc, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.fc), ], "matrix"))
format(sigtab.D11.fc$padj, scientific = TRUE)
sigtab.D11.fc$newp <- format(round(sigtab.D11.fc$padj, digits = 3), scientific = TRUE)
sigtab.D11.fc$Treatment <- ifelse(sigtab.D11.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D11.fc'
deseq.D11.fc <- ggplot(sigtab.D11.fc, aes(x=reorder(rownames(sigtab.D11.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 11')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
deseq.D11.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.fc'
sigtab.D11.fc$OTU <- rownames(sigtab.D11.fc)
sigtab.D11.fc$comp <- 'D11_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D11.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D11.fc)
    
########### Day 11 Nasal IM vs NON  ############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "11_NW_IM")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.ic' table
FS1b.D11.nw.De$Set
res.D11.ic = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IM","11_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D11.ic = res.D11.ic[which(res.D11.ic$padj < .05), ]
sigtab.D11.ic = cbind(as(sigtab.D11.ic, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.ic), ], "matrix"))
format(sigtab.D11.ic$padj, scientific = TRUE)
sigtab.D11.ic$newp <- format(round(sigtab.D11.ic$padj, digits = 3), scientific = TRUE)
sigtab.D11.ic$Treatment <- ifelse(sigtab.D11.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D11.ic'
deseq.D11.ic <- ggplot(sigtab.D11.ic, aes(x=reorder(rownames(sigtab.D11.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 11')+ coord_flip() +
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4EF', '#999999'))
deseq.D11.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.ic'
sigtab.D11.ic$OTU <- rownames(sigtab.D11.ic)
sigtab.D11.ic$comp <- 'D11_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D11.ic' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D11.ic)
    
######### Day 11 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "11_NW_IM")
#Output:
#IF = 7
sum(meta2$Set == "11_NW_IF")
#Output:
#IM = 7
    
#Extract results from 'FS1b.D11.nw.De' DESeq object and organize into 'sigtab.D11.if' table
res.D11.if = results(FS1b.D11.nw.De, contrast=c("Set","11_NW_IM","11_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D11.if
sigtab.D11.if = res.D11.if[which(res.D11.if$padj < .05), ]
sigtab.D11.if = cbind(as(sigtab.D11.if, "data.frame"), as(tax_table(FS1b.D11.nw)[rownames(sigtab.D11.if), ], "matrix"))
format(sigtab.D11.if$padj, scientific = TRUE)
sigtab.D11.if$newp <- format(round(sigtab.D11.if$padj, digits = 3), scientific = TRUE)
sigtab.D11.if$Treatment <- ifelse(sigtab.D11.if$log2FoldChange >=0, "IM", "IF")
sigtab.D11.if
    
#Plot 'sigtab.D11.if'
deseq.D11.if <- ggplot(sigtab.D11.if, aes(x=reorder(rownames(sigtab.D11.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF at Nasal Site on Day 11')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4E9'))
deseq.D11.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D11.if'
sigtab.D11.if$OTU <- rownames(sigtab.D11.if)
sigtab.D11.if$comp <- 'D11_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D11.if' to 'final.sigtab'
final.sigtab <- rbind(final.sigtab, sigtab.D11.if)
    
################################################## Day 14 ############################################################
    
########### D14 Nasal IF vs NON################
    
#Convert Phyloseq Data 'FS1b.genus' to DESeq2 object 'FS1b.D14.nw.De'
FS1b.D14.nw <- subset_samples(FS1b.genus, Day == 14 & Tissue == 'NW')
sample_sums(FS1b.D14.nw)
FS1b.D14.nw <- prune_taxa(taxa_sums(FS1b.D14.nw) > 1, FS1b.D14.nw)
rowSums(FS1b.D14.nw@otu_table)
FS1b.D14.nw.De <- phyloseq_to_deseq2(FS1b.D14.nw, ~ Set)
FS1b.D14.nw.De <- DESeq(FS1b.D14.nw.De, test = "Wald", fitType = "parametric")
FS1b.D14.nw.De$Set
    
########### D14 Nasal IF vs NON################
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.fc' table
FS1b.D14.nw.De$Set
res.D14.fc = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IF","14_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D14.fc = res.D14.fc[which(res.D14.fc$padj < .05), ]
sigtab.D14.fc = cbind(as(sigtab.D14.fc, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.fc), ], "matrix"))
sigtab.D14.fc
format(sigtab.D14.fc$padj, scientific = TRUE)
sigtab.D14.fc$newp <- format(round(sigtab.D14.fc$padj, digits = 3), scientific = TRUE)
sigtab.D14.fc$Treatment <- ifelse(sigtab.D14.fc$log2FoldChange >=0, "IF", "NON")
    
#Plot 'sigtab.D14.fc'
deseq.D14.fc <- ggplot(sigtab.D14.fc, aes(x=reorder(rownames(sigtab.D14.fc), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.fc), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IF Relative to NON at Nasal Site on Day 14')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IF", "NON"), values = c('#E69F00', '#999999'))
deseq.D14.fc
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.fc'
sigtab.D14.fc$OTU <- rownames(sigtab.D14.fc)
sigtab.D14.fc$comp <- 'D14_nasal_IFvsNON'
    
#Add all genera that were differentially abundant between IF and NON treatment groups from 'sigtab.D14.fc' to 'final.sigtab' 
final.sigtab <- rbind(final.sigtab, sigtab.D14.fc)
    
#Write 'final.sigtab' into csv file
write.csv(final.sigtab, file= "Nasal_FinalDiffAbundGenus.csv")
    
########### Day 14 Nasal IM vs NON  ############
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe:
sum(meta2$Set == "14_NW_IM")
#Output:
#IF = 6
sum(meta2$Set == "14_NW_NON")
#Output:
#IM = 8
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.ic' table
FS1b.D14.nw.De$Set
res.D14.ic = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IM","14_NW_NON"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
sigtab.D14.ic = res.D14.ic[which(res.D14.ic$padj < .05), ]
sigtab.D14.ic = cbind(as(sigtab.D14.ic, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.ic), ], "matrix"))
format(sigtab.D14.ic$padj, scientific = TRUE)
sigtab.D14.ic$newp <- format(round(sigtab.D14.ic$padj, digits = 3), scientific = TRUE)
sigtab.D14.ic$Treatment <- ifelse(sigtab.D14.ic$log2FoldChange >=0, "IM", "NON")
    
#Plot 'sigtab.D14.ic'
deseq.D14.ic <- ggplot(sigtab.D14.ic, aes(x=reorder(rownames(sigtab.D14.ic), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ic), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to NON at Nasal Site on Day 14')+ coord_flip() +
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c("IM", "NON"), values = c('#56B4E9', '#999999'))
deseq.D14.ic
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.ic'
sigtab.D14.ic$OTU <- rownames(sigtab.D14.ic)
sigtab.D14.ic$comp <- 'D14_nasal_IMvsNON'
    
#Add all genera that were differentially abundant between IM and NON treatment groups from 'sigtab.D14.ic' to 'final.nonsigtab' 
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.ic)
    
######### Day 14 Nasal IM vs IF ###################
    
#Determine number of pigs per group from "Set" column in 'meta2' dataframe: 
sum(meta2$Set == "14_NW_IM")
#Output:
#IF = 6
sum(meta2$Set == "14_NW_IF")
#Output:
#IM = 7
    
#Extract results from 'FS1b.D14.nw.De' DESeq object and organize into 'sigtab.D14.if' table
res.D14.if = results(FS1b.D14.nw.De, contrast=c("Set","14_NW_IM","14_NW_IF"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
res.D14.if
sigtab.D14.if = res.D14.if[which(res.D14.if$padj < .05), ]
sigtab.D14.if = cbind(as(sigtab.D14.if, "data.frame"), as(tax_table(FS1b.D14.nw)[rownames(sigtab.D14.if), ], "matrix"))
format(sigtab.D14.if$padj, scientific = TRUE)
sigtab.D14.if$newp <- format(round(sigtab.D14.if$padj, digits = 3), scientific = TRUE)
sigtab.D14.if$Treatment <- ifelse(sigtab.D14.if$log2FoldChange >=0, "IM", "IF")
sigtab.D14.if
    
#Plot 'sigtab.D14.if'
deseq.D14.if <- ggplot(sigtab.D14.if, aes(x=reorder(rownames(sigtab.D14.if), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
      geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.if), y=0, label = paste(Family,Genus, sep = ' ')), size=5)+ labs(x="Family Genus")+
      theme(axis.text.x=element_text(color = 'black', size = 13),
            axis.text.y=element_text(color = 'black', size=13, face = 'italic'), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ ggtitle('Differentially Abundant OTUs in IM Relative to IF Group at Nasal Site on Day 14')+ coord_flip() + 
      theme(plot.title = element_text(size = 20, hjust=0.5), legend.text = element_text(size=12), legend.title = element_text(size=13)) +
      scale_fill_manual(labels = c('IF', 'IM'), values = c('#E69F00', '#56B4E9'))
deseq.D14.if
    
#Add OTU and treatment group comparisons columns to 'sigtab.D14.if'
sigtab.D14.if$OTU <- rownames(sigtab.D14.if)
sigtab.D14.if$comp <- 'D14_nasal_IMvsIF'
    
#Add all genera that were differentially abundant between IM and IF treatment groups from 'sigtab.D14.if' to 'final.nonsigtab'
final.nonsigtab <- rbind(final.nonsigtab, sigtab.D14.if)
    
#Write 'final.nonsigtab' into csv file
write.csv(final.nonsigtab, file= "Nasal_FinalDiffAbundGenus_NonSignificantDays.csv")
    
#######################################################################################################
    
######### Plots of Differentially Abundant Nasal Families and Genera Combined for Each Pairwise Comparison of Treatment Groups ########
    
#IF and NON Log2fold plot Part A
final_fc <- sigtab.D0.fc
final_fc <- rbind(final_fc, sigtab.D4.fc, sigtab.D7.fc, sigtab.D11.fc, sigtab.D14.fc)
final_fc$Family_Genus <- paste(final_fc$Family, final_fc$Genus) #create new column with Family_Genus
fc_plot <- ggplot(final_fc, aes(x=Family_Genus, log2FoldChange, fill = comp)) +
      geom_bar(stat='identity') +
      labs(x="Family Genus", y = "Total log2 Fold Change") +
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=7), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ 
      coord_flip() +
      ggtitle('Differentially Abundant Nasal Wash Families and Genera between IF and NON Groups') + 
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13))
fc_plot
write.csv(final_fc, file= "FS1b_FinalDiffAbundNasalGenus_FC.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_FC.csv in a spreadsheet editor by removing all genera except for Blautia, Lachnospiraceae_unclassified,Roseburia, Lactobacillus, Acinetobacter, Actinobacillus, and Streptococcus.
#Saved modified file as FS1b_FinalDiffAbundNasalGenus_FC_final.csv
    
#IF and NON Log2fold plot Part B
fc <- read.csv('FS1b_FinalDiffAbundNasalGenus_FC_final.csv', header = TRUE, sep = ",")
head(fc[,1:10])
colnames(fc)
fc$DayComp <- sub('_[A-Za-z]+', '\\2', fc$comp)
unique(fc$DayComp)
fc$Day <- sub('_[A-Za-z]+', '\\2', fc$DayComp)
unique(fc$Day) #D4"  "D7"  "D11"
fc$Day = factor(fc$Day, levels=c("D4","D7", "D11"))
levels(fc$Day) #"D4"  "D7"  "D11"
(fc_logfoldplot <- ggplot(data=fc, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, ncol = 2, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 4)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IF"), values = c('#999999', '#E69F00')))
fc_logfoldplot <- fc_logfoldplot + theme(strip.text = element_text(size= 15, face='italic'))
fc_logfoldplot
write.csv(final_fc, file= "FS1b_FinalDiffAbundNasalGenus_FC_final.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_FC_final.csv by removing Roseburia genus and saving as 
#FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv
    
#IF and NON Log2fold Plot Part C
fc1 <- read.csv('FS1b_FinalDiffAbundNasalGenus_FCnoRoseburia_final.csv', header = TRUE, sep = ",")
head(fc1[,1:10])
colnames(fc1)
fc1$DayComp <- sub('_[A-Za-z]+', '\\2', fc1$comp)
unique(fc1$DayComp)
fc1$Day <- sub('_[A-Za-z]+', '\\2', fc1$DayComp)
unique(fc1$Day) #D4"  "D7"  "D11"
fc1$Day = factor(fc1$Day, levels=c("D4","D7", "D11"))
levels(fc1$Day) #"D4"  "D7"  "D11"
(fc1_logfoldplot <- ggplot(data=fc1, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, ncol = 2, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 3)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IF"), values = c('#999999', '#E69F00')))
fc1_logfoldplot <- fc1_logfoldplot + theme(strip.text = element_text(size= 13, face='italic'))
fc1_logfoldplot
    
#IM and NON Log2fold Plot Part A
final_ic <- sigtab.D0.ic
final_ic <- rbind(final_ic, sigtab.D4.ic, sigtab.D7.ic, sigtab.D11.ic, sigtab.D14.ic)
final_ic$Family_Genus <- paste(final_ic$Family, final_ic$Genus) #create new column with Family_Genus
ic_plot <- ggplot(final_ic,  aes(x=Family_Genus, log2FoldChange, fill = comp)) +
      geom_bar(stat='identity') +
      labs(x="Family Genus", y = "Total log2 Fold Change") +
      theme(axis.text.x=element_text(color = 'black', size = 10),
            axis.text.y=element_text(color = 'black', size=7), 
            axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15))+ 
      coord_flip() +
      ggtitle('Differentially Abundant Nasal Wash Families and Genera between IM and NON Groups') + 
      theme(plot.title = element_text(size = 20), legend.text = element_text(size=12), legend.title = element_text(size=13))
ic_plot
write.csv(final_ic, file= "FS1b_FinalDiffAbundNasalGenus_IC.csv")
    
#Modified FS1b_FinalDiffAbundNasalGenus_IC.csv in a spreadsheet editor by removing all genera except for Actinobacillus and Streptococcus.
#Saved modified file as FS1b_FinalDiffAbundNasalGenus_IC_final.csv
    
#IM and NON Log2fold Plot Part B
ic <- read.csv('FS1b_FinalDiffAbundNasalGenus_IC_final.csv', header = TRUE, sep = ",")
head(ic[,1:10])
colnames(ic)
ic$DayComp <- sub('_[A-Za-z]+', '\\2', ic$comp)
unique(ic$DayComp)
ic$Day <- sub('_[A-Za-z]+', '\\2', ic$DayComp)
unique(ic$Day) #"D0"  "D4"  "D7"  "D11" "D14"
ic$Day = factor(ic$Day, levels=c("D0", "D4","D7", "D11", "D14"))
levels(ic$Day) #"D0"  "D4"  "D7"  "D11" "D14"
(ic_logfoldplot <- ggplot(data=ic, aes(x=Day, y=log2FoldChange, fill=Treatment)) +
        geom_bar(stat = 'identity', position="dodge") +
        facet_wrap(~Genus, scales = "free") + ylab('log2-fold change') +
        theme_gray()+
        theme(plot.title = element_text(hjust = 2)) +
        theme(axis.line = element_line()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=13),
              axis.text.y = element_text(size=13), 
              axis.title.x = element_text(size=15), 
              axis.title.y = element_text(size=15),
              legend.text=element_text(size=15), 
              legend.title=element_text(size=15)) +
        scale_fill_manual(labels = c("NON", "IM"), values = c('#999999', '#56B4E9')))
ic_logfoldplot <- ic_logfoldplot + theme(strip.text = element_text(size= 15, face='italic'))
ic_logfoldplot
    
#Combine plots 'fc1_logfoldplot' and 'ic_logfoldplot'
ggtwo=plot_grid(fc1_logfoldplot, ic_logfoldplot, labels = c("A", "B"))
ggtwo
    
#Save 'ggtwo' as a .tiff for publication, 500dpi
ggsave("NasalDESeq.tiff", plot=ggtwo, width = 15, height = 5, dpi = 500, units =c("in"))

#################################################################################################################################################################################################

#EIGHTH SECTION: Nasal and Tonsil Microbiota: Genus Abundance
    
#Purpose: This code generates a list of percent total genera found in each treatment group per day for each tissue and creates a bar graph plot of the data 
    
#Files needed:
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared
#FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy
#FS1babundsingleton2000metadata.csv

#Load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library("ggsci")
    
otu <- import_mothur(mothur_shared_file = 'FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'FS1bfinal.outsingletons.abund.opti_mcc.0.03.cons.taxonomy')
shared <- read.table('FS1bfinal.outsingletons.abund.opti_mcc.0.03.subsample.shared', header = TRUE)
meta <- read.table('FS1babundsingleton2000metadata.csv', header = TRUE, sep = ",")
head(meta)
colnames(meta)[1] <- 'group' 
#Rename first column of "meta" as "group" temporarily. Will use "group" to set as rownames later and remove the "group" column
meta$Day<- gsub("D", "", meta$Day) #Remove "D"
meta$group <- as.character(meta$group)
head(meta)
phy_meta <- sample_data(meta) 
rownames(phy_meta) <- phy_meta$group
head(phy_meta)
phy_meta <- phy_meta[,-1]
head(phy_meta)
    
#Create phyloseq-class objects with "otu" and "taxo"
FS1b <- phyloseq(otu, taxo)
FS1b <- merge_phyloseq(FS1b, phy_meta)  #This combines the 'phy_meta' metadata with 'FS1b' phyloseq object
colnames(tax_table(FS1b)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_sums(FS1b) #Calculate the sum of all OTUs for each sample. Almost all samples have 2000 sequences
FS1b <- prune_taxa(taxa_sums(FS1b) > 2, FS1b)  #Removes OTUs that occur less than 2 times globally
FS1b.genus <- tax_glom(FS1b, 'Genus')
phyla_tab <- as.data.frame(t(FS1b.genus@otu_table)) #Transpose 'Fs1b.genus' by "otu_table"
head(phyla_tab)
FS1b.genus@tax_table[,6]
colnames(phyla_tab) <- FS1b.genus@tax_table[,6] #Replace column names in phyla_tab from Otuxxxx with Genus names
phyla_tab2 <- phyla_tab/rowSums(phyla_tab) #Calculate the proportion of specific phyla per phyla column in 'phyla_tab'
head(phyla_tab2)
phyla_tab2$group <- rownames(phyla_tab2) #Create new column called "group" in 'phyla_tab2' containing rownames
head(phyla_tab2)
fobar <- merge(meta, phyla_tab2, by = 'group') #Merge 'meta' with 'phyla_tab2' by "group"
head(fobar)
fobar.gather <- fobar %>% gather(Genus, value, -(group:Treatment))  #This converts 'fobar' to long-form dataframe. This is handy for using ggplot faceting functions, check out tidyverse tutorials
#This also created new columns "Genus", "value"; it added columns "group" through "Treatment" before "Genus" and "value"
head(fobar.gather)
    
#Check to see if there is an extra "group" column. If so, run the next set of commands (up to "head(fobar2)") and remove appropriate column
which(colnames(phyla_tab2)=="group") #Results say column 237 is "group" column
phyla_tab3 <- phyla_tab2[,-237] #Drop the 237th column
phyla_tab4 <- phyla_tab3[,colSums(phyla_tab3)>0.1] #Keep the columns that have greater than 0.1 value
phyla_tab4$group <- rownames(phyla_tab4) #Rename rownames as "group"
fobar2 <- merge(meta, phyla_tab4, by = 'group')
head(fobar2)
    
#To see how many TT are in meta$Tissue: 
length(which(meta$Tissue== "TT")) 
#Output:
#35
    
fobar2.gather$day <- NULL
    
#Reorder days 0-14 in 'fobar2.gather' plot
levels(sample_data(fobar2.gather)$Day)
fobar2.gather$Day <- factor(fobar2.gather$Day, levels=c("D0", "D4", "D7", "D11", "D14"))
head(fobar2.gather$Day)
    
#Create "All" column with "Day", "Treatment" and "Tissue" in 'fobar2.gather'
fobar2.gather$All <- paste(fobar2.gather$Day, fobar2.gather$Treatment, fobar2.gather$Tissue, sep = '_')
    
#Count the number of unique items in 'fobar2.gather'. We're interested in the total unique number of genera
fobar2.gather %>% summarise_each(funs(n_distinct)) #90 total unique genera
fobar2.gather <- fobar2.gather %>% group_by(All) %>% mutate(value2=(value/(length(All)/90))*100) #90 refers to number of Genera
    
#Label nasal genera that have less than 2% abundance as "Other"
nasalgen$More.than.2=as.character(nasalgen$Genus)
str(nasalgen$More.than.2)
nasalgen$More.than.2[nasalgen$Percent.abundance<2]<-"Other"
    
#Rename treatment group names
nasalgen$Treatment2 = nasalgen$Treatment
nasalgen$Treatment2 <- as.character(nasalgen$Treatment2)
nasalgen$Treatment2[nasalgen2$Treatment2 == 'Control'] <- "NON"
nasalgen$Treatment2[nasalgen$Treatment2 == 'Injected'] <- "IM"
nasalgen$Treatment2[nasalgen$Treatment2 == 'Infeed'] <- "IF"
write.csv(nasalgen, file = "Nasal_GenusPercentAbundanceAbove2percent.csv")
    
#To make sure each day's bar added up to 100%, modify Nasal_GenusPercentAbundanceAbove2percent.csv in a spreadsheet editor and save as "Nasal_genus.csv"
    
#Create nasal genera plot
nasalgen2 = read.csv("Nasal_genus.csv", header = TRUE)
nasalgen2.plot <- ggplot(data=nasalgen2, aes(x=Treatment2, y=Percent.abundance, fill=Genus)) +
      geom_bar(stat = 'identity') +
      facet_grid(~Day) + ylab('Relative abundance (%) at nasal site') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
            axis.title.x = element_blank(),
            strip.text = element_text(size= 15),
            axis.text.y = element_text(size=14), 
            axis.title.y = element_text(size=14), 
            legend.text=element_text(size=14), 
            legend.title=element_text(size=14)) +
      scale_fill_igv(name = "Genus") +
      theme(legend.direction = "vertical")
nasalgen2.plot <- nasalgen2.plot + guides(fill= guide_legend(ncol = 1))
nasalgen2.plot <- nasalgen2.plot + theme(legend.text = element_text(face = 'italic'))
nasalgen2.plot
    
#Label tonsil genera that have less than 2% abundance as "Other"
tonsilgen$More.than.2=as.character(tonsilgen$Genus)
str(tonsilgen$More.than.2)
tonsilgen$More.than.2[tonsilgen$Percent.abundance<2]<-"Other"
write.csv(tonsilgen, file = "Tonsil_GenusPercentAbundanceAbove2percent.csv")
    
#To make sure each day's bar add up to 100%, modify Tonsil_GenusPercentAbundanceAbove2percent.csv in a spreadsheet editor and save as "Tonsil_genus.csv"
    
#Tonsil
tonsilgen2 = read.csv("Tonsil_genus.csv", header = TRUE)
    
#Renaming injected and infeed to parenteral and medicated feed for tonsil
tonsilgen2$Treatment2 = tonsilgen2$Treatment.group
tonsilgen2$Treatment2 <- as.character(tonsilgen2$Treatment2)
tonsilgen2$Treatment2[tonsilgen2$Treatment2 == 'Control'] <- "NON"
tonsilgen2$Treatment2[tonsilgen2$Treatment2 == 'Injected'] <- "IM"
tonsilgen2$Treatment2[tonsilgen2$Treatment2 == 'Infeed'] <- "IF"
    
#Create tonsil genera plot
tonsilgen2.plot <- ggplot(data=tonsilgen2, aes(x=Treatment2, y=Percent.abundance, fill=Genus)) +
      geom_bar(stat = 'identity') +
      facet_grid(~Day) + ylab('Relative abundance (%) at tonsil site') +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
            axis.title.x = element_blank(),
            strip.text = element_text(size= 15),
            axis.text.y = element_text(size=14), 
            axis.title.y = element_text(size=14), 
            legend.text=element_text(size=14), 
            legend.title=element_text(size=14)) +
      scale_fill_igv(name = "Genus")
tonsilgen2.plot <- tonsilgen2.plot + theme(legend.text = element_text(face = 'italic'))
tonsilgen2.plot
    
#Combine nasal and tonsil genera plots
nw2tt2.combo <- plot_grid(nasalgen2.plot,tonsilgen2.plot, labels = "AUTO", rel_widths = c(1.2, 1))
nw2tt2.combo
    
#Save 'nw2tt2.combo' as a .tiff for publication, 500dpi
ggsave("NasalTonsilGenera.tiff", plot=nw2tt2.combo, width = 15, height = 7, dpi = 500, units =c("in"))

#################################################################################################################################################################################################

#NINTH SECTION: Nasal Oxytetracycline Levels
    
#Purpose: This code generates a box & whisker plot depicting the nasal oxytetracycline levels measured from each treatment group
    
#Files needed:
#Nasal_results.csv
#Metadata.csv


#Load libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
    
nasal <- read.csv('Nasal_results.csv', stringsAsFactors = FALSE)  
meta <- read.csv('Metadata.csv', stringsAsFactors = FALSE)      
inj <- meta$Pig.ID.Numbers[grep('IM', meta$NON.treatment)]  #Creates a vector of pig numbers that have 'IM' in their treatment column
oral <- meta$Pig.ID.Numbers[grep('IF', meta$NON.treatment)] #Creates a vector of pig numbers that have 'IF' in their treatment column
control <- meta$Pig.ID.Numbers[grep('NON', meta$NON.treatment)] #Creates a vector of pig numbers that have 'NON' in their treatment column
nasal$Time.Point <- as.numeric(gsub('D', '', nasal$Time.Point)) #Replaces the 'D' in the time column with '' (nothing)
nasal
nasal.melt <- melt(nasal, id.vars = c(1,2), measure.vars = c(3)) #Converts to long dataframe format for easy plotting
nasal.melt                                                       #Just checking on the new long dataframe
colnames(nasal.melt) <- c('pig', 'day', 'tissue', 'concentration')  #Changing the column names to something more understandable
nasal.melt$concentration <- gsub('NF', 0, nasal.melt$concentration) #Replaces 'NF' with 0
nasal.melt$concentration <- gsub('<', '', nasal.melt$concentration) #Replaces '<' with '' (nothing)
nasal.melt$concentration <- as.numeric(nasal.melt$concentration)    #Forces the concentration column to be numeric 
nasal.melt <- na.exclude(nasal.melt)                                #Removes NAs which were introduced in the previous line
nasal.melt$treatment <- NA                                          #Creates a new column 'treatment'
nasal.melt
nasal.melt$treatment[nasal.melt$pig %in% control] <- 'NON'   #Assigns 'NON' to pigs whose numbers are in the control vector we generated earlier
nasal.melt$treatment[nasal.melt$pig %in% inj] <- 'IM'        #Assigns 'IM' to pigs whose numbers are in the control vector we generated earlier
nasal.melt$treatment[nasal.melt$pig %in% oral] <- 'IF'       #Assigns 'IF' to pigs whose numbers are in the control vector we generated earlier
nasal.melt$dayXtreat <- paste(nasal.melt$day, nasal.melt$treatment, sep = '_') #Creates a dayXtreatment column; it's just the day and treatment columns pasted together
nasal.melt$day <- factor(nasal.melt$day)    #Makes the day column a factor (this means it is categorical data, not continuous data)
nasal.melt
    
#This section uses tidyverse type pipes, filters, and ggplot to produce boxplots
nasalplot <- ggplot(nasal.melt, aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
      geom_boxplot() + ylab('Nasal oxytetracycline level (ng/mL)') + xlab('Day')
nplot <-  nasalplot + theme_gray(base_size = 14) + guides(fill=guide_legend(title="Treatment group")) +
      scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
      theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size=14),
            legend.text=element_text(size=14))
nplot
    
#Save 'nplot' as a .tiff for publication, 500dpi
ggsave("Nasal_AntibioticLevels.tiff", plot=nplot, width = 10, height = 5, dpi = 500, units =c("in"))