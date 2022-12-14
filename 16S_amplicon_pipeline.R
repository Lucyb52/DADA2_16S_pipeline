
library("dada2"); packageVersion("dada2")
library(ggplot2)
setwd("/Users/lucyburns/Downloads/Course_Docs_F22/Bioinformatics/Assignment 3")
path<-"/Users/lucyburns/Downloads/Course_Docs_F22/Bioinformatics/Assignment 3"
#this lists my files
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#looking at quality of F reads in a plot
plotQualityProfile(fnFs[1:2])
#looking at quality of reverse reads
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#make a judgement call on where to truncate reads - but make sure there is still overlap, 
#Chose 280 and 260

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = c(18,20)) # On Windows set multithread=FALSE
head(out)

#de-replicate (will trim some adn make less computationally heavy downstream)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#run error rates
errF <- learnErrors(filtFs, multithread=TRUE)
#for reverse
errR <- learnErrors(filtRs, multithread=TRUE)

#plot errors to check
plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspect the dada-class object
dadaFs[[1]]
dadaRs[[1]]

#Merge paired reads
#if merge doesn't work then maybe you trimmed too much during quality step
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table 
seqtab <- makeSequenceTable(mergers)
#returns lengths vary error
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim<-removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE,verbose = TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#transpose seqtab.nochim data if you want to look at it as a column
flipped_seqtab.nochim<-as.data.frame(t(seqtab.nochim))

#look at the reads that have made it through the steps
getN <- function(x) sum (getUniques(x))
track <-cbind(out,sapply(dadaFs,getN),sapply(dadaRs,getN),sapply(mergers,getN),rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy, download Silva training, make sure to move to working directory (path)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
list()

#optional additional species assignment training
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

#inspect taxa assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save as a csv file (not in dada2 protocol)
write.csv (taxa, file = "Ass3.csv_taxa.csv")
write.csv (seqtab.nochim, file="Ass3.seqtab.nochim.csv")
write.csv (flipped_seqtab.nochim, file="Ass3Flipped_seqtab.nochim.csv")

#save ASV and count as one sheet flipped seqtab no chim file with your taxa data
OTUabund<-cbind(flipped_seqtab.nochim, taxa)
write.csv(OTUabund, file="OTUabund2.csv")

#install phyloseq package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", force=TRUE)

BiocManager::install("phyloseq",force=TRUE)

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(Biostrings); packageVersion("Biostrings")

#construct a dataframe from our filenames:
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(samples.out)
rownames(samdf) <- samples.out

#make a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

#table with ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#plot
plot_bar(ps)
plot_bar(ps, fill="Phylum") +labs(y= "Abundance", x = "Sample") + ggtitle ("Relative abundance of phyla in a pumice sample")

#take away bars
p=plot_bar(ps)
p + geom_bar(aes(fill=Phylum), stat="identity", position="stack")

#transform to relative abundance
relative<- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
relative
#plot
Phylum_graph <- plot_bar(relative, fill="Phylum") +ylab ("Relative Abundance")
Phylum_graph
#get rid of bars
Phylum_graph + geom_bar(aes(fill=Phylum), stat="identity", position="stack") +labs(y= "Relative Abundance", x = "Sample") + ggtitle ("Relative abundance of phyla in a pumice sample")

#do it for order
plot_bar(ps)
plot_bar(ps, fill="Order") +labs(y= "Abundance", x = "Sample") + ggtitle ("Relative abundance of order in a pumice sample")

#plot relative order abundance
Order_graph <- plot_bar(relative, fill="Order") +ylab ("Relative Abundance")
Order_graph
#take out black bars
Order_graph + geom_bar(aes(fill=Order), stat="identity", position="stack") +labs(y= "Relative Abundance", x = "Sample") + ggtitle ("Relative abundance of order in a pumice sample")

#save where we are
saveRDS(ps, "ps.rds")

#Tax_glom will conglomerate all the reads of an identical taxa together ( here, we want to count everything that is the same phylum)
ps_phylum <- tax_glom(ps, "Phylum")

#Now we are just taking this, and turning it into relative abundance as we did above
ps1_phylum_relabun <- transform_sample_counts(ps_phylum, function(ASV) ASV/sum(ASV))

#The psmelt function of phyloseq will make a dataframe of our phyloseq data for us, we also need it to be a factor
taxa_abundance_table_phylum <- psmelt(ps1_phylum_relabun)
taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)

#graph in ggplot
ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample,y=Abundance*100,))+geom_col(position="stack", stat="identity")+aes(fill=Phylum)+ 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Samples", y="Phylum Relative Abundance (%)", size = "Relative Abundance(%)")  + 
  theme(axis.text.x = element_text(face="plain", color="Black", size=10, angle=90),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.text.y = element_text(face="bold", color="Black", size=10, angle=0)) 
#diff figure
ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample,y=Phylum ))+geom_point(scales="free_x")+(aes(size = Abundance,colour=Phylum))+ggtitle("Permafrost microbiota ")+scale_size_area()+ 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Samples", y="Phylum Relative Abundance (%)", size = "Relative Abundance(%)")  + 
  theme(axis.text.x = element_text(face="plain", color="Black", size=10, angle=90),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.text.y = element_text(face="bold", color="Black", size=10, angle=0)) 
#fix size?
taxa_abundance_table_phylum[taxa_abundance_table_phylum == 0] <- NA

ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample,y=Phylum ))+geom_point(scales="free_x")+(aes(size = Abundance,colour=Phylum))+ggtitle("Permafrost microbiota ")+scale_size_area()+ 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Samples", y="Phylum Relative Abundance (%)", size = "Relative Abundance(%)")  + 
  theme(axis.text.x = element_text(face="plain", color="Black", size=10, angle=90),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.text.y = element_text(face="bold", color="Black", size=10, angle=0)) 
#plot again
ggplot(data=taxa_abundance_table_phylum,mapping=aes(x=Sample,y=Phylum ))+geom_point(scales="free_x")+(aes(size = Abundance,colour=Phylum))+ggtitle("Pumice microbiota")+scale_size_area()+ 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Samples", y="Phylum Relative Abundance (%)", size = "Relative Abundance(%)")  + 
  theme(axis.text.x = element_text(face="plain", color="Black", size=8, angle=90),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), #tickmark aesthetics
        axis.text.y = element_text(face="bold", color="Black", size=6, angle=0)) 
