detach("package:seqinr", unload=TRUE)
library(plyr)
library(msa)

orthos <- read.csv("data/LRAD3_orthologues.csv", stringsAsFactors = F)
orthos$ensemblID <- gsub("LDLRAD3.*$","", orthos$Ensembl.identifier...gene.name, ignore.case = T)
orthos <- orthos[!grepl(" ", orthos$ensemblID),]
orthos$key <- gsub("^.*\\(", "",orthos$Species)
orthos$key <- gsub("\\)$","", orthos$key)

ds <- listDatasets(ensembl)
ds$key <- gsub(" genes.*$", "", as.character(ds$description))

orthos <- merge(orthos,ds, by = "key")

library(biomaRt)
ensembl <- useEnsembl("ensembl")
seqs <- ddply(orthos, "ensemblID", function(x) {
    ensembl <- useDataset(as.character(x$dataset), ensembl)
    prot <- getSequence(id = x$ensemblID, 
                        type = "ensembl_gene_id", 
                        seqType = "peptide",
                        mart = ensembl)
    x$sequence <- prot$peptide[1]
    x
})

ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)
humanSeq <- getSequence(id = "LDLRAD3", 
                       type ="hgnc_symbol",
                       seqType = "peptide",
                       mart = ensembl)

library(seqinr)

write.fasta(sequences = as.list(seqs$sequence), names = seqs$Species, "data/LRAD3_orthologues.fasta")

msaSeqs <- readAAStringSet("data/LRAD3_orthologues.fasta")

alignment <- msa(msaSeqs)

print(alignment, show="complete")
# msaPrettyPrint(alignment, output="pdf", showNames="none",
#                showLogo="none", askForOverwrite=FALSE, verbose=FALSE)