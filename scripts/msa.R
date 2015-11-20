library(plyr)
library(msa)

orthos <- read.csv("data/LRAD3_orthologues.csv", stringsAsFactors = F)
orthos$ensemblID <- gsub("LDLRAD3.*$","", orthos$Ensembl.identifier...gene.name, ignore.case = T)
orthos <- orthos[!grepl(" ", orthos$ensemblID),]
orthos$key <- gsub("^.*\\(", "",orthos$Species)
orthos$key <- gsub("\\)$","", orthos$key)

library(biomaRt)

ensembl <- useEnsembl("ensembl")

ds <- listDatasets(ensembl)
ds$key <- gsub(" genes.*$", "", as.character(ds$description))
orthos <- merge(orthos,ds, by = "key")

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
hSap <- data.frame(Species = "Human (Homo sapien)", sequence = humanSeq$peptide[1])
seqs <- merge(hSap, seqs, all = T)

library(seqinr)

seqList <- as.list(seqs$sequence)
species <- seqs$Species
write.fasta(sequences = seqList, names = species, "data/LRAD3_orthologues.fasta")

msaSeqs <- readAAStringSet("data/LRAD3_orthologues.fasta")

alignment <- msa(msaSeqs, order = "input")

print(alignment, show="alignment")
msaPrettyPrint(alignment, output="pdf", showNames="left",
                showLogo="none", askForOverwrite=FALSE, verbose=F, 
               alFile = "test.fasta",
               consensusThreshold = 50,
               shadingMode = "similar")

detach("package:seqinr", unload=TRUE)