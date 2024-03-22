install.packages("qPCRtools")
install.packages("Biostrings")
install.packages("magrittr")
library(qPCRtools)
library(magrittr)
library(Biostrings)

# Reading the RNA sequences from the fasta files
rna_seq1 <- readDNAStringSet("RNA-sequence1.fna")
rna_seq2 <- readDNAStringSet("RNA-sequence2.fna")

# Reverse transcribing the RNA sequences to DNA 
dna_seq1 <- reverseComplement(rna_seq1)
dna_seq2 <- reverseComplement(rna_seq2)

cat("RNA-sequence1.fna Forward Strand:", as.character(dna_seq1), "\n")
cat("RNA-sequence2.fna Forward Strand:", as.character(dna_seq2), "\n")

# Calculating di-nucleotide and tri-nucleotide frequencies for RNA-sequence1.fna
di_freq1 <- oligonucleotideFrequency(rna_seq1, width=2)
tri_freq1 <- oligonucleotideFrequency(rna_seq1, width=3)

print(di_freq1)
print(tri_freq1)

# Calculating di-nucleotide and tri-nucleotide frequencies for RNA-sequence2.fna
di_freq2 <- oligonucleotideFrequency(rna_seq2, width=2)
tri_freq2 <- oligonucleotideFrequency(rna_seq2, width=3)

print(di_freq2)
print(tri_freq2)

# Printing the results for RNA-sequence1.fna
cat("RNA-sequence1.fna Di-nucleotide Frequencies:\n")
print(di_freq1, quote=FALSE, right=FALSE)
cat("\nRNA-sequence1.fna Tri-nucleotide Frequencies:\n")
print(tri_freq1, quote=FALSE, right=FALSE)

# Printing the results for RNA-sequence2.fna
cat("\nRNA-sequence2.fna Di-nucleotide Frequencies:\n")
print(di_freq2, quote=FALSE, right=FALSE)
cat("\nRNA-sequence2.fna Tri-nucleotide Frequencies:\n")
print(tri_freq2, quote=FALSE, right=FALSE)


# Calculating total number of oligonucleotides
total1 <- sum(di_freq1)
total2 <- sum(di_freq2)
total3 <- sum(tri_freq1)
total4 <- sum(tri_freq2)

# Calculating percentages of di-nucleotide and tri-nucleotide frequencies for RNA-sequence1.fna
di_perc1 <- di_freq1 / total1 * 100
tri_perc1 <- tri_freq1 / total3 * 100

#Calculating percentages of di-nucleotide and tri-nucleotide frequencies for RNA-sequence2.fna
di_perc2 <- di_freq2 / total2 * 100
tri_perc2 <- tri_freq2 / total4 * 100

print(di_perc1)
print(tri_perc1)

print(di_perc2)
print(tri_perc2)

# Identifying di and tri-nucleotides that differ by at least 3-fold between the two sequences
di_fc <- di_perc2 / di_perc1
tri_fc <- tri_perc2 / tri_perc1

print(di_fc)
print(tri_fc)

di_3fold <- names(di_fc[abs(log2(di_fc)) >= log2(3)])
tri_3fold <- names(tri_fc[abs(log2(tri_fc)) >= log2(3)])

# Checking if any dinucleotides or trinucleotides differ by 3X in percentage
if (is.null(di_3fold) && is.null(tri_3fold)) {
  cat("No dinucleotides or trinucleotides differ by 3X in percentage between the two RNA sequences.")
} else {
  # Printing the results
  if (!is.null(di_3fold)) {
    cat("Dinucleotides that differ by 3X in percentage:\n")
    for (di in di_3fold) {
      cat(di, "\tRNA Sequence 1: ", di_perc1[di], "%\tRNA Sequence 2: ", di_perc2[di], "%\n")
    }
  }
  
  if (!is.null(tri_3fold)) {
    cat("\nTrinucleotides that differ by 3X in percentage:\n")
    for (tri in tri_3fold) {
      cat(tri, "\tRNA Sequence 1: ", tri_perc1[tri], "%\tRNA Sequence 2: ", tri_perc2[tri], "%\n")
    }
  }
}









