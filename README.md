# RNA-to-DNA-Reverse-Transcription-and-Nucleotide-Frequency-Analysis
This project aims to showcase the reverse transcription of  RNA sequences into DNA sequences and analyze the nucleotide frequencies, including di-nucleotide and tri-nucleotide frequencies.

# Install Required Packages
Before running the program, make sure to install the necessary R packages:<br>
```
install.packages("qPCRtools")
install.packages("Biostrings")
install.packages("magrittr")
library(qPCRtools)
library(magrittr)
library(Biostrings)

```

# Input RNA Sequences
Provide the RNA sequences in fasta format, RNA-sequence1.fna and RNA-sequence2.fna.

# Run the Program
Execute the R script provided in the repository to perform the following steps: <br>
Read RNA sequences from fasta files.<br>
Reverse transcribe RNA sequences to DNA.<br>
Calculate di-nucleotide and tri-nucleotide frequencies.<br>
Compare nucleotide compositions between sequences.<br>
Identify differences exceeding a threefold threshold.<br>

