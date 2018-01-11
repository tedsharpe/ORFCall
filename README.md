# ORFCall

This is a program to interpret the results of a MITE-seq experiment.

To install:

Clone repo and type make.  This is a very simple, one-source-file project with a two-line handcrafted makefile.
You'll need some reasonably modern version of gcc and zlib.

To run:

`ORFCall ref.fasta codons_used.dat reads1.fastq [reads2.fastq...]
 or
ORFCall -p ref.fasta codons_used.dat reads1a.fastq reads1b.fastq [reads2a.fastq reads2b.fastq...]`

The ref.fasta file is a single-contig fasta file describing the amplicon for the MITE-seq experiment.
You should mark the open reading frame by delimiting it with square brackets like this:

`> my.MITE-seq.ref
TTTGGGCCCAAA[AAGAATAAC]TTTGGGCCCAAA`

The codons_used.dat file is a tab-delimited file that starts with a header line enumerating all possible
codon sequences.  The header must be followed by one line for each codon in the ORF, with the value 0, 1, or 2
for each possible codon sequence.  The value 0 means that the planned perturbations do not include the codon.
The value 1 indicates that this codon sequence is expected to appear as a planned perturbation.  And the value 2
means that this is the wild-type codon sequence.
Like this:

`AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA ... TTT
0   1   2   0   1   1   0   0   1   1   0   0   1   0   0   1   1       0
1   0   0   2   0   1   1   0   1   0   1   0   0   1   0   1   1       1
0   2   0   0   1   0   0   1   1   0   0   0   1   1   0   1   0       1`

If the fastq files end with ".gz", they will be unzipped on the fly.
