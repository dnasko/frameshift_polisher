![alt text](https://github.com/dnasko/frameshift_polisher/blob/master/images/polisher_logo.png?raw=true "The Frameshift Polisher")

A BLAST-driven pipeline for polishing out pesky frameshifts -- particularly in PacBio reads.

Preamble
--------

Even after error correction, PacBio reads posses a *few* pesky errors. These errors are enough to shift the frame of subsequent peptide sequences. A labmate (Eric Sakowski) was working on a certain marker gene and noticed that the majority of PacBio reads were essentially useless due to these lingering frameshifts. We needed a way to quickly and accurately correct these frameshifts. He would see long stretches where the beginning of the peptides were in frame, then long stretches where it was in another frame, then it may or may not switch back to its original frame.

If we have a read, that has been mostly error corrected, and we know what this gene *ought* to be, then by performing a 6-frame translation and aligning all six frames with other genes similar to it we should be able to see what regions are in what frames.

### What is it capable of ?

As the term "polisher" implies this program is meant to only fix a few lingering frame shifts that are present in a given sequence. This is an inheritied restriction from BLAST as significant peptide alignments will only be found if they are long enough (i.e. not broken apart by frequent frameshifts).

Usage
-----

**Input**: FASTA of nucleotide sequences

**Output**: FASTA of frameshift corrected peptide sequences

`$ perl frameshift_polisher -fasta /Path/to/reads.fasta -db /Path/to/db -work /Path/to/working/directory`


How it works
------------

Okay, so you have a large set of nucleotide sequences that you know possess a few frameshifts. Polisher will take this FASTA file and perform a BLASTX against a user-defined BLASTable datadase. Polisher will then parse the resutls of this BLASTX to frameshift correct the sequences.

The BLASTable database is a peptide BLAST database of (you guessed it) peptide sequences. It can be as small as a set of 2,000 *Cas* peptides or as large as NCBI's NR peptide database. The larger the database, the longer the blast will take, but the more thorough the search will be.

### Brief example

You recieve the following BLAST results whereby your sequence of interest "gene_01" has hit the same subject sequence "subj_01" from your database 4 times. Notice the hits are in 3 frames relative to the subject sequence. The polisher will parse these results to piece together a translated query sequence that will ultimatly remain in the same frame.

| Query   | Subject | Perc ID | q.start | q.end | s.start | s.end | e.value | q.frame |
| ------- | ------- | -------:| -------:| -----:| -------:| -----:| -------:| -------:|
| gene_01 | subj_01 |   82.76 |     455 |   282 |     224 |   281 |   2e-44 |      -3 |  
| gene_01 | subj_01 |   51.02 |     591 |   448 |     178 |   226 |   2e-44 |      -2 |  
| gene_01 | subj_01 |   83.33 |     102 |    31 |     341 |   364 |   2e-44 |      -2 |  
| gene_01 | subj_01 |   80.00 |     247 |   188 |     279 |   340 |   2e-44 |      -1 |  

By merging these four peptide alignmetn together The Polisher will have effectivly cliped out any frameshifts, be them insertions, deletions, or miscalled bases.


