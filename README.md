The Frameshift Polisher
=======================

A BLAST-driven pipeline for polishing out pesky frameshifts -- particularly in PacBio reads.

Preamble
--------

Even after error correction, PacBio reads posses a *few* pesky errors. These errors are enough to shift the frame of subsequent peptide sequences. A labmate (Eric Sakowski) was working on a certain marker gene and noticed that the majority of PacBio reads were essentially useless due to these lingering frameshifts. We needed a way to quickly and accurately correct these frameshifts. He would see long stretches where the beginning of the peptides were in frame, then long stretches where it was in another frame, then it may or may not switch back to its original frame.

If we have a read, that has been mostly error corrected, and we know what this gene *ought* to be, then by performing a 6-frame translation and aligning all six frames with other genes similar to it we should be able to see what regions are in what frames.

What is it capable of ?
-----------------------

As the "polisher" portion of the name implies this program is meant to only fix a few lingering frame shifts that are present in a given sequence. This is an inheritied restriction from BLAST.

Sounds mysterious, how dow it work?

The Frameshift Polisher
-----------------------

Input: FASTA of nucleotide sequences

Output: FASTA of frameshift corrected peptide sequences

Okay, so you have a large set of nucleotide sequences that you know possess a few frameshifts. Polisher will take this FASTA file and perform a BLASTX against a user-defined BLASTable datadase. Polisher will then parse the resutls of this BLASTX to frameshift correct the sequences.

The BLASTable database is a peptide BLAST database of (you guessed it) peptide sequences. It can be as small as a set of 2,000 *Cas* peptides or as large as NCBI's NR peptide database. The larger the database, the longer the blast will take, but the more thorough the search will be.

`query_id_1   subject_id_1 51.26    119	  57	   1  736  380     170     287     4e-44    114    -1      0`
