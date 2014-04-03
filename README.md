The Frameshift Polisher
=======================

A BLAST-driven pipeline for polishing out pesky frameshifts -- particularly in PacBio reads.

Preamble
--------

Even after error correction, PacBio reads posses a *few* pesky errors. These errors are enough to shift the frame of subsequent peptide sequences. A labmate (Eric Sakowski) was working on a certain marker gene and noticed that the majority of PacBio reads were essentially useless due to these lingering frameshifts. We needed a way to quickly and accurately correct these frameshifts. He would see long stretches where the beginning of the peptides were in frame, then long stretches where it was in another frame, then it may or may not switch back to its original frame.

If we have a read, that has been mostly error corrected, and we know what this gene *ought* to be, then by performing a 6-frame trnalsiation and aligning all six frames with other genes similar to it we should be able to see what regions are in what frames.

What It Does?
-------------

As the "polisher" portion of the name implies this program is meant to only fix a few lingering frame shifts that are present in a given sequence.