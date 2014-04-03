The Frameshift Polisher
=======================

A BLAST-driven pipeline for polishing out pesky frameshifts -- particularly in PacBio reads.

Description
-----------

Even after error correction, PacBio reads posses a *few* nagging errors. These errors are enough to shift the frame of subsequent peptide sequences. A labmate (Eric Sakowski) was working on a certain marker gene and noticed that the majority of PacBio reads were essentially useless due to these lingering frameshifts. We needed a way to quickly and accurately correct these frameshifts, this pipeline will do just that.
