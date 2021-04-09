# Fork of Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273   					 

User sanxiyn on Hacker News [said](https://news.ycombinator.com/item?id=26630537):

> Both sequences express the same protein.
>
> Sequences are different because they are differently codon optimized. See https://en.wikipedia.org/wiki/Codon_usage_bias, especially "Effect on transcription or gene expression" section. 

The file `proteindiff.py` has the goal of checking the veracity of this statement. I'm not 100% confident that it currently is extracting the correct [ORF](https://en.wikipedia.org/wiki/Open_reading_frame).  It currently translates to 1050 amino acids and I would expect closer to 4,000 (bases) / 3 (bases/codon) = 1333 (amino acids).  I'm more sure the start codon is being found correctly. Within the 1050 amino acids that are translated from the DNA sequences, I am finding 20 amino acids different between them.  A majority of these are represented in ascii as `*` characters, which denote a stop codon.  Additionally, Biopython warns me that I'm invoking translation on a non power of three length string.  Both of these lead me to suspect the reading frame may be wrong.  So take all of this with a grain of salt until I investigate further.

Now on to the original readme:

# Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273   					 
##version 0.2Beta 03/30/21: (update intended to (i) clarify the clinical and research importance of sequence information and strand topology measurements, and (ii) clarify that the mRNA sequence is not a recipe to produce vaccine)##

Dae-Eun Jeong, Matthew McCoy, Karen Artiles, Orkan Ilbay, Andrew Fire*, Kari Nadeau, Helen Park, Brooke Betts, Scott Boyd, Ramona Hoh, and Massa Shoura*

Departments of Pathology, Genetics, Pediatrics, and Medicine, Stanford University School of Medicine and Veterans Affairs Palo Alto Medical Center
*Correspondence: afire@stanford.edu and/or massa86@stanford.edu

Abstract: RNA vaccines have become a key tool in moving forward through the challenges raised both in the current pandemic and in numerous other public health and medical challenges.  With the rollout of vaccines for COVID-19, these synthetic mRNAs have become broadly distributed RNA species in numerous human populations.  Despite their ubiquity, sequences are not always available for such RNAs.  Standard methods facilitate such sequencing.  In this note, we provide experimental sequence information for the RNA components of the initial Moderna (https://pubmed.ncbi.nlm.nih.gov/32756549/) and Pfizer/BioNTech (https://pubmed.ncbi.nlm.nih.gov/33301246/) COVID-19 vaccines, allowing a working assembly of the former and a confirmation of previously reported sequence information for the latter RNA.

Objective: Sharing of sequence information for broadly used therapeutics has substantial benefit in design of improved clinical tools and precise diagnostics.  As examples of such applications, (i) Any medical or public health analysis relying on high throughput sequencing data to track SARS-COV2 and its variants would benefit from knowledge of vaccine sequences in order to distinguish RNA sequencing reads coming from the vaccine from those of viral origin, and (ii) Diagnostic labs designing nucleic acid surveillance tests (like PCR or LAMP assays) can benefit from vaccine sequence information to avoid confusion between vaccinated and infected test subjects when analyzing assay results. 

Description: For this work, RNAs were obtained as discards from the small portions of vaccine doses that remained in vials after immunization; such portions would have been required to be otherwise discarded and were analyzed under FDA authorization for research use.  To obtain the small amounts of RNA needed for characterization, vaccine remnants were phenol-chloroform extracted using TRIzol Reagent (Invitrogen), with intactness assessed by Agilent 2100 Bioanalyzer before and after extraction.  

Although our analysis mainly focused on RNAs obtained as soon as possible following discard, we also analyzed samples which had been refrigerated (~4 ℃) for up to 42 days with and without the addition of EDTA.  Interestingly a substantial fraction of the RNA remained intact in these preparations.  We note that the formulation of the vaccines includes numerous key chemical components which are quite possibly unstable under these conditions-- so these data certainly do not suggest that the vaccine as a biological agent is stable.  But it is of interest that chemical stability of RNA itself is not sufficient to preclude eventual development of vaccines with a much less involved cold-chain storage and transportation. 

For further analysis, the initial RNAs were fragmented by heating to 94℃, primed with a random hexamer-tailed adaptor, amplified through a template-switch protocol (Takara SMARTer Stranded RNA-seq kit), and sequenced using a MiSeq instrument (Illumina) with paired end 78-per end sequencing.  As a reference material in specific assays, we included RNA of known concentration and sequence (from bacteriophage MS2).

From these data, we obtained partial information on strandedness and a set of segments that could be used for assembly.  This was particularly useful for the Moderna vaccine, for which the original vaccine RNA sequence was not available at the time our study was carried out.  Contigs encoding full-length spikes were assembled from the Moderna and Pfizer datasets.  The Pfizer/BioNTech data [Figure 1] verified the reported sequence for that vaccine (https://berthub.eu/articles/posts/reverse-engineering-source-code-of-the-biontech-pfizer-vaccine/), while the Moderna sequence [Figure 2] could not be checked against a published reference. 

RNA preparations lacking dsRNA are desirable in generating vaccine formulations as these will minimize an otherwise dramatic biological (and nonspecific) response that vertebrates have to double stranded character in RNA (https://www.nature.com/articles/nrd.2017.243).  Numerous recent advances have resulted in approaches to minimize dsRNA (e.g. https://pubmed.ncbi.nlm.nih.gov/30933724/, https://pubmed.ncbi.nlm.nih.gov/31900329/); nonetheless measurement remains a continued necessity.   In the sequence data that we analyzed, we found that the vast majority of reads were from the expected sense strand.  In addition, the minority of antisense reads appeared different from sense reads in lacking the characteristic extensions expected from the template switching protocol.  Examining only the reads with an evident template switch (as an indicator for strand-of-origin), we observed that both vaccines overwhelmingly yielded sense reads (>99.99%).  Independent sequencing assays and other experimental measurements are ongoing and will be needed to determine whether this template-switched sense read fraction in the SmarterSeq protocol indeed represents actual dsRNA content in the original material.

This work provides an initial assessment of two RNAs that are now a part of the human ecosystem and that are likely to appear in numerous other high throughput RNA-seq studies in which a fraction of the individuals may have previously been vaccinated.

ProtoAcknowledgements: We thank our colleagues here for their help and suggestions (Nimit Jain, Emily Greenwald, Nelson Hall, Lamia Wahba, William Wang, Amisha Kumar, Sameer Sundrani, David Lipman, Marc Salit), and additionally acknowledge numerous colleagues who have discussed and educated us (directly and indirectly) in areas of RNA synthesis enzymology, vaccine design, and software engineering.

Figure 1: Spike-encoding contig assembled from BioNTech/Pfizer BNT-162b2 vaccine.  Although the full coding region is included, the nature of the methodology used for sequencing and assembly is such that the assembled cDNA-derived contig could lack some sequence from the ends of the vaccine RNA and would lack any indication of base modifications.  Within the assembled sequence, this hypothetical sequence shows a perfect match to the corresponding sequence from documents available online derived from manufacturer communications with the World Health Organization [as reported by https://berthub.eu/articles/posts/reverse-engineering-source-code-of-the-biontech-pfizer-vaccine/].. The 5’ end for the assembly matches the start site noted in these documents, while the read-based assembly lacks an interrupted polyA tail (A30(GCATATGACT)A70) that is expected to be present in the mRNA.  

Figure 2: Spike-encoding contig assembled from Moderna mRNA-1273 vaccine.  This is a partial sequence of the vaccine RNA.  Although the full coding region is included, the assembled contig could lack some sequence from the ends of the vaccine RNA.  




