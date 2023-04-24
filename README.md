# pair_maps01
maps Programs for GeF-seq (Genome Footprinting) analysis with paired-end information

MacOSX Version compiled under MacOS Sierra version 10.12.3

Download the folowing files
1. mkindex
2. mpsmap
3. pair_map
4. mapsview_b01.app
5. mapsview_tutorial.pdf

Open terminal, and place the download file to where you want, for instance, to your home/current directory.
mv ~/Downloads/pair_map .
mv ~/Downloads/mpsmap .
mv ~/Downloads/mkindex .

You may also need to set path to this directory, if you will carry out data analysis in other directory.

------------------------------------------------------
Tutorial for a paired end GeF-seq analysis: 
------------------------------------------------------

Please refer to the following URL for the General description of Gef-seq analysis and the output files. Gef-seq analysis using maps software 
Following is an example of a GeF-seq analysis using pair_map. 

1. Preparation: 

Suppose you have four kinds of data files for analysis. 
1. Genome sequence file in fasta format (genome.fasta). 
2. Genome annotation file in genbank format (genome.gb). 
3. Genome sequencing data file by Illumina sequencer in fastq format (genome_seq.fastq). 
4. Paired end GeF-seq sequencing data files from Illumina sequencer in fastq format (gef_seq_R1.fastq and gef_seq_R2.fastq). 
Read length of about 75 base pairs is assumed.

First, create an index file for reference genome. 
% ./mkindex -11 genome.fasta
The option -11 specifies the length of the index as 11.
This will create an index file genome.fasta.indx11


2. Obtaining the ORI-TER bias parameters: 

Follow the sections 2. to 4. of the Gef-seq description page to obtain the ORI-TER bias parameter.


3. Paiered end Gef-seq mapping: 

pair_map program replace the analysis by pmapsr program used for the regular Gef-seq analysis. In this paired end sequence analysis, we assume the initial base of each reads of a read pair correspond to the both ends of the sequencing segments, that is, the sequence segment bound by the DNA binding proteins. 
To analyze, first, map the GeF-seq data on the reference as following.
% mpsmap -a 8 -i 11 -r 36 -h 35 -map -job r36h35i11 genome.fasta gef_seq_R1.fastq
% mpsmap -a 8 -i 11 -r 36 -h 35 -map -job r36h35i11 genome.fasta gef_seq_R2.fastq
Then you will get the hit position files (gef_seq_R1.fastq_r36h35i11.hit, gef_seq_R2.fastq_r36h35i11.hit) and several other files.

Note that we map both R1 and R2 read files separately.

4. Paired end Gef-seq analysis: 

Using the mapped data, carry out the Gef-seq analysis using the following command.
% pair_map -job r36h35i11 -pbo 0 -pbt 2018418 -pbs 0.193873 -rgb genome.gb genome.fasta gef_seq_R1.fastq gef_seq_R2.fastq
-pbo, -pbt, -pbs options specifies the positions of ori, ter, and the scaling facter for the adjustment of ORI-TER bias, and -rgb option specifies the genbank genome anotation file.
The last four arguments specify the reference genbank file, reference fasta file, R1 read file, and R2 read file, respectively.
This will provide the following outputs. 

READ_R1.fastq_i11r36h35.peb : mapping depth in peb format (for mapsview)
READ_R1.fastq_i11r36h35.bed : mapping depth in bed format
READ_R1.fastq_i11r36h35_e.peb : edge intensity in peb format (for mapsview)
READ_R1.fastq_i11r36h35_fw_e.peb : 5' edge intensity in peb format (for mapsview)
READ_R1.fastq_i11r36h35_bw_e.peb : 3' edge intensity in peb format (for mapsview)
READ_R1.fastq_i11r36h35_edgp.fasta : list of areas between edges above threshold intensity
READ_R1.fastq_i11r36h35_rbox.fasta : list of areas specified by more than threshold number(default:10) read pair segments
READ_R1.fastq_i11r36h35.ps : Graphic representation of the analysis results in postscript format 

The postscript(.ps) file can be viewed using Preview application on Macintosh, or can be converted using ps2pdf command on Unix system.
The description of format for other files (peb and bed) are found in Gef-seq description page. 
To visualize the result, use mapsview application (mapsview_b01.app). 
mapsview_tutorial.pdf includes a brief description of how to use it.

maps programs are distributed "as is", free of charge, and without warranty of any kind.
We appreciate it if you acknowledge use of maps in any reports or publications.

References:

O. Chumsakul, D. P. Anantsri, T. Quirke, T. Oshima, K. Nakamura, S. Ishikawa, and M. M. Nakano. 
"Genome-wide Identification of ResD, NsrR, and Fur Binding in Bacillus subtilis during Anaerobic Fermentative Growth by in vivo Footprinting." 
Journal of Bacteriology, (2017). 

Onuma Chumsakul, K. Nakamura, T. Kurata, T. Sakamoto, J. L. Hobman, N. Ogasawara, T. Oshima, S. Ishikawa. 
"High-resolution mapping of in vivo genomic transcription factor binding sites using in situ DNAse footprinting and high throuput sequencing"
DNA Research, 20:325-338 (2013).

K.Nakamura,T.Oshima,T.Morimoto,S.Ikeda,H.Yoshikawa,Y.Shiwa,S.Ishikawa,M.C.Linak,A.Hirai,H.Takahashi,Md.Altaf-Ul-Amin,N.Ogasawara and S.Kanaya,
"Sequence-specific error profile of Illumina sequencers",
Nucleic Acids Research, 39, Pp. e90 (2011).

Copyright 2017. Kensuke Nakamura, All rights reserved.
