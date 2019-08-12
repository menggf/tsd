# tsd
Target sequencing discovery using long reads
# usage:

 Command: perl LongAssembly.pl -d output_dir -s seq.fq -l min_fragment_length
                  -i insert_seq.fa -G genome_ref_prefix -p cores -f overlap_score
                  -r min_reads_num -R ranking_type -o min_merge_overlap -e path_bwa
                  -m -h
				  
    LongAssembly.pl is the interface of TSD to identify the complex genomic structure 
	of targeted sequence using long reads, e.g. PacBio sequencing data.
    It is especially designed for the sequences with complex genomic structure, 
	e.g. the virus integrated sequences. 
	
	Before usage, users need to install (1) BWA and build the genome reference index.
	and (2) R for drawing the plots.
	Only two parameters are mandatory: -s and -G.
	
    Usage:
        -d: the output directory (default: .)
        -s: the sequencing files in fastq format;
        -G: the bwa genome reference prefix;
        -e: the path of bwa
        -l: the minimum fragment length for bwa alignment (default: 200)
        -f: the minimum score for fragment similarity (default: 0.7)
        -i: the extra sequences, e.g. the virus sequences integrated in the genome. (Optional but recommend)
        -p: the threads number for parallel computation
        -r: the minimum reads number for final output
        -R: the ranking type, either "length" or "fragment" (default: fragment)
        -o: the minimum overlapped fragment for merging the reads
        -m: merge the reads or not?
        -h: helps
	
# Output：

    The output of TSD is many files stored in output directory (set by "-d"). Important oness include:
    (1) output_dir/result.txt: the fragment location of long reads;
    (2) output_dir/report.txt: the genomic strcuture of SVs;
    (3) output_dir/seq_pic/*.pdf: the reads plot for SVs. One *.pdf is one SV. In each plot, the complex SVs
        are displayed as a line of fragments. 
    
# Note:

   TSD support continual analysis. That is to say, if some errors happened during previous analysis, TSD can
   restart the analysis by automaticsly detecting the temporary output files. For example, if TSD find \*.sam
   files, TSD will not run bwa alignment again. However, if bwa alignment failed in previous analysis, the 
   imcomplete\*.sam files in output_directory should be deleted, or to delete the whole output_directory!!! 

# Q & A：
    1. When to use TSD?
     A: When (a) the studied genome has complex recombination or rearrangement, or the genome is injected with 
     transgene vectors or infected by viruses (e.g. HBV), the host genome usually carries complex structure,
     e.g. multiple rearrangements in a single integration, which makes it impossible to uncover the complex 
     organization stucture using NGS or read it out directly using long reads (e.g. PacBio).  (b) Long reads 
     is sequenced (e.g. PacBio Platform), TSD can be used to identify the composition of the complex structure.
     One example is displayed in Figure 1 in the submitted paper (link), where the HBV viruses have complex 
     rearrangement before their integration in the human genome. In this example, the integrated HBV includes
     6 HBV fragments, spanning about 3000 bp. In PacBio sequencing data, no long read covers the whole region. 
     TSD recovers the HBV rearrangement profile by assemblying multiple PacBio reads. 
     
    2. What is a DNA fragment?
     A: A DNA fragment is a a piece of DNA, which can be unique mapped to the genome or known DNA, e.g. 
     part of virus genome. The DNA fragments are usually generated during the DNA rearrangement. The 
     complex stuctural variants are composed of the combination of DNA fragments.
     
    3. Error: cannot find bwa in alignment
     A: Before usage, bwa should be installed and its location added into the $PATH variable in Linux system.
    
    4. "-i" option is necessay for TSD?
     A: "-i" option is an optioinal setting, to specifiy the targeted sequences. TSD is designed for targeted 
     sequence discovery which allow identifying the genomic structure of targeted sequence. Therefore, "-i" 
     is highly recommended for TSD and it has many benefits: (1) to reduce the analysis time comsumption; 
     (2) less-redundent output. (3) good visualization to the output results.
    
    5.How to mask the homologous regions？
     A: The targeted sequences should firstly be checked using blast. For low repeat regions, e.g. SINEs, 
     Repeatmasker can help to mask such regions. If transgene vector carries homogous sequeces, the bases
     of homogous regions can be replaced with "N" based on their blast results.
     
    6. How to interprete the outputs?
     A: One demo interpretation is avialable in sampleData/*. Please read it.
     
    7. Who to contact for any problem? 
     A: Send email to Guofeng Meng(menggf@gmail.com)
