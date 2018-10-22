# tsd
Target sequencing discovery using long reads
# usage:

 Command: perl LongAssembly.pl -d output_dir -s seq.fq -l min_fragment_length
                  -i insert_seq.fa -G genome_ref_prefix -p cores -f overlap_score
                  -r min_reads_num -R ranking_type -o min_merge_overlap -e path_bwa
                  -m -h
				  
    LongAssembly is a computational tool to identify the complex genomic structure 
	of target sequencing using long reads, e.g. PacBio sequencing data.
    It is especially designed for the sequences with complex genetic structure, 
	e.g. the virus integrated sequences. 
	
	Before usage, users need to install BWA and build the genome reference index.
	Only two parameters are mandatory: -s and -G.
	
    Usage:
        -d: the output directory (default: .)
        -s: the sequencing files in fastq format;
        -G: the bwa genome reference prefix;
        -e: the path of bwa
        -l: the minimum fragment length for bwa alignment (default: 200)
        -f: the minimum score for fragment similarity (default: 0.7)
        -i: the extra sequences, e.g. the virus sequences integrated in the genome. (Optional)
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
    (3) output_dir/seq_pic/*.pdf: the reads plot for SVs. One *.pdf is one SV;
    
# Q & A：
    1. When to use TSD?
     A: When (a) the studied genome (e.g. Hela genome), transgene vector or virus (e.g. HBV) has complex recombination or rearrangement.  (b) Long reads is sequenced (e.g. PacBio Platform), TSD can be used to identify the composition of the complex structure. 
     
    2. What is a DNA fragment?
     A: A DNA fragment is a a piece of DNA, which can be unique mapped to the genome or known DNA, e.g. part of virus genome. The complex stuctural variants are composed of the combination of DNA fragments.
     
    3. Who to contact for any problem? 
     A: Send email to Guofeng Meng(menggf@gmail.com)
