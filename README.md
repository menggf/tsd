# tsd
Target sequencing discovery using long reads
# usage:

 Command: perl LongAssembly.pl -d output_dir -s seq.fq -l min_fragment_length
                  -i insert_seq.fa -G genome_ref_prefix -p cores -f overlap_score
                  -r min_reads_num -R ranking_type -o min_merge_overlap -e path_bwa
                  -m -h
				  
    LongAssembly is a de novo assembly tool for long reads, e.g. PacBio sequencing data.
    It is designed for the sequences with complex genetic structure, e.g. the virus 
    integrated sequences.
	
    Usage:
        -d: the output directory (default: .)
        -s: the sequencing files in fastq format
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
The output of TSD is many files stored in output directory (set by "-d"). 
	
    Contact: Guofeng Meng(menggf@gmail.com)
