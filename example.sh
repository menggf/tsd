
# TSD need user to install BWA and build the reference index firstly. 
# Please replace the -G with your own one !!!

perl LongAssembly.pl -d . -s sampleData/input.fastq -G /mnt/work/genome/hg38.fa -i sampleData/hbv.fa

