
# TSD need user to install BWA and build the reference index firstly. 
# Please replace the -G with your own one

perl LongAssembly.pl -d . -s input.fq -G /mnt/work/genome/hg38.fa -i virus.fa

