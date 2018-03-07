#!/usr/bin/perl -w
use strict;

my $output_dir=shift @ARGV; #output directory
my $seq_file=shift @ARGV; # a fastq file
my $genome_ref=shift @ARGV; # bwa genome reference
my $cores=shift @ARGV; # thread number
my $min_fragment_length=shift @ARGV; # the minimum fragment length for alignment
my $insert_seq=shift @ARGV; # target sequencing

srand();
if(!-e "$output_dir/temp_files"){
	mkdir("$output_dir/temp_files");
}
if(!-e "$output_dir/seq_pic"){
	mkdir("$output_dir/seq_pic");
}
if(!-e "$output_dir/insert_ref" and $insert_seq ne "na"){
	mkdir("$output_dir/insert_ref");
}

if(!-e "$output_dir/aln_genome.sam"){
  	my $cmd="bwa mem -t $cores -x pacbio -Y -v 1 -o ${output_dir}/aln_genome.sam $genome_ref $seq_file";
  	system($cmd);
}

my $insert_ref="$output_dir/insert_ref/insert";  
if($insert_seq ne "na"){
  #building the insert reference
  if(!-e "${insert_ref}.bwt"){
  	my $cmd="bwa index -p $insert_ref $insert_seq";
  	system($cmd);
  }
  if(!-e "$output_dir/aln_insert.sam"){
  	my $cmd="bwa mem -t $cores -x pacbio -Y -v 1 -o ${output_dir}/aln_insert.sam $insert_ref $seq_file";
  	system($cmd);
  }
}

my %all1;
my %all2;
my %range;
my %ids;
open GENOME, "${output_dir}/aln_genome.sam" or die;
while(my $str=<GENOME>){
	next if($str=~/^@/);
	my @str=split(/\t/,$str);
	next if($str[2] eq "*"  or $str[1] > 128);
	$all1{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
	$ids{$str[0]}=1;
  last;
}
while(my $str=<GENOME>){
	my @str=split(/\t/,$str);
	next if($str[2] eq "*"  or $str[1] > 128);
	$all1{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
	$ids{$str[0]}=1;
}
close GENOME;

if($insert_seq ne "na"){
	open INSERT, "${output_dir}/aln_insert.sam" or die;
	while(my $str=<INSERT>){
		next if($str=~/^@/);
		my @str=split(/\t/,$str);
		next if($str[2] eq "*" or $str[1] > 128);
		$all2{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
		$ids{$str[0]}=1;
		last;
	}
	while(my $str=<INSERT>){
		my @str=split(/\t/,$str);
		next if($str[2] eq "*" or $str[1] > 128);
		$all2{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
		$ids{$str[0]}=1;
	}
	close INSERT;
}
print "Finishing reading alignment\n";

my @ids=keys(%ids);
undef %ids;

# 0id / 1tag / 2chr / 3pos / 4cigar / 5seq
open(TEMP,">$output_dir/temp_files/temp_seq_0.fa") or die;
foreach my $id (@ids){
	if(defined($all1{$id}) and !defined($all2{$id})){
		#print "1\t$id\n";
		my @str1=split(/\t/,$all1{$id});
		my $strand1=1;
		#$strand1=-1 if($str1[1]==16);
		my ($lefts, $match_len,$match_rlen,$rights)=&cigar_len($str1[4]);
		my $from=$str1[3];
		my $to=$str1[3] + $match_len;
		$range{$id}="$str1[2]\t$strand1\t$from\t$to"; 
		if($lefts > $min_fragment_length){
			my $seq_left=substr($str1[5], 0, $lefts);
			print TEMP ">$str1[0]#left\n$seq_left\n";
		}
		if($rights > $min_fragment_length){
			my $seq_right=substr($str1[5], $lefts + $match_rlen-1 ,$rights);
			print TEMP ">$str1[0]#right\n$seq_right\n";
		}    
	}
	elsif(!defined($all1{$id}) and defined($all2{$id})){
		#print "2\t$id\n";
		my @str2=split(/\t/,$all2{$id});
		my $strand2=1;
		#$strand2=-1 if($str2[1]==16);
		my ($lefts, $match_len,$match_rlen,$rights)=&cigar_len($str2[4]);
		#print "$lefts\t$match_len\t$match_rlen\t$rights\t$str2[4]\t",length($str2[5]),"\n";
		my $from=$str2[3];
		my $to=$str2[3] + $match_len;
		$range{$id}="$str2[2]\t$strand2\t$from\t$to"; 
		if($lefts > $min_fragment_length){
			my $seq_left=substr($str2[5], 0, $lefts);
			print TEMP ">$str2[0]#left\n$seq_left\n";
		}
		if($rights > $min_fragment_length){
			my $seq_right=substr($str2[5], $lefts + $match_rlen-1 ,$rights);
			print TEMP ">$str2[0]#right\n$seq_right\n";
		}   
	}
	elsif(defined($all1{$id}) and defined($all2{$id})){
		#print "3\t$id\n";
		my @str1=split(/\t/,$all1{$id});
		my @str2=split(/\t/,$all2{$id});
		my $strand2=1;
		$strand2=-1 if($str2[1]==16);
		my $strand1=1;
		$strand1=-1 if($str1[1]==16);
		if($strand1==-1 and $strand2==-1){
			$strand1=1;
			$strand2=1;
		}
		elsif($strand1==-1 and $strand2==1){
			$strand1=1;
			$strand2=-1;
		}
		my ($lefts1, $match_len1,$match_rlen1,$rights1)=cigar_len($str1[4]);
		my ($lefts2, $match_len2,$match_rlen2,$rights2)=cigar_len($str2[4]);
		my $from2=$str2[3];
		my $to2=$str2[3] + $match_len2;
		my $from1=$str1[3];
		my $to1=$str1[3] + $match_len1;	
		if($strand2 ==-1){
			my $tmp=$lefts2;
			$lefts2=$rights2;
			$rights2=$tmp;
		}
		if($lefts1 <= $lefts2){
			if($lefts1 > $min_fragment_length){
				my $seq_right=substr($str1[5], 0 ,$lefts1);
				print TEMP ">$str1[0]#left\n$seq_right\n";
			}
			my $gap=$lefts2 -$lefts1 - $match_rlen1;
			$range{$str2[0]}="$str1[2]\t$strand1\t$from1\t$to1\t1000\t$gap\t$str2[2]\t$strand2\t$from2\t$to2"; 
			if( $gap > $min_fragment_length){
				my $seq_mid=substr($str1[5], $lefts1 + $match_rlen1-1, $gap);
				print TEMP ">$str2[0]#1000\n$seq_mid\n";
			}
			if($rights2 > $min_fragment_length){
				my $seq_right=substr($str1[5], $lefts2 +$match_rlen2-1 ,$rights2);
				print TEMP ">$str1[0]#right\n$seq_right\n";
			}
		}
		if($lefts1 > $lefts2){
			if($lefts2 > $min_fragment_length){
				my $seq_right=substr($str1[5], 0 ,$lefts2);
				print TEMP ">$str1[0]#left\n$seq_right\n";
			}
			my $gap=$lefts1 -$lefts2- $match_rlen2;
			$range{$str2[0]}="$str2[2]\t$strand2\t$from2\t$to2\t1000\t$gap\t$str1[2]\t$strand1\t$from1\t$to1"; 
			if( $gap > $min_fragment_length){
				my $seq_mid=substr($str1[5], $lefts2 + $match_rlen2-1 , $gap);
				print TEMP ">$str2[0]#1000\n$seq_mid\n";
			}
			if( $rights1 > $min_fragment_length){
				my $seq_right=substr($str1[5], $lefts1 +$match_rlen1-1 ,$rights1);
				print TEMP ">$str1[0]#right\n$seq_right\n";
			}
		}
	}
}
close(TEMP);
undef %all2;
undef %all1;

my $pp=0;
my @rsd;
for(my $i=0;$i< 1000;$i=$i+6){
  push @rsd, $i+4;
}
while(-s "${output_dir}/temp_files/temp_seq_$pp.fa" !=0){
	print "Begin $pp-round alignment..\n";
	my $cmd1="bwa mem -t $cores -x pacbio -Y -v 1 -o ${output_dir}/temp_files/temp_genome_$pp.sam $genome_ref ${output_dir}/temp_files/temp_seq_$pp.fa";
	system($cmd1);
	if($insert_seq  ne "na"){ 
		my $cmd2="bwa mem -t $cores -x pacbio -Y -v 1 -o ${output_dir}/temp_files/temp_insert_$pp.sam $insert_ref ${output_dir}/temp_files/temp_seq_$pp.fa";
		system($cmd2);
	}
	print "Finishing $pp-round alignment.\n";
	open GENOME, "${output_dir}/temp_files/temp_genome_$pp.sam" or die;
	my %all1;
	my %all2;
	my %ids;
	while(my $str=<GENOME>){
		next if($str=~/^@/);
		my @str=split(/\t/,$str);
		next if($str[2] eq "*" or $str[1] > 128);
		$all1{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
		$ids{$str[0]}=1;
	}
	close GENOME;
	if($insert_seq ne "na"){ 
		open INSERT, "${output_dir}/temp_files/temp_insert_$pp.sam" or die;
		while(my $str=<INSERT>){
			next if($str=~/^@/);
			my @str=split(/\t/,$str);
			next if($str[2] eq "*" or $str[1] > 128);
			$all2{$str[0]}="$str[0]\t$str[1]\t$str[2]\t$str[3]\t$str[5]\t$str[9]";
			$ids{$str[0]}=1;
		}
		close INSERT;
	}
	my @ids=keys(%ids);
	undef %ids;
	$pp++;
	open(TEMP,">${output_dir}/temp_seq_$pp.fa") or die;
	foreach my $str (@ids){
		die "$str" if($str !~/(\S+)#(\S+)$/);
		my $id=$1;
		my $tag=$2;	
		my @rr=split(/\t/,$range{$id});
		my %has;
		if(@rr > 4){
			for(my $i=4;$i< @rr;$i=$i+6){
				$has{$rr[$i]}=1;
			}
		}
		if(defined($all1{$str}) and !defined($all2{$str})) {
			my @str1=split(/\t/,$all1{$str});
			#my @str2=split(/\t/,$all2{$str});
			my $strand=1;
			$strand=-1 if($str1[1]==16);
			my ($lefts, $match_len, $match_rlen, $rights)=&cigar_len($str1[4]);
			my $from=$str1[3];
			my $to=$str1[3] + $match_len;
			my $sq=$str1[5];
			if($strand==-1){
				my $tmp=$lefts;
				$lefts=$rights;
				$rights=$tmp;
				$sq=&revs($sq);
			}
			if($tag eq "left"){
				my $cc=int(rand(1000));
				while(defined($has{$cc})){
					$cc=int(rand(1000));
				}
				$has{$cc}=1;
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					print TEMP ">$id#left\n$seq_left\n";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen-1 ,$rights);
					print TEMP ">$id#$cc\n$seq_right\n";
				}
				$range{$id}="$str1[2]\t$strand\t$from\t$to\t$cc\t$rights\t$range{$id}";
			}
			elsif($tag eq "right"){
				my $cc=int(rand(1000));
				while(defined($has{$cc})){
					$cc=int(rand(1000));
				}
				$has{$cc}=1;
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					print TEMP ">$id#$cc\n$seq_left\n";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen-1 ,$rights);
					print TEMP ">$id#right\n$seq_right\n";
				}
				$range{$id}="$range{$id}\t$cc\t$lefts\t$str1[2]\t$strand\t$from\t$to";
			}
			else{
				my $temp1="0\t$lefts";
				my $temp2="0\t$rights";
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					my $cc=int(rand(1000));
					while(defined($has{$cc})){
						$cc=int(rand(1000));
					}
					$has{$cc}=1;
					print TEMP ">$id#$cc\n$seq_left\n";
					$temp1="$cc\t$lefts";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen ,$rights);
					my $cc=int(rand(1000));
					while(defined($has{$cc})){
						$cc=int(rand(1000));
					}
					$has{$cc}=1;
					print TEMP ">$id#$cc\n$seq_right\n";
					$temp2="$cc\t$rights";
				}
				my $wh=-1;
				foreach my $i (@rsd){
					if($rr[$i]==$tag){
						$wh=$i;
						last;
					}
				}
				next if($wh==-1);
				my $ll=join("\t", @rr[0..($wh-1)]);
				my $rr=join("\t", @rr[($wh+2)..(scalar(@rr)-1)]);
				$range{$id}="$ll\t$temp1\t$str1[2]\t$strand\t$from\t$to\t$temp2\t$rr";
			}
		}
		if(!defined($all1{$str}) and defined($all2{$str})) {
			#my @str1=split(/\t/,$all1{$str});
			my @str2=split(/\t/,$all2{$str});
			my $strand=1;
			$strand=-1 if($str2[1]==16);
			my ($lefts, $match_len,$match_rlen, $rights)=cigar_len($str2[4]);
			my $from=$str2[3];
			my $to=$str2[3] + $match_len;
			my $sq=$str2[5];
			if($strand==-1){
				my $tmp=$lefts;
				$lefts=$rights;
				$rights=$tmp;
				$sq=revs($sq);
			}
			if($tag eq "left"){
				my $cc=int(rand(1000));
				while(defined($has{$cc})){
					$cc=int(rand(1000));
				}
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					print TEMP ">$id#left\n$seq_left\n";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen-1,$rights);
					print TEMP ">$id#$cc\n$seq_right\n";
				}
				$range{$id}="$str2[2]\t$strand\t$from\t$to\t$cc\t$rights\t$range{$id}";
			}
			elsif($tag eq "right"){
				my $cc=int(rand(1000));
				while(defined($has{$cc})){
					$cc=int(rand(1000));
				}
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					print TEMP ">$id#$cc\n$seq_left\n";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen-1 ,$rights);
					print TEMP ">$id#right\n$seq_right\n";
				}
				$range{$id}="$range{$id}\t$cc\t$lefts\t$str2[2]\t$strand\t$from\t$to";
			}
			else{
				my $temp1="0\t$lefts";
				my $temp2="0\t$rights";
				if($lefts > $min_fragment_length){
					my $seq_left=substr($sq, 0, $lefts);
					my $cc=int(rand(1000));
					while(defined($has{$cc})){
						$cc=int(rand(1000));
					}
					print TEMP ">$id#$cc\n$seq_left\n";
					$temp1="$cc\t$lefts";
				}
				if($rights > $min_fragment_length){
					my $seq_right=substr($sq, $lefts + $match_rlen ,$rights);
					my $cc=int(rand(1000));
					while(defined($has{$cc})){
						$cc=int(rand(1000));
					}
					print TEMP ">$id#$cc\n$seq_right\n";
					$temp2="$cc\t$rights";
				}
				my $wh=-1;
				foreach my $i (@rsd){
					if($rr[$i] == $tag){
						$wh=$i;
						last;
					}
				}
				next if($wh==-1);
				my $ll=join("\t", @rr[0..($wh-1)]);
				my $rr=join("\t", @rr[($wh+2)..(scalar(@rr)-1)]);
				$range{$id}="$ll\t$temp1\t$str2[2]\t$strand\t$from\t$to\t$temp2\t$rr";
			}
		}
		if(defined($all1{$str}) and defined($all2{$str})) {
			my @str1=split(/\t/,$all1{$str});
			my @str2=split(/\t/,$all2{$str});
			my $strand2=1;
			$strand2=-1 if($str2[1]==16);
			my $strand1=1;
			$strand1=-1 if($str1[1]==16);
			my ($lefts1, $match_len1,$match_rlen1,$rights1)=cigar_len($str1[4]);
			my ($lefts2, $match_len2,$match_rlen2,$rights2)=cigar_len($str2[4]);
			my $from2=$str2[3];
			my $to2=$str2[3] + $match_len2;
			my $from1=$str1[3];
			my $to1=$str1[3] + $match_len1;
			my $sq=$str1[5];
			if($strand1==-1){
				my $tmp=$lefts1;
				$lefts1=$rights1;
				$rights1=$tmp;
				$sq=&revs($sq);
			}
			if($strand2==-1){
				my $tmp=$lefts2;
				$lefts2=$rights2;
				$rights2=$tmp;
			}
			
			if($tag eq "left"){
				if($lefts1 <= $lefts2){
					my $temp1="0\t$lefts1";
					my $temp3="0\t$rights2";
					if($lefts1 > $min_fragment_length){
						my $seq_right=substr($sq, 0 ,$lefts1);
						print TEMP ">$id#left\n$seq_right\n";
					}
					
					if($rights2 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts2 +$match_rlen2-1 ,$rights2);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp3="$cc\t$rights2";
						print TEMP ">$id#$cc\n$seq_right\n";
					}
					my $gap=$lefts2 -$lefts1 - $match_rlen1;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts1 + $match_rlen1-1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					$range{$id}="$str1[2]\t$strand1\t$from1\t$to1\t$temp2\t$str2[2]\t$strand2\t$from2\t$to2\t$temp3\t$range{$id}";
				}
				if($lefts1 > $lefts2){
					my $temp1="0\t$lefts2";
					my $temp3="0\t$rights1";
					if($lefts2 > $min_fragment_length){
						my $seq_left=substr($sq, 0 ,$lefts2);
						print TEMP ">$id#left\n$seq_left\n";
					}
					if($rights1 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts1 +$match_rlen1-1 ,$rights1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp3="$cc\t$rights1";
						print TEMP ">$id#$cc\n$seq_right\n";
					}
					my $gap=$lefts1 -$lefts2 - $match_rlen2;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts2 + $match_rlen2-1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					$range{$id}="$str2[2]\t$strand2\t$from2\t$to2\t$temp2\t$str1[2]\t$strand1\t$from1\t$to1\t$temp3\t$range{$id}";
				}
			}
			elsif($tag eq "right"){
				if($lefts1 <= $lefts2){
					my $temp1="0\t$lefts1";
					my $temp3="0\t$rights2";
					if($lefts1 > $min_fragment_length){
						my $seq_right=substr($sq, 0 ,$lefts1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						print TEMP ">$id#$cc\n$seq_right\n";
						$temp1="$cc\t$lefts1";
					}
					
					if($rights2 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts2 +$match_rlen2-1 ,$rights2);
						print TEMP ">$id#right\n$seq_right\n";
					}
					my $gap=$lefts2 -$lefts1-$match_rlen1;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts1 +$match_rlen1 -1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					$range{$id}="$range{$id}\t$temp1\t$str1[2]\t$strand1\t$from1\t$to1\t$temp2\t$str2[2]\t$strand2\t$from2\t$to2";
				}
				if($lefts1 > $lefts2){
					my $temp1="0\t$lefts2";
					my $temp3="0\t$rights1";
					if($lefts2 > $min_fragment_length){
						my $seq_right=substr($sq, 0 ,$lefts2);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						print TEMP ">$id#$cc\n$seq_right\n";
						$temp1="$cc\t$lefts2";
					}				
					if($rights1 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts1 +$match_rlen1-1 ,$rights1);
						print TEMP ">$id#right\n$seq_right\n";
					}
					my $gap=$lefts1 -$lefts2 -  $match_rlen2;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts2 +$match_rlen2 -1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					$range{$id}="$range{$id}\t$temp1\t$str2[2]\t$strand2\t$from2\t$to2\t$temp2\t$str1[2]\t$strand1\t$from1\t$to1";
				}
			}
			else{
				if($lefts1 <= $lefts2){
					my $temp1="0\t$lefts1";
					my $temp3="0\t$rights2";
					if($lefts1 > $min_fragment_length){
						my $seq_right=substr($sq, 0 ,$lefts1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						print TEMP ">$id#$cc\n$seq_right\n";
						$temp1="$cc\t$lefts1";
					}
					if($rights2 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts2 +$match_rlen2-1 ,$rights2);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						print TEMP ">$id#right\n$seq_right\n";
						$temp3="$cc\t$rights2";
					}
					my $gap=$lefts2 -$lefts1-$match_rlen1;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts1 +$match_rlen1 -1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					my $wh=-1;
					foreach my $i (@rsd){
						if($rr[$i]==$tag){
							$wh=$i;
							last;
						}
					}
					next if($wh==-1);
					my $ll=join("\t", @rr[0..($wh-1)]);
					my $rr=join("\t", @rr[($wh+2)..(scalar(@rr)-1)]);
					$range{$id}="$ll\t$temp1\t$str1[2]\t$strand1\t$from1\t$to1\t$temp2\t$str2[2]\t$strand2\t$from2\t$to2\t$temp3\t$rr";
				}
				if($lefts1 > $lefts2){
					my $temp1="0\t$lefts2";
					my $temp3="0\t$rights1";
					if($lefts2 > $min_fragment_length){
						my $seq_right=substr($sq, 0 ,$lefts2);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						print TEMP ">$id#$cc\n$seq_right\n";
						$temp1="$cc\t$lefts2";
					}				
					if($rights1 > $min_fragment_length){
						my $seq_right=substr($sq, $lefts1 +$match_rlen1-1 ,$rights1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp3="$cc\t$rights1";
						print TEMP ">$id#$cc\n$seq_right\n";
						print TEMP ">$id#right\n$seq_right\n";
					}
					my $gap=$lefts1 -$lefts2-$match_rlen2;
					my $temp2="0\t$gap";
					if( $gap > $min_fragment_length){
						my $seq_mid=substr($sq, $lefts2 +$match_rlen2 -1, $gap+1);
						my $cc=int(rand(1000));
						while(defined($has{$cc})){
							$cc=int(rand(1000));
						}
						$has{$cc}=1;
						$temp2="$cc\t$gap";
						print TEMP ">$id#$cc\n$seq_mid\n";
					}
					my $wh=-1;
					foreach my $i (@rsd){
						if($rr[$i]==$tag){
							$wh=$i;
							last;
						}
					}
					next if($wh==-1);
					my $ll=join("\t", @rr[0..($wh-1)]);
					my $rr=join("\t", @rr[($wh+2)..(scalar(@rr)-1)]);
				
					$range{$id}="$ll\t$temp1\t$str2[2]\t$strand2\t$from2\t$to2\t$temp2\t$str1[2]\t$strand1\t$from1\t$to1\t$temp3\t$rr";
				}
			}
		}
	}
	close(TEMP);
}

open OUT,">$output_dir/result.txt" or die;
foreach my $id (keys %range){
	print OUT "$id\t$range{$id}\n";
}
close OUT;

exit(0);

sub cigar_len{
	my $s=shift @_;
	my $ll=0;
	my $rr=0;
	my $mm=0;
	my $rmm=0;
	$ll=$1 if($s=~/^(\d+)[HS]/);
	$rr=$1 if($s=~/(\d+)[HS]$/);
	while($s=~/(\d+)(\w)/g){
		if($2 eq "D"  or $2 eq "M" or $2 eq "N" or $2 eq "X" or $2 eq "="){
			$mm += $1;
		}
		if($2 eq "M" or $2 eq "I" or $2 eq "X" or $2 eq "="){
			$rmm += $1;
		}
	}
	return($ll, $mm, $rmm, $rr);
}
sub revs{
	my $s=shift @_;
    $s=reverse($s);
	$s=~tr/ATCGatcg/TAGCtagc/;
	return($s);
}
