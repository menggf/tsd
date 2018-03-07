#!/usr/bin/perl -w
use strict;
use align;

my $output_dir = shift @ARGV; 
my $RANK_TYPE = shift @ARGV; 
my $min_reads = shift @ARGV; 
my $min_fragment_length = shift @ARGV; #the fragment length for further split
my $min_merge_overlap = shift @ARGV;  #
my $merge = shift @ARGV; #merge the reads or not
my $cutoff_similarity = shift @ARGV;  #overlap scores cutoff

open DATA,"$output_dir/result.txt" or die; # read the fragments
my %sites; 
my %all; 
my %len; 
while(my $str=<DATA>){
	next if($str!~/insert/); # !!!!!!!!!!!!!!!!!!!
	chomp $str;
	my @str=split(/\t/,$str);
	next if(@str <= 6);
	my $id=shift @str;
	$all{$id}=join("\t", @str);
	my $nn=scalar(@str);
	$len{$id}=0;
	for(my $i=0;$i<$nn; $i=$i+6){
		$sites{"$str[$i]\t$str[$i+2]\t$str[$i+3]"}="";
		if($RANK_TYPE eq "length"){
			$len{$id}+=($str[$i+3]-$str[$i+2]);
			if($i+5 <$nn){
				$len{$id}+=$str[$i+5];
			}
		}
		else{
			$len{$id}=($nn+2)/6;
		}
	}
}
close DATA;

system("Rscript cluster_fragment.r $output_dir"); # temp

open DATA,"$output_dir/new_sites.txt" or die;
my %ann;
my %ann_chr;
my %ann_len;
while(my $str=<DATA>){
	chomp $str;
	my @str=split(/\t/,$str);
	$sites{"$str[0]\t$str[1]\t$str[2]"}=$str[3];
	$ann{$str[3]}="$str[0]\t$str[4]\t$str[5]";
	$ann_chr{$str[3]}=$str[0];
	$ann_len{$str[3]}=$str[5]-$str[4]+1;
}
close DATA;

my %dd1; # score 1
my %dd2; # score 2
my @ss=keys %ann;
my %freq;
for(my $i=0;$i<@ss;$i++){
	my @t1=split(/\t/,$ann{$ss[$i]});
	my $from1=$t1[1];
	my $to1=$t1[2];
	next if($to1 - $from1 < 100);
	for(my $j=$i;$j<@ss;$j++){
		next if($i==$j);
		next if($ann_chr{$ss[$i]} ne $ann_chr{$ss[$j]});
		my @t2=split(/\t/,$ann{$ss[$j]});
		my $from2=$t2[1];
		my $to2=$t2[2];
		next if($to2 - $from2 < 100);
		next if($to1 < $from2 or $to2 < $from1);
		
		my $overlap=0;
		if($from1 < $from2 ){
			$overlap=min($to1, $to2)-$from2+1;
		}
		elsif($from1 >= $from2){
			$overlap= min($to1,$to2)-$from1 + 1;
		}
		my $r1=$overlap/min($to2-$from2 +1, $to1-$from1 +1);
		my $r2=$overlap/max($to2-$from2 +1, $to1-$from1 +1);
		#print "$from1\t$to1\t$from2\t$to2\t$r1\t$r2\n";
		next if($r1 < $cutoff_similarity);
		$dd1{"$ss[$i]\t$ss[$j]"}=$r1;
		$dd2{"$ss[$i]\t$ss[$j]"}=$r2;
		if($r2 > $cutoff_similarity){
			$freq{$ss[$i]}++;
			$freq{$ss[$j]}++;
		}
	}
}



my %pattern; #reformat into site string
my %canseed; #can be used as seed or not
foreach my $id (keys %all){
	#next if($id ne "m54079_171127_175044/34145083/26887_32368");
	my $str=$all{$id};
	$canseed{$id}=1;
	next if($str!~/insert/); #!!!!!!!!!!!!!!!!!!!!!!!!
	my @str=split(/\t/,$str);
	my $nn=scalar(@str);
	my $s="";
	my @x;
	my $fr=-1;
	
	for(my $i=0;$i<$nn; $i=$i+6){
		my $a=$sites{"$str[$i]\t$str[$i+2]\t$str[$i+3]"};
		if($s eq ""){
			$s="$a\t$str[$i+1]";
		}
		else{
			$s="$s\t$a\t$str[$i+1]";
		}
		if($i+5 < $nn){
			$s="$s\t$str[$i+5]";
			$canseed{$id}=0 if($str[$i+5] > 400);
		}
	}
	# find the repeated sequencing regions
	#my @s=split(/\t/,$s);
	#$nn=(scalar(@s)+1)/3;
	#my @matrix;
	#$matrix[0][0]  = 0;
	#for ( my $j = 1 ; $j <= $nn ; $j++ ) {
    #    $matrix[0][$j]   = 0;
    #}
    #for ( my $i = 1 ; $i <= $nn ; $i++ ) {
    #    $matrix[$i][0]   = 0;
    #}
	#my $wh_x=0;
	#my $wh_y=0;
	#my $max_score=0;
	#for(my $q=1; $q <= $nn; $q++){
	#	for(my $p=1; $p <= $nn; $p++){
	#		my ($u, $d, $m)=(0,0,0);
	#		my $a=($q-1)*3;
	#		my $b=($p-1)*3;
	#		if($s[$a] eq $s[$b] ){ #match
	#			$m=1;
	#		}
	#		else{
	#			$m=-100;
	#		}
	#		$m=0 if($q == $p);
	#		$u=100;
	#		$d=100;
	#		$matrix[$q][$p]=mymax($matrix[$q-1][$p-1] + $m, $matrix[$q][$p-1] - $d, $matrix[$q-1][$p] - $u, 0);
	#		$matrix[$q][$p]=0 if($p ==$q);
	#		if($matrix[$q][$p] > $max_score){
	#			$max_score=$matrix[$q][$p];
	#			$wh_x=$a;
	#			$wh_y=$b;
	#		}
	#	}
	#}
	#for(my $q=0; $q <= $nn; $q++){
	#	for(my $p=0; $p <= $nn; $p++){
	#		print " ",$matrix[$q][$p]
	#	}
	#	print "\n";
	#}
	#if($max_score >=4){
	#	$pattern{$id}=join("\t", @s[($wh_x-$max_score*3)..($wh_x+1)]);
	#	if($RANK_TYPE ne "length"){
	#		$len{$id}=$max_score
	#	}
	#}
	#else{
		$pattern{$id}=$s;
	#}
}
close DATA;

my %used;
my @reads=sort{$len{$b} <=> $len{$a}}(keys %len);
open OUT,">$output_dir/report.txt" or die;
my $counter=1;
for(my $i=0;$i<@reads;$i++){
	#next if($reads[$i] ne "m54079_171127_175044/34341462/3993_16767");
	next if(!$canseed{$reads[$i]});
	print "reads $i $counter\n";
	next if(defined($used{$reads[$i]}));
	my %first; # where is the location
	my $has=0; #how many supporting reads
	my @tp=split(/\t/,$pattern{$reads[$i]}); # the seed sequences
	my $np=(scalar(@tp)+1)/3;
	my @cs=@tp; # the consensus
	my $nn=(scalar(@cs)+1)/3;
	next if($nn < 2);
	my @ccs; # consensus 
	my %tmp;
	my $extended=0;
	my @extend;
	my @ext_dir;
	my @info=split(/\t/,$all{$reads[$i]});
	for(my $j=0;$j<@reads;$j++){
		#next if($reads[$j] ne "m54079_171127_175044/32571837/27563_33463");
		next if($i==$j or defined($tmp{$reads[$j]})); #re-do the match by ignore used ones
		my @pr=split(/\t/, $pattern{$reads[$j]}); #partner sequences
		my $mm=(scalar(@pr)+1)/3;
		next if($mm < 2);
		$nn=(scalar(@cs)+1)/3;
		my $output=align->new(\@cs, \@pr, \%dd1, \%dd2, \%ann_len,$min_fragment_length,$cutoff_similarity);
		if($output->{max_score} < 2 and $output->{max_score} != -1){ # reverse the partner fragments
			@pr=&revfrag(@pr);
			$output=align->new(\@cs, \@pr, \%dd1, \%dd2, \%ann_len,$min_fragment_length, $cutoff_similarity);
		}
		next if($output->{max_score} < 2);
		next if($output->{used} == 0);
		my $a=$output->{a};
		my $b=$output->{b};
		my @aa;
		my @bb;
		foreach (@$a){ push @aa, $_;}
		foreach (@$b){ push @bb, $_;}
		next if($aa[0] !=0 and $bb[0]!=0);
		my $na=scalar(@aa);
		next if($aa[$na-1] != ($nn-1) and $bb[$na-1]!=($mm-1));
		if($aa[0]==0 and $bb[0]!=0 and $output->{max_score} >= $min_merge_overlap and $merge){ #extend the head of seed
			$extended=1;
			my @pinfo=split(/\t/,$all{$reads[$j]}); #partner
			for(my $zz=3*$bb[0]-1;$zz >=0;$zz--){
				unshift @cs, $pr[$zz];
			}
			for(my $zz=6*$bb[0]-1;$zz >=0;$zz--){
				unshift @info, $pinfo[$zz];
			}			
			my $sf=$bb[0];
			if(@extend!=0){
				foreach (@extend){
					$_+=$sf;
				}
			}
			push @extend, $sf;
			push @ext_dir, "l";
			foreach (@aa){
				$_+=$sf;
			}
			for(my $zz=$sf-1;$zz >= 0;$zz--){
				unshift @aa, $zz;
				unshift @bb, $zz;
			}
			foreach my $s (@ccs){
				my @s=split(/#/, $s);
				my @as=split(/\t/,$s[0]);
				foreach(@as){
					$_+= $sf;
				}
				$s[0]=join("\t",@as);
				$s=join("#", @s);
			}
			my $as=join("\t", @aa);
			my $bs=join("\t", @bb);
			#push @ccs, "$as#$bs";
			$tmp{$reads[$j]}=1;
			$has++;
			$j=0;
		}
		elsif($aa[$na-1]==($nn-1) and $bb[$na-1]!= ($mm-1) and $output->{max_score} >= $min_merge_overlap and $merge){ #extend the seed  at tail
			$extended=1;
			my @pinfo=split(/\t/,$all{$reads[$j]}); #partner
			push @extend, $nn;
			push @ext_dir, "r";
			for(my $zz = ( $bb[$na-1] + 1 ) *3 -1;$zz < scalar(@pr); $zz++){
				push @cs, $pr[$zz];
			}
			push @info, 0;
			for(my $zz = ( $bb[$na-1] + 1 ) *6 -1; $zz < scalar(@pinfo); $zz++){
				push @info, $pinfo[$zz];
			}
			for(my $zz=$bb[$na-1]+1 ;$zz < $mm; $zz++){
				push @aa, $aa[scalar(@aa)-1]+1;
				push @bb, $zz;
			}
			my $as=join("\t", @aa);
			my $bs=join("\t", @bb);
			#push @ccs, "$as#$bs";
			$j=0;
			$tmp{$reads[$j]}=1;
			$has++;
		}
		else{
			my $as=join("\t", @aa);
			my $bs=join("\t", @bb);
			push @ccs, "$as#$bs";
			$tmp{$reads[$j]}=1;
			$has++;
		}
	}
	if($has > $min_reads){
		my %ccs;
		my @cv;
		for(my $xx=0; $xx< (scalar(@cs)+1)/3;$xx++){
			push @cv, 0;
		}
		foreach my $s (@ccs){
			my @s=split(/#/,$s);
			my @aa=split(/\t/,$s[0]);
			foreach (@aa){
				next if($_ < 0);
				$cv[$_]++;
			}
			$s=~s/#/\t/;
			$ccs{$s}++;
		}
		my $start=0;
		my $end=scalar(@cv)-1;
		for(my $p=0;$p<=$end; $p++){
			$start=$p;
			last if($cv[$p] >= 3);
		}
		for(my $p=$end;$p >= 0;$p--){
			$end=$p;
			last if($cv[$p] >= 3);
		}
		next if($end - $start < 3);
		open(my $tmp,">","$output_dir/temp_files/temp_$i.txt") or die;
		print $tmp "$reads[$i]\t",(scalar(@cs[(3*$start)..(3*$end+1)])+1)/3,"\t$start\t$end\n";
		print $tmp  join("\t", @info),"\n";
		
		my $pp=0;
		foreach my $s (sort{$ccs{$b} <=> $ccs{$a}}keys %ccs){ #counted frequency
			print $tmp "$ccs{$s}\t$s\n";
			$pp++;
			last if($pp > 100);
		}
		close $tmp;
		my $ext=join(",",@extend); #extend information
		my $ext2=join(",",@ext_dir); # extended direction
		if($extended){
			system("Rscript plot_fragment.r $output_dir/temp_files/temp_$i.txt $counter $reads[$i]/E $ext $ext2");
		}
		else{
			system("Rscript plot_fragment.r $output_dir/temp_files/temp_$i.txt $counter $reads[$i] $ext $ext2");
		}
		
		#unlink("temp_files/temp_$i.txt");
		foreach (keys %tmp){
			$used{$_}=1;
		}
		@cv=@cv[$start..$end];
		print OUT "seq$counter\t$reads[$i]\t$nn\t$has\t",mymin(@cv),"\t", coverage(@cv),"\t",join(",",@cv);
		#my @all=split(/\t/,$all{$reads[$i]});
		for(my $yy=$start*6; $yy <($end*6+4); $yy++){
			next if(($yy-4)%6==0);
			print OUT "\t",$info[$yy];
		}
		print OUT "\n";
		$counter++;
	}
}
close OUT;
exit(0);

sub min{
  my ($a, $b)=@_;
  if($a < $b){
    return $a;
  }
  else{
    return $b;
  }
}
sub max{
  my ($a, $b)=@_;
  if($a > $b){
    return $a;
  }
  else{
    return $b;
  }
}
sub mymax{
  my @x=@_;
  my $re=-1000;
  foreach (@x){
	$re=$_ if($re < $_);
  }
  return $re;
}
sub mymin{
  my @x=@_;
  my $re=10000000;
  foreach (@x){
	$re=$_ if($re > $_);
  }
  return $re;
}
sub coverage{
  my @x=@_;
  my $cc=0;
  foreach (@x){
	$cc++ if( $_ > 0);
  }
  return $cc/scalar(@x);
}
sub mean{
  my @x=@_;
  my $sum=0;
  foreach (@x){
	$sum=$_;
  }
  return($sum/scalar(@x));
}
sub revfrag{ #reverse the pattern
	my @x=@_;
	my $nn=scalar(@x);
	my @re;
	for(my $i=0;$i<$nn;$i=$i+3){
		unshift @re, -1*$x[$i+1];
		unshift @re, $x[$i];
		if($i < $nn - 2){
			unshift @re, $x[$i+2];
		}
	}
	return @re;
}
