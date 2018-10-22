####################################################################################
#  assmbelyCF.pm @ TSD
#  Info: Assemble the concensus fragment of PacBio reads
#  Output: The location iformation of SVs
#
#  Copyright (c) 2017- Guofeng Meng
#
#  EMAIL: menggf@gmail.com


#!/usr/bin/perl -w
use strict;
use TSD::align;
package TSD::assemblyCF;


sub new{
	my $class=shift;
	my ($output_dir, $RANK_TYPE, $min_reads, $min_fragment_length, $min_merge_overlap, $merge,$cutoff_similarity, $lib, $insert_seq)=@_;
  my $self={};
  #my $output_dir = shift @ARGV; 
  #my $RANK_TYPE = shift @ARGV; 
  #my $min_reads = shift @ARGV; 
  #my $min_fragment_length = shift @ARGV; #the fragment length for further split
  #my $min_merge_overlap = shift @ARGV;  #
  #my $merge = shift @ARGV; #merge the reads or not
  #my $cutoff_similarity = shift @ARGV;  #overlap scores cutoff
  
  open DATA,"$output_dir/result.txt" or die; # read the fragments
  open OUT,">$output_dir/sites.txt" or die;
  my %sites; 
  my %all; 
  my %len; 
  while(my $str=<DATA>){
    if($insert_seq ne "na"){
  	  next if($str!~/insert/); # !!!!!!!!!!!!!!!!!!!
    }
  	chomp $str;
  	my @str=split(/\t/,$str);
  	next if(@str <= 6);
  	my $id=shift @str;
  	$all{$id}=join("\t", @str);
  	my $nn=scalar(@str);
  	$len{$id}=0;
  	for(my $i=0;$i<$nn; $i=$i+6){
  		$sites{"$str[$i]\t$str[$i+2]\t$str[$i+3]"}="";
  		print OUT "$str[$i]\t$str[$i+2]\t$str[$i+3]\n";
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
  close OUT;
  
  
  system("Rscript $lib/TSD/cluster_fragment.r $output_dir"); # temp
  #print "ok\n";
  open DATA,"$output_dir/new_sites.txt" or die;
  my %ann_from;
  my %ann_to;
  my %ann_chr;
  my %ann_len;
  while(my $str=<DATA>){
  	chomp $str;
  	my @str=split(/\t/,$str);
  	$sites{"$str[0]\t$str[1]\t$str[2]"}=$str[3];
  	$ann_from{$str[3]} = $str[4];
  	$ann_to{$str[3]}   = $str[5];
  	$ann_chr{$str[3]}  = $str[0];
  	$ann_len{$str[3]}  = $str[5]-$str[4]+1;
  }
  close DATA;
  
  my %dd1; # score 1
  my %dd2; # score 2
  my @ss=keys %ann_chr;
  my %freq;
  for(my $i=0;$i<@ss;$i++){
  	my $from1=$ann_from{$ss[$i]};
  	my $to1=$ann_to{$ss[$i]};
  	next if($to1 - $from1 < 100);
  	for(my $j=$i;$j<@ss;$j++){
  		next if($i==$j);
  		next if($ann_chr{$ss[$i]} ne $ann_chr{$ss[$j]});
  		my $from2=$ann_from{$ss[$j]};
  		my $to2=$ann_to{$ss[$j]};
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
  foreach my $id (keys %all){
  	my $str=$all{$id};
    if($insert_seq ne "na"){
  	  next if($str!~/insert/); # !!!!!!!!!!!!!!!!!!!
    }
  	my @str=split(/\t/,$str);
  	my $nn=scalar(@str);
  	my $s="";
  	my @x;
  	my $fr=-1;
  	my $tg=1;
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
  			if($str[$i+5] > 200 or $str[$i+5] < -100){
  				$tg=0;
  				last;	
  			}
  		}
  	}
  	$pattern{$id}=$s if($tg);
  }
  close DATA;
  
  my %used; # it has been used as supporting reads,not used as seed again
  my @reads=sort{$len{$b} <=> $len{$a}}(keys %pattern); #rank the reads
  open OUT,">$output_dir/report.txt" or die;
  my $counter=1;
  for(my $i=0;$i<@reads;$i++){
  	#next if($reads[$i] ne "m54079_171127_175044/20185230/18772_29123");
  	#print "reads $i $counter\n";
  	next if(defined($used{$reads[$i]})); #skip it if it has included in previous steps
  	my %first; # where is the location
  	my $has=0; #how many supporting reads
  	my @tp=split(/\t/,$pattern{$reads[$i]}); # the seed sequences
  	my $np=(scalar(@tp)+1)/3;
  	my @cs=@tp; # the consensus
  	my $nn=(scalar(@cs)+1)/3;
  	next if($nn < 3);
  	my @ccs;
  	my %tmp; #it has been included or not?
  	my $extended=0;
  	my @extend;
  	my @ext_dir;
  	my @info=split(/\t/,$all{$reads[$i]}); # the extended raw information
  	for(my $j=0;$j<@reads;$j++){
  		#next if($reads[$j] ne "m54079_171127_175044/57999529/30541_35747");
  		next if($i==$j or defined($tmp{$reads[$j]}));
  		my @pr=split(/\t/, $pattern{$reads[$j]}); #partner sequences
  		my $mm=(scalar(@pr)+1)/3;
  		next if($mm < 2);
  		$nn=(scalar(@cs)+1)/3;
  		my $output=TSD::align->new(\@cs, \@pr, \%dd1, \%dd2, \%ann_len,\%ann_from, \%ann_to, $min_fragment_length,$cutoff_similarity);
  		my $rev=0;
  		if($output->{max_score} < 2 and $output->{max_score} != -1){ # reverse the partner fragments
  			@pr=&revfrag(@pr);
  			$output=TSD::align->new(\@cs, \@pr, \%dd1, \%dd2, \%ann_len,\%ann_from, \%ann_to, $min_fragment_length, $cutoff_similarity);
  			$rev=1;
  		}
  		next if($output->{max_score} < 2);
  		next if($output->{used} == 0);
  		my $a=$output->{a};
  		my $b=$output->{b};
  		my @aa;
  		my @bb;
  		foreach (@$a){ push @aa, $_;}
  		foreach (@$b){ push @bb, $_;}
  		
  		my $na=scalar(@aa);
  		#print join("/",@aa),"\n";
  		#print join("/",@bb),"\n";
  		next if($aa[$na-1] != ($nn-1) and $bb[$na-1]!=($mm-1));
  		next if($aa[0] !=0 and $bb[0]!=0);
  		my $my_min_merge_overlap=$min_merge_overlap;
  		$my_min_merge_overlap=3 if($nn > 5 and $my_min_merge_overlap < 3);
  		if($aa[0]==0 and $bb[0]!=0 and $output->{max_score} >= $my_min_merge_overlap and $merge){
  			$extended=1;
  			my @pinfo=split(/\t/,$all{$reads[$j]}); #partner
  			for(my $zz=3*$bb[0]-1;$zz >=0;$zz--){
  				unshift @cs, $pr[$zz];
  			}
  			if($rev){
  				print join("\t", @pinfo),"\n";
  				@pinfo=&revfrag2(@pinfo);
  				print join("\t", @pinfo),"\n";
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
  		elsif($aa[$na-1]==($nn-1) and $bb[$na-1]!= ($mm-1) and $output->{max_score} >= $my_min_merge_overlap and $merge){
  			$extended=1;
  			my @pinfo=split(/\t/,$all{$reads[$j]}); #partner
  			if($rev){
  				print join("\t", @pinfo),"\n";
  				@pinfo=&revfrag2(@pinfo);
  				print join("\t", @pinfo),"\n";
  			}
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
  		next if($end - $start < 2);
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
  			system("Rscript $lib/TSD/plot_fragment.r $output_dir $output_dir/temp_files/temp_$i.txt $counter $reads[$i]/E $ext $ext2");
  		}
  		else{
  			system("Rscript $lib/TSD/plot_fragment.r $output_dir $output_dir/temp_files/temp_$i.txt $counter $reads[$i] $ext $ext2");
  		}
  		
  		#unlink("temp_files/temp_$i.txt");
  		foreach (keys %tmp){
  			$used{$_}=1;
  		}
  		@cv=@cv[$start..$end];
  		print OUT "seq$counter\t$reads[$i]\t",scalar(@cv), "\t$has\t",mymin(@cv),"\t", coverage(@cv),"\t",join(",",@cv);
  		#my @all=split(/\t/,$all{$reads[$i]});
  		for(my $yy=$start*6; $yy <=($end*6+3); $yy++){
  			next if(($yy-4)%6==0);
  			print OUT "\t",$info[$yy];
  		}
  		print OUT "\n";
  		$counter++;
  		
  	}
  	#print join("\n",keys %used),"\n";
  }
  close OUT;
  bless($self,$class);
	return($self);
}


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
sub revfrag2{ #reverse the pattern
	my @x=@_;
	my $nn=scalar(@x);
	my @re;
	for(my $i=0;$i<$nn;$i=$i+6){
		unshift @re, -1*$x[$i+2];
		unshift @re, -1*$x[$i+3];
		unshift @re, -1*$x[$i+1];
		unshift @re, $x[$i];
		if($i < $nn - 4){
			unshift @re, $x[$i+5];
			unshift @re, $x[$i+4];
		}
	}
	return @re;
}
1;

