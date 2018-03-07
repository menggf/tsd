#!/usr/bin/perl -w
use strict;

open DATA,"result.txt" or die;
my %sites;
my %all;
my %len;
while(my $str=<DATA>){
	next if($str!~/insert/);
	chomp $str;
	my @str=split(/\t/,$str);
	next if(@str <=6);
	my $id=shift @str;
	$all{$id}=join("\t", @str);
	my $nn=scalar(@str);
	$len{$id}=0;
	for(my $i=0;$i<$nn; $i=$i+6){
		$sites{"$str[$i]\t$str[$i+2]\t$str[$i+3]"}="";
		$len{$id}+=($str[$i+3]-$str[$i+2]);
		if($i+5 <$nn){
			$len{$id}+=$str[$i+5];
		}
	}
}
close DATA;


system("Rscript deal.r");
open DATA,"new_sites.txt" or die;
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
	for(my $j=$i;$j<@ss;$j++){
		next if($i==$j);
		next if($ann_chr{$ss[$i]} ne $ann_chr{$ss[$j]});
		my @t2=split(/\t/,$ann{$ss[$j]});
		my $from2=$t2[1];
		my $to2=$t2[2];
		next if($to1 < $from2 or $from1 < $from2);
		next if($to2 < $from1 or $from2 < $from1);
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
		next if($r1 < 0.85);
		$dd1{"$ss[$i]\t$ss[$j]"}=$r1;
		$dd2{"$ss[$i]\t$ss[$j]"}=$r2;
		if($r2 > 0.85){
			$freq{$ss[$i]}++;
			$freq{$ss[$j]}++;
		}
	}
}

my %pattern; #reformat into site string
my %has_gap; # used to fill the gap
foreach my $id (keys %all){
	my $str=$all{$id};
	next if($str!~/insert/);
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
			if($str[$i+5] > 200){
				my $ts=$i/6;
				if(defined($has_gap{$id})){
					$has_gap{$id}.="\t$ts";
				}
				else{
					$has_gap{$id}=$ts;
				}
			}
		}
	}
	$pattern{$id}=$s;
}
close DATA;

my %ass;
my %used;
open OUT,">report.txt" or die;
print "Begin alignments\n";
my @reads=sort{$len{$b} <=> $len{$a}}(keys %len);
for(my $i=0;$i<@reads;$i++){
	print "read $i\n";
	#next if($i!=0);
	next if(defined($used{$reads[$i]}));
	my %first; # where is the location
	$first{$reads[$i]}=0;
	my $has=0; #how many supporting reads
	my @tp=split(/\t/,$pattern{$reads[$i]}); # the seed sequences
	my @cs=@tp; # the consensus
	next if(@tp < 9);
	for(my $j=0;$j<@reads;$j++){
		next if($i==$j);
		my @pr=split(/\t/, $pattern{$reads[$j]}); #partner sequences
		my @al=&align(\@cs, \@pr, \%dd1, \%dd2, \%ann_len);
		next if(scalar(@al)==0);
		print "$pattern{$reads[$i]}\n#\n$pattern{$reads[$j]}\n";
		my %can; # defined all the possible match fragments
		for(my $p=0;$p<@pr;$p=$p+3){ # pre determined the all the possible fragment combination
			for(my $q=0;$q < @cs;$q=$q+3){
				if($cs[$q] eq $pr[$p]){
					$can{"$q\t$p"}=1;
				}
				elsif(defined($dd1{"$cs[$q]\t$pr[$p]"}) or defined($dd1{"$pr[$p]\t$cs[$q]"})){
					if($q==0 or $p==0 or $q == @cs-2 or $p== @pr-2 ){
						$can{"$q\t$p"}=1;
					}
					elsif(defined($dd2{"$cs[$q]\t$pr[$p]"}) and $dd2{"$cs[$q]\t$pr[$p]"} > 0.85){
						$can{"$q\t$p"}=1;
					}
					elsif(defined($dd2{"$pr[$p]\t$cs[$q]"}) and $dd2{"$pr[$p]\t$cs[$q]"} > 0.85){
						$can{"$q\t$p"}=1;
					}
				}
			}
		}
		next if(scalar(keys %can)==0);
		my %breaks1;
		my %breaks2;
		for(my $p=0;$p< @cs -3;$p=$p+3){
			$breaks1{$p}=$cs[$p+2];
		}
		for(my $p=0;$p< @pr -3;$p=$p+3){
			$breaks2{$p}=$pr[$p+2];
		}
		my $tag=0;
		
		for(my $q=0;$q < @cs;$q=$q+3){ #start from consensus
			if(defined($can{"$q\t0"})){
				my $nostop=1;
				my $y=0; #a pointer to position in partner
				my $support=1; # number of overlap fragment
				if(scalar(@cs) - $q *3 >= scalar(@pr)){ # consensus sequence is longer
					for(my $x=$q + 3;$x < scalar(@pr) + $q*3; $x=$x+3){
						$y=$y+3;
						if(!defined($can{"$x\t$y"})){
							
							$nostop=0;
							last;
						}
						$support++;
					}
				}
				else{
					for(my $x=$q+3;$x < scalar(@cs);$x=$x+3){
						$y=$y+3;
						if(!defined($can{"$x\t$y"})){
							$nostop=0;
							last;
						}
						$support++;
					}
				}
				if($nostop and $support >= max(2, 0.5 * min((scalar(@pr)+1)/3,(scalar(@cs)+1)/3))){
					$tag=1;
					$used{$reads[$j]}=1;
					$first{$reads[$j]}=$q/3;
					if(scalar(@cs) - $q *3 < scalar(@pr)){
						push @cs, @pr[(scalar(@pr)-$support *3 +1)..(scalar(@pr)-1)];
					}
					if(!defined($ass{$reads[$i]})){
						$ass{$reads[$i]}=$ass{$reads[$j]};
					}
					else{
						$ass{$reads[$i]}.="\t$ass{$reads[$j]}";
					}
					$has++;
					last;
				}
			}
		}
		next if($tag);
	
		for(my $p=3;$p < @pr;$p=$p+3){ #start from new sequence
			if(defined($can{"0\t$p"})){
				my $nostop=1;
				my $y=0;
				my $support=1; # number of overlap fragment
				if(scalar(@pr) -$p *3 >= scalar(@cs)){ # new sequence is longer
					for(my $x=$p+3;$x < scalar(@cs) + $p*3;$x=$x+3){
						$y=$y+3;
						if(!defined($can{"$y\t$x"})){
							$nostop=0;
							last;
						}
						$support++;
					}
				}
				else{ #consensus sequencing is longer
					for(my $x=$p+3;$x < scalar(@pr);$x=$x+3){
						$y=$y+3;
						if(!defined($can{"$y\t$x"})){
							$nostop=0;
							last;
						}
						$support++;
					}
				}
				if($nostop and $support >= max(2,0.5*min( (scalar(@pr)+1)/3,(scalar(@cs)+1)/3) ) ){
					$tag=1;
					$used{$reads[$j]}=1;
					if(!defined($ass{$reads[$i]})){
						$ass{$reads[$i]}=$ass{$reads[$j]};
					}
					else{
						$ass{$reads[$i]}.="\t$ass{$reads[$j]}";
					}
					if(scalar(@pr) -$p *3 >= scalar(@cs)){
						push @cs, @pr[(scalar(@pr)-$support *3 -$p*3 +1)..(scalar(@pr)-1)];
					}
					unshift @cs,$pr[0..($p-1)];
					foreach (keys %first){
						$first{$_}+=$p/3;
					}
					$first{$reads[$j]}=0;
					$has++;
					last;
				}
			}
		}
	}
	if($has >2){
		print "$i\t A:$has\t",scalar(@cs),"\n";
		print OUT "Sites $reads[$i]\n";
		print OUT "Consensus\t0\t",join("\t", @cs),"\n";
		foreach my $rds (@reads){
			next if(!defined($first{$rds}));
			my $tmp="";
			for(my $i=0;$i<$first{$rds};$i++){
				$tmp.="\t\t\t\t\t\t";
			}
			print OUT "$rds\t$first{$rds}\t$tmp\t$all{$rds}\n";
		}
	}	
}
close OUT;

exit(0);

sub align{
	my ($cs, $pr, $dd1, $dd2, $site_len)=@_;
	my %can; # defined all the possible match fragments
	my @al;
	for(my $p=0;$p<@$pr;$p=$p+3){ # pre determined the all the possible fragment combination
		for(my $q=0;$q < @$cs;$q=$q+3){
			if($$cs[$q] eq $$pr[$p]){
				$can{"$q\t$p"}=1;
			}
			elsif(defined($$dd1{"$$cs[$q]\t$$pr[$p]"}) or defined($$dd1{"$$pr[$p]\t$$cs[$q]"})){
				if($q==0 or $p==0 or $q == @$cs-2 or $p== @$pr-2 ){
					$can{"$q\t$p"}=1;
				}
				elsif(defined($$dd2{"$$cs[$q]\t$$pr[$p]"}) and $$dd2{"$$cs[$q]\t$$pr[$p]"} > 0.85){
					$can{"$q\t$p"}=1;
				}
				elsif(defined($$dd2{"$$pr[$p]\t$$cs[$q]"}) and $$dd2{"$$pr[$p]\t$$cs[$q]"} > 0.85){
					$can{"$q\t$p"}=1;
				}
			}
		}
	}
	return(@al) if(scalar(keys %can) <=2);
	my %breaks1;
	my %breaks2;
	for(my $p=0;$p< @$cs -3;$p=$p+3){
		$breaks1{$p}=$$cs[$p+2];
	}
	for(my $p=0;$p< @$pr -3;$p=$p+3){
		$breaks2{$p}=$$pr[$p+2];
	}
	my @matrix;
	$matrix[0][0]  = 0;
	
	my $mm=(@$cs+1)/3;
	my $nn=(@$pr+1)/3;
	for ( my $j = 1 ; $j <= $nn ; $j++ ) {
        $matrix[0][$j]   = 0;
    }
    for ( my $i = 1 ; $i <= $mm ; $i++ ) {
        $matrix[$i][0]   = 0;
    }
	
	my $max_score=-1;
	my $wh_x=0;
	my $wh_y=0;
	for(my $q=1; $q <= $mm; $q++){
		for(my $p=1; $p <= $nn; $p++){
			my ($u, $d, $m)=(0,0,0);
			my $a=($q-1)*3;
			my $b=($p-1)*3;
			if(defined($can{"$a\t$b"})){
				$m=1;
			}
			else{
				$m=-1;
			}
			if(!$breaks1{($q-1)*3-1} and $q != $mm ){
				$u=1;
			}
			if($breaks2{($p-1)*3-1} and $p !=$nn){
				$d=1;
			}
			$matrix[$q][$p]=mymax($matrix[$q-1][$p-1] + $m, $matrix[$q][$p-1] - $d, $matrix[$q-1][$p] - $u, 0);
			if($matrix[$q][$p] > $max_score){
				$max_score=$matrix[$q][$p];
				$wh_x=$q;
				$wh_y=$p;
			}
		}
	}
	return(@al) if($wh_x < $mm and $wh_y < $nn);
	return(@al) if($max_score < 2);
	
	my $i=$wh_x;
	my $j=$wh_y;
	my @xx;
	my @yy;
	push @xx,$i-1;
	push @yy,$j-1;
	while($i >= 1 and $j >= 1){
		last if($matrix[$i][$j]==0);
		if($matrix[$i-1][$j-1] == $matrix[$i][$j] - 1){
			$i=$i-1;
			$j=$j-1;
			unshift @xx, $i-1;
			unshift @yy, $j-1;
		}
		elsif($matrix[$i-1][$j] == $matrix[$i][$j] - 1){
			$i=$i-1;
			unshift @xx, $i-1;
			unshift @yy, -1;
		}
		elsif($matrix[$i][$j-1] == $matrix[$i][$j] - 1){
			$j=$j-1;
			unshift @xx, -1;
			unshift @yy, $j-1;
		}
		else{
			last;
		}
	}
	return(@al) if($xx[0]!=0 and $yy[0] !=0);
	print "A: ",join("\t", @xx),"\n";
	print "B: ",join("\t", @yy),"\n";
	return((@xx,@yy));
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
sub mean{
  my @x=@_;
  my $sum=0;
  foreach (@x){
	$sum=$_;
  }
  return($sum/scalar(@x));
}
