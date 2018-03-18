#!/usr/bin/perl -w
use strict;
package align;

sub new{
	my $class=shift;
	my ($cs, $pr, $dd1, $dd2, $site_len, $site_from, $site_to, $min_fragment_length,$cutoff_similarity)=@_;
	my $self={};
	my %can; # defined all the possible match fragments
	my @al;
	my $max_score=-1;
	my $mm=(@$cs+1)/3;
	my $nn=(@$pr+1)/3;
	for(my $p=0;$p<@$pr;$p=$p+3){ # pre determined the all the possible fragment combination
		for(my $q=0;$q < @$cs;$q=$q+3){
			next if($$cs[$q+1] != $$pr[$p+1]);
			if($$cs[$q] eq $$pr[$p] ){ # same name
				$can{"$q\t$p"}=1;
			}
			elsif(defined($$dd1{"$$cs[$q]\t$$pr[$p]"}) or defined($$dd1{"$$pr[$p]\t$$cs[$q]"})){
				if($q==0 or $p==0 or $q == @$cs-2 or $p== @$pr-2 ){
					$can{"$q\t$p"}=1;
				}
				elsif(defined($$dd2{"$$cs[$q]\t$$pr[$p]"}) and $$dd2{"$$cs[$q]\t$$pr[$p]"} > $cutoff_similarity){
					$can{"$q\t$p"}=1;
				}
				elsif(defined($$dd2{"$$pr[$p]\t$$cs[$q]"}) and $$dd2{"$$pr[$p]\t$$cs[$q]"} > $cutoff_similarity){
					$can{"$q\t$p"}=1;
				}
			}
		}
	}
	if(scalar(keys %can) <= 2){
		$self={
			a=>\@al,
			b=>\@al,
			max_score=>$max_score,
			n1=>$mm,
			n2=>$nn,
			used=>0
		};
		bless($self,$class);
		return($self);
	}
	my %breaks1;
	my %breaks2;
	for(my $p=0;$p< @$cs -3;$p=$p+3){
		$breaks1{$p}=$$cs[$p+2];
	}
	for(my $p=0;$p< @$pr -3;$p=$p+3){
		$breaks2{$p}=$$pr[$p+2];
	}
	#dynamic programming
	my @matrix;
	$matrix[0][0]  = 0;
	
	for ( my $j = 1 ; $j <= $nn ; $j++ ) {
        $matrix[0][$j]   = 0;
    }
    for ( my $i = 1 ; $i <= $mm ; $i++ ) {
        $matrix[$i][0]   = 0;
    }
	# begin DP
	my $wh_x=0;
	my $wh_y=0;
	for(my $q=1; $q <= $mm; $q++){
		for(my $p=1; $p <= $nn; $p++){
			my ($u, $d, $m)=(1,1,0);
			my $a=($q-1)*3;
			my $b=($p-1)*3;
			if(defined($can{"$a\t$b"})){ #match
				$m=1;
			}
			if(defined($breaks1{$a-1}) and $breaks1{$a-1}> 200 and abs($breaks1{$a-1} - $$site_len{$$cs[$q-1]}) <200  ){ #not allow gap
				$u=0;
			}
			if(defined($breaks2{$b-1}) and $breaks2{$b-1}> 200 and abs($breaks2{$b-1} - $$site_len{$$pr[$p-1]}) <200 ){ #not allow gap
				$d=0;
			}
			$matrix[$q][$p]=mymax($matrix[$q-1][$p-1] + $m, $matrix[$q][$p-1] - $d, $matrix[$q-1][$p] - $u, 0);
			if($m==0 and $u==1 and $d==1){
				$matrix[$q][$p]=0;
			}
			if($matrix[$q][$p] > $max_score){
				$max_score=$matrix[$q][$p];
				$wh_x=$q;
				$wh_y=$p;
			}
		}
	}
	#print "1: ",join(" ", @$cs),"\n";
	#print "2: ",join(" ", @$pr),"\n";
	#for(my $q=1; $q <= $mm; $q++){
	#	for(my $p=1; $p <= $nn; $p++){
	#		print " $matrix[$q][$p]"
	#	}
	#	print "\n";
	#}
	
	if(($wh_x < $mm and $wh_y < $nn) or $max_score < 2){
		$self={
			a=>\@al,
			b=>\@al,
			max_score=>$max_score,
			n1=>$mm,
			n2=>$nn,
			used=>0
		};
		bless($self,$class);
		return($self);
	}
	
	my $i=$wh_x;
	my $j=$wh_y;
	#print "$wh_x\t$wh_y\n";
	my @xx;
	my @yy;
	push @xx,$i-1;
	push @yy,$j-1;
	while($i > 1 and $j > 1){
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
	
	
	$self={
		a=>\@xx,
		b=>\@yy,
		max_score=>$max_score,
		n1=>$mm,
		n2=>$nn,
		used=>1
	};
	bless($self,$class);
	return($self);
}
sub mymax{
  my @x=@_;
  my $re=-1000;
  foreach (@x){
	$re=$_ if($re < $_);
  }
  return $re;
}
1;
