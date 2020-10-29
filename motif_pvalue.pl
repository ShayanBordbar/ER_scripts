#!/usr/bin/perl

# implementation based on the recursive p-value calculation in the paper MONKEY

use strict;
use warnings;
use lib "/home/sinhas/sinhas/lib/";
use Bio::Seq;
use Bio::SeqIO;

if ($#ARGV+1 != 2) {
    print STDERR "Usage: motif_pvalue.pl wtmx_file out_dir\n";
    exit 1;
}

#my $BINCNT = 10000.0;
my $BINCNT = 1000.0;
my @chars = ("A","C","G","T");
my $PCOUNT = 0.05;	# default pseudo count
my $pwms_f = $ARGV[0];
my $out_d = $ARGV[1];

# make output directory
system("mkdir -p $out_d");

# read background nucleotide frequencies
my %hsbkg = ("A"=>0.3,"C"=>0.2,"G"=>0.2,"T"=>0.3);  
#my %hsbkg = ("A"=>0.30755,"C"=>0.19195,"G"=>0.19188,"T"=>0.30862);  
#my %hsbkg = ();
#$hsbkg{'A'} = 0.25;
#$hsbkg{'C'} = 0.25;
#$hsbkg{'G'} = 0.25;
#$hsbkg{'T'} = 0.25;

# Debug
printf STDERR "Backgroupd(A,C,G,T): %f %f %f %f\n",$hsbkg{"A"},$hsbkg{"C"},$hsbkg{"G"},$hsbkg{"T"};

my $mname = "";
my $mpcount = $PCOUNT;    
my @arpwm = (); # Array of hash table for each pwn position
my %hs_pwms = ();
open(FILE, "$pwms_f");
while(<FILE>) {
    chomp;
    if (substr($_,0,1) eq ">") {
        my @ar = split(/\s+/, substr($_,1));
		$mname = $ar[0];
        if ($#ar == 2) { $mpcount = $ar[2]; } 
        print STDERR "Motif:$ar[0] size:$ar[1] pseudocount:$mpcount\n";
    } elsif (substr($_,0,1) eq "<") {
        # copy array and store it into hashtable
		my @artmp = @arpwm;
		$hs_pwms{$mname} = \@artmp;

        @arpwm = ();
        $mpcount = $PCOUNT;
        $mname = "";
    } else {
        my @ar = split(/\s+/);
        my ($nA2,$nC2,$nG2,$nT2) = ($ar[0]+$mpcount,$ar[1]+$mpcount,$ar[2]+$mpcount,$ar[3]+$mpcount);
        my $sum = $nA2+$nC2+$nG2+$nT2;
        my %hstmp = ("A"=>$nA2/$sum,"C"=>$nC2/$sum,"G"=>$nG2/$sum,"T"=>$nT2/$sum);
        push(@arpwm, \%hstmp); 
    }
}
close(FILE);

foreach my $mname (keys %hs_pwms) {
	my $artmp = $hs_pwms{$mname};
    @arpwm = @$artmp;
    print STDERR "[Motif:$mname]\n";

    # calculate p-value
    # 1. Make matrix M_ic (log likelihood score for each nucleotide(c) at each position(i)
    my %M_ic = ();
    my $min = 0; 
    my $ccnt = 0;
    foreach my $nt (@chars) {
        for(my $i=1; $i<=@arpwm; $i++) {
            # convert llr to int value
            $M_ic{$nt}{$i} = int($BINCNT*(log($arpwm[$i-1]{$nt}) - log($hsbkg{$nt})));
		    if ($ccnt == 0 || $M_ic{$nt}{$i} < $min) { $min = $M_ic{$nt}{$i}; }
		    $ccnt++;
	    }
    }
    print STDERR "Minimum value in M_ic matrix:$min\n";
	my $totalmin = $min*@arpwm;

    # 2. Transform all values in M_ic to non-negative values
    my $max_val = 0;
    #foreach my $nt (@chars) {
    #	for(my $i=0; $i<@arpwm; $i++) {
    for(my $i=1; $i<=@arpwm; $i++) {
	    my $max = 0;
	    my $cnt = 0;
	    foreach my $nt (@chars) {
		    $M_ic{$nt}{$i} -= $min;
		    if ($cnt == 0 || $M_ic{$nt}{$i} > $max) { $max = $M_ic{$nt}{$i}; }
		    $cnt++;
	    }

	    $max_val += $max;
    }
    print STDERR "Maximum value in M_ic matrix(after transformation):$max_val\n";

    # 3. Make pdf matrix for all scores from 0 to max_val for each position
    my %pdf = (); # pdf values 
    # initialization
    for(my $i=0; $i<=$max_val; $i++) {
	    for(my $j=1; $j<=@arpwm; $j++) { $pdf{$i}{$j} = -1.0; }
   }

    # write result to file
    my $out_f = $out_d."/".$mname.".txt";
	open(FILE, ">$out_f");

    # Background nuc frequency
	print FILE ">Bkg\n";
    foreach my $nt (@chars) {
        print FILE "$nt\t$hsbkg{$nt}\n";
	}
	print FILE "<\n";

    # M_ic matrix
	print FILE ">Mic\n";
    for(my $i=1; $i<=@arpwm; $i++) {
	    foreach my $nt (@chars) {
		    my $val = $M_ic{$nt}{$i};
			print FILE "$val\t";
	    }
		print FILE "\n";
    }


    print FILE "<\n";

    print FILE ">Pvalue\n";
   my %pvalue = (); # cdf values
    my $msize = scalar(@arpwm);
    for(my $i=$max_val; $i>=0; $i--) {
	    if ($i == $max_val) {
            $pvalue{$i} = prob_iS($i, \%M_ic, \%hsbkg, $msize, \@chars, \%pdf);
	    } else {
		    my $newp = prob_iS($i, \%M_ic, \%hsbkg, $msize, \@chars, \%pdf);
		    $pvalue{$i} = $newp + $pvalue{$i+1};
#print STDERR "$i $newp $pvalue{$i+1}\n";			
	    }

	    printf FILE "%d\t%.10f\n",$i+$totalmin,$pvalue{$i};
    }
    print FILE "<\n";

	close(FILE);
}

#############################################################################
sub prob_iS {
    my $scr = shift;
    my $Mic = shift;
	my $hsb = shift;
    my $pos = shift;
    my $chars = shift;
	my $pdf = shift;

#print STDERR "[scr:$scr]\n";

    if ($scr < 0) { return 0.0; }

    # base condition
    if ($pos == 0) {
        if ($scr == 0) { return 1; }
		else { return 0; }
	}

    # reuse computation
	#if (!$$pdf{$scr}{$pos}) { print STDERR "==>$scr $pos\n"; }
	if ($$pdf{$scr}{$pos} != -1) { return $$pdf{$scr}{$pos}; } 

    # recursive call
	$$pdf{$scr}{$pos} = 0.0;
    foreach my $nt (@$chars) {
	    my $scr_p = $scr - $$Mic{$nt}{$pos};
	    my $pr1 = prob_iS($scr_p,$Mic,$hsb,$pos-1,$chars,$pdf);	
	    my $pr2 = $$hsb{$nt};
		my $tmp = $pr1*$pr2;
#print STDERR "$scr_p $pos $nt $pr1 $pr2 $tmp\n";	
	    $$pdf{$scr}{$pos} += $pr1 * $pr2;
	}

    #print STDERR "$scr $pos $$pdf{$scr}{$pos}\n";
    return $$pdf{$scr}{$pos};
}
