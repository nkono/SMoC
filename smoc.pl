use strict;
use Getopt::Long;

my %opts = (word => 100, cover => 5, error => 1, limit => 8000, output => "SMoC");
GetOptions(\%opts, qw( seed|s=s read|r=s output|o=s word|w=i cover|c=i error|e=i limit|l=i help|h) ) or exit 1;

if($opts{help}) {
    die &show_help();
}

my $dieflag;
foreach my $field ( qw( seed read ) ){
    $dieflag .= "Error: [$field] file is required.\n" if ! exists $opts{$field};
}
die $dieflag if $dieflag;

my $f1 = $opts{read};
my $seed = $opts{seed};
my $w = $opts{word};
my $cov = $opts{cover};
my $error = $opts{error};
my $limit = $opts{limit};
my $name = $opts{output};
my $e = 10;

print "Seed (-seed):\t".$seed."\n";
print "Read (-read):\t".$f1."\n";
print "K-mer (-word):\t".$w."\n";
print "error (-error):\t".$error."\n";
print "Limit (-limit):\t".$limit."\n";

my $output = $name."_w".$w."_c".$cov."_contig.fasta";
my $output_res = $name."_w".$w."_c".$cov."_contig.res";

open(OUT, "> $output");
open(RES, "> $output_res");

my @reads;

my $flag = 1;
my @reads;
open(FILE, $f1);
while(<FILE>){
    if(/^\@/){
	$flag = 1;
    }elsif($flag == 1){
	chomp();
	push(@reads, $_);
	$flag = 0;
    }
}
close(FILE);

my $seq;
open(FILE, "$seed");
while(<FILE>){
    chomp;
    next if />/;
    $seq .= lc($_);
}
close FILE;

my %store;
my %pattern;
my $round = 0;
&extender($seq, "root");

sub extender{
    my $seq = shift;
    my $id = shift;
    my $tmp = shift;
    my $compare = uc(substr($seq, -$w, $w));
    my $revcomp = uc(complement(substr($seq, -$w, $w)));
    my $after = {};
    for my $read (@reads){
	my $readseq = '';
	if ($read =~ /$compare/){
	    $readseq = $read;
	}elsif($read =~ /$revcomp/){
	    $readseq = complement($read);
	}

	if(length($readseq)){
	    my $mismatch = 0;
	    my $j = 0;
	    my $index = index($readseq, $compare);
	    for (my $i = $index - 1; $i >= 0 ;$i --){
		$j ++;
		next if !substr($seq, -$w - $j, 1);
		$mismatch ++ if (lc(substr($seq, -$w - $j, 1)) ne lc(substr($readseq, $i, 1)));
	    }
	    unless($mismatch >= $error){
		for (0..$e){
		    last if (length($readseq) < $index + $w + $_);
		    $after->{$_}->{substr($readseq, $index + $w + $_, 1)} ++;
		}
	    }
	}
    }

    &dp($after);

    my $stop;
    $pattern{$id}{SEQ} = $seq;
    for (0..$e){
	my @sort = sort {$after->{$_}->{$b} <=> $after->{$_}->{$a}} keys %{$after->{$_}};
	if(length($sort[0]) < 1){
	    $stop = "There were no appropriate bases.";
	    last;
	}
	my $sum = ($after->{$_}->{"A"} + $after->{$_}->{"T"} + $after->{$_}->{"G"} + $after->{$_}->{"C"});
	my $per = ($after->{$_}->{$sort[0]} / $sum) * 100;
	my $per2nd = (($after->{$_}->{$sort[0]} + $after->{$_}->{$sort[1]}) / $sum) * 100;
	print RES "cov:".$sum."\tper:".$per." ".$sort[0];
	if($sort[1]){
	    print RES "\t2nd per:".$per2nd." ".$sort[0]." ".$sort[1]."\n";
	}else{
	    print RES "\n";
	}
	
	# stop for covrage
	if($sum < $cov){
	    $stop = "Coverage was low.";
	    last;
	}

	if($per >= 70){ # top hit
	    $pattern{$id}{SEQ} .= $sort[0];
	    $pattern{$id}{EXTEND} .= $sort[0];
	}elsif($per2nd >= 90){ # pattern
	    $pattern{$id."-1"}{SEQ} = $pattern{$id}{SEQ}.$sort[0];
	    $pattern{$id."-2"}{SEQ} = $pattern{$id}{SEQ}.$sort[1];
	    $stop = "pattern";
	    last;
	}else{
	    $stop = "There were no appropriate bases.";
	    last;
	}
    }

    my $cons = uc(substr($pattern{"root"}{SEQ}, -1 * $w, $w));
    print RES "CONS\t".$cons."\n";
    my $newseq = $pattern{$id}{SEQ};
    print RES $tmp."\t".$stop."\n";
    print RES seqsplit($id, $newseq);

    my $rootnum = $id =~ s/-/-/g;
    my ($p1_ext, $p2_ext);

    if(length($newseq) > 5000){
	print RES "//\n";
	print RES "Because: long\n";
	if($tmp =~ /tmp/){
	    my $lastN = $1 if $tmp =~ /tmp(\w)/;
	    my $returnseq = $lastN.$pattern{$id}{EXTEND};
	    print RES "Return: ".$returnseq."\nReturn0: ".$lastN."\nReturn1: ".$pattern{$id}{EXTEND}."\n";
	    return $returnseq;
	}else{
	    print OUT seqsplit($id."|LONG", $newseq);
	}
    }elsif($limit > 1 && -s $output > $limit){ # size $limit
	print RES "//\n";
	print RES "Because: file size\n";
	if($tmp =~ /tmp/){
	    my $lastN = $1 if $tmp =~ /tmp(\w)/;
	    my $returnseq = $lastN.$pattern{$id}{EXTEND};
	    print RES "Return: ".$returnseq."\nReturn0: ".$lastN."\nReturn1: ".$pattern{$id}{EXTEND}."\n";
	    return $returnseq;
	}else{
	    print OUT seqsplit($id."|SIZE", $newseq);
	}
    }elsif($rootnum > 10){ # variety
	print RES "//\n";
	print RES "Because: variety\n";
	if($tmp =~ /tmp/){
	    my $lastN = $1 if $tmp =~ /tmp(\w)/;
	    my $returnseq = $lastN.$pattern{$id}{EXTEND};
	    print RES "Return: ".$returnseq."\nReturn0: ".$lastN."\nReturn1: ".$pattern{$id}{EXTEND}."\n";
	    return $returnseq;
	}else{
	    print OUT seqsplit($id."|VARIETY", $newseq);
	}
    }elsif($stop ne "pattern" && length($newseq) != length($seq)){
	print RES "Extended\n";
	&extender($newseq, $id, $tmp);
    }elsif($stop){
	print RES "//\n";
	if($tmp =~ /tmp/){
	    my $lastN = $1 if $tmp =~ /tmp(\w)/;
	    my $returnseq = $lastN.$pattern{$id}{EXTEND};
	    print RES "Return: ".$returnseq."\nReturn0: ".$lastN."\nReturn1: ".$pattern{$id}{EXTEND}."\n";
	    return $returnseq;
	}elsif($stop eq "pattern"){
	    my ($id1, $id2) = ($id."-1", $id."-2");
	    my $endflag;

	    my $addN1 = substr($pattern{$id1}{SEQ}, -1, 1);
	    $p1_ext = &extender($pattern{$id1}{SEQ}, $id1, "tmp".$addN1);
	    print RES "Return tmp".$id1."\n";
	    print RES "=== CONS CHECK ===\n";
	    print RES "newseq: ".$newseq."\n";
	    print RES "p1_ext: ".$p1_ext."\n";
	    print RES "FULL: ".$newseq.$p1_ext."\n";
	    print RES "CONS?: ".substr($newseq.$p1_ext, -1 * $w, $w)."\n";
	    if(uc(substr($newseq.$p1_ext, -1 * $w, $w)) eq $cons){ # consensus stop
		$endflag += 1;
		print RES "Because: end\n";
		print RES "p1_ext\t".$p1_ext."\n";
		print OUT seqsplit($id."|CONS", $newseq.$p1_ext);
	    }elsif(length($newseq.$p1_ext) > 60000){ # long stop
		$endflag += 1;
		print RES "Because: too long\n";
		print RES "p1_ext\t".$p1_ext."\n";
		print OUT seqsplit($id."|LONG", $newseq.$p1_ext);
	    }
	    $pattern{$id1}{EXTEND} = $addN1.$pattern{$id1}{EXTEND};

	    my $addN2 = substr($pattern{$id2}{SEQ}, -1, 1);
	    $p2_ext = &extender($pattern{$id2}{SEQ}, $id2, "tmp".$addN2);
	    print RES "Return tmp".$id2."\n";
	    if(uc(substr($newseq.$p2_ext, -1 * $w, $w)) eq $cons){ # consensus stop
		$endflag += 2;
		print RES "Because: end\n";
		print RES "p2_ext\t".$p2_ext."\n";
		print OUT seqsplit($id."|CONS", $newseq.$p2_ext);
	    }elsif(length($newseq.$p2_ext) > 60000){
		$endflag += 2;
		print RES "Because: too long\n";
		print RES "p2_ext\t".$p2_ext."\n";
		print OUT seqsplit($id."|LONG", $newseq.$p2_ext);
	    }
	    $pattern{$id2}{EXTEND} = $addN2.$pattern{$id2}{EXTEND};	    

	    $tmp = "";
	    
	    print RES "========== LOOP CHECK ==========\n";
	    print RES $pattern{$id}{EXTEND}."\n";
	    print RES $p1_ext."\n";
	    print RES $p2_ext."\n";
	    print RES "========== LOOP DUMP  ==========\n";
	    for my $tmp1 (keys %store){
		for my $tmp2 (keys %{$store{$tmp1}}){
		    for my $tmp3 (keys %{$store{$tmp1}{$tmp2}}){
			print RES $tmp1."\n";
			print RES $tmp2."\n";
			print RES $tmp3."\n";
			print RES "\n";
		    }
		}
	    }
	    print RES "\n";
	    

	    if($store{$pattern{$id}{EXTEND}}{$p1_ext}{$p2_ext} || $store{$pattern{$id}{EXTEND}}{$p2_ext}{$p1_ext}){ # loop stop
		print RES "Because: loop\n";
		print OUT seqsplit($id."|LOOP", $newseq);
	    }else{
		$store{$pattern{$id}{EXTEND}}{$p1_ext}{$p2_ext} ++ if length($pattern{$id}{EXTEND}) > $w;
		&extender($newseq.$p1_ext, $id1) if $endflag == 0 || $endflag == 2;
		&extender($newseq.$p2_ext, $id2) if $endflag == 0 || $endflag == 1;
	    }
	}else{
	    print RES "Because: ".$stop."\n";
	    print OUT seqsplit($id."|".$stop, $newseq);
	}
    }else{
	print RES "Here\n";
	&extender($newseq, $id, $tmp);
    }
}


sub complement{
    my $nucl = reverse(shift);
    $nucl =~ tr/[atgcATGC]/[tacgTACG]/;
    return $nucl;
}

sub seqsplit{
    my $name = shift;
    my $nucl = shift;

    my $seqsplit .= ">".$name."\n";
    for(my $i = 0; $i < length($nucl); $i = $i + 60){
        $seqsplit .= substr($nucl, $i, 60)."\n";
    }

    return $seqsplit;
}

sub dp{
    my $dptmp = shift;
    print RES "\n";
    foreach my $p1 (sort {$a <=> $b} keys %{$dptmp}){
        foreach my $p2 (keys %{$dptmp->{$p1}}){
            print RES $p1." ".$p2."\t".$dptmp->{$p1}->{$p2}."\n";
        }
        print RES "\n";
    }
    print RES "\n";
}




###################################
## Help
sub show_help {
    my $help_doc = <<EOF;
Program: SMoC (Spidroin Motif Collection)
Version: 1.0
Usage:   perl SMoC --seed <fasta> --read <fastq> <command> [options]
Command: -s, --seed FILE   Seed file (format: FASTA), required
         -r, --read FILE   Read file (format: FASTQ), required
         -o, --output STR  Output name (default: SMoC)
	 -w, --word INT    K-mer size (default: 100)
	 -c, --cover INT   Minimum coevrage (default: 5)
	 -e, --error INT   Maximum mismatch (default: 1)
         -l, --limit INT   Limit size (default: 8000)

License: GNU General Public License
         Copyright (C) 2021
         Institute for Advanced Biosciences, Keio University, JAPAN
Author:  Nobuaki Kono, Kazuharu Arakawa
EOF
return $help_doc;
}

