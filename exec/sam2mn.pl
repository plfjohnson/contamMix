#!/usr/bin/perl -w

# Copyright 2012,2013 Philip Johnson.
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published
#     by the Free Software Foundation, either version 3 of the License,
#     or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License at <http://www.gnu.org/licenses/> for
#     more details.

# Based off earlier (more general) sam2mtdiff.pl + implement baseQ thresholding
#
# Inputs: take BAM/SAM alignment against consensus, fasta-formatted
# multiple alignment (mafft?) of full-length mtDNA genomes (including consensus), id for consensus in multiple alignment (by default uses reference id in BAM).
#
# Output: one line per read with the following columns:
#  <read id>
#  <read len>
#  <# matches to consensus>
#  <# non-matches to consensus>
#  <# matches to genome 1>
#  <# non-matches to genome 1>
#  ...
#  <# matches to genome n>
#  <# non-matches to genome n>
#
# And the last line contains an estimate of error

use strict;
use Getopt::Long;

use constant { #name SAM fields
    ID => 0,
    FLAG => 1,
    REFID => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    MRNM => 6,
    MPOS => 7,
    ISIZE => 8,
    SEQ => 9,
    QUAL => 10,
    TAG => 11
};

my ($bamFN, $contamFN, $consId, $baseQthreshold);
my $mapQthreshold = 30;
my $transverOnly;
my $trimBases = 0; # number of bases to trim off each end of read (simple hack to minimize aDNA damage)
my $warnId = 1; #warn if reference id does not match supplied consId
if (!GetOptions('bam=s' => \$bamFN,
                'ref=s' => \$consId,
                'fa=s' => \$contamFN,
                'mapq=i' => \$mapQthreshold,
                'baseq=i'=> \$baseQthreshold,
                'trimBases=i' => \$trimBases,
                'transverOnly'=> \$transverOnly) ||
    !defined($bamFN  &&  $contamFN)) {
    print "usage: $0 -bam <adna.bam> -fa <genomes.fa> [-ref <consensus id>] [...]\n".
        "where ... could be:\n".
        "\t-mapq <30>  minimum map quality (defaults to 30)\n".
        "\t-baseq <>  minimum base quality (defaults to no filter)\n".
        "\t-transverOnly  (only count transversions)\n".
        "\t-trimBases <0>  num bases to trim/ignore at ends of each read (think aDNA damage)\n".
        "\n";
    exit 1;
}

## load all potential contaminate mt genomes into memory (311 * 16e3 = ~5MB -- no prob)
my (@genomes, $alignedLength, $alignedConsensus);
{
    my ($acc, $seq);
    open INFA, $contamFN or die "$! attempting to open '$contamFN'\n";
    while (<INFA>) {
        if (/^>/) {
            push @genomes, [$acc, $seq] if defined $seq;
            ($acc) = />(.+?)\s/;
            if ($acc =~ /|/) { # try to pick out actual accession
                foreach my $id (split /\|/, $acc) {
                    $acc = $id if length($id)>2  &&  $id =~ /^\w/;
                }
            }
            $seq = [];
        } else {
            chomp;
            $_ = uc $_;
            tr/ACGTN-/N/c; #restrict to ACGTN-
            push @$seq, split //, $_;
        }
    }
    push @genomes, [$acc, $seq] if defined $seq;
    close INFA;
    
    $alignedLength = scalar(@{$genomes[0][1]});
}

# extra check that genomes are properly multiply-aligned
for (my $j = 1;  $j < @genomes;  ++$j) {
    if ($alignedLength != scalar(@{$genomes[$j][1]})) {
        die "authentic + potential contaminant genomes do not appear to be aligned -- they are different lengths!\n"
    }
}


# propogate a N in any genome to all ungapped genomes & note status of
# following at positions across all genomes:
my @fixedPositions = (0) x $alignedLength;
my @transverPositions = (0) x $alignedLength;
my @tiPositions = (0) x $alignedLength;
my @indelPositions = (0) x $alignedLength;
for (my $i = 0;  $i < $alignedLength;  ++$i) {
    my $foundN = 0;
    my %bases = ($genomes[0]->[1][$i] => 1);
    for (my $j = 0;  $j < @genomes;  ++$j) {
        if ($genomes[$j]->[1][$i] eq 'N') {
            $foundN = 1;
            last;
        } elsif ($genomes[$j]->[1][$i] ne $genomes[0]->[1][$i]) {
            $bases{$genomes[$j]->[1][$i]} = 1;
        }
    }
    if ($foundN) {
        for (my $j = 0;  $j < @genomes;  ++$j) {
            if ($genomes[$j]->[1][$i] ne '-') {
                $genomes[$j]->[1][$i] = 'N';
            } else {
                $indelPositions[$i] = 1;
            }
        }
    } else {
        if (exists($bases{'-'})) {
            $indelPositions[$i] = 1;
        } elsif (keys(%bases) == 1) {
            $fixedPositions[$i] = 1;
        } elsif (keys(%bases) == 2  &&
                 !((exists($bases{'C'})  &&  exists($bases{'T'}))  ||
                   (exists($bases{'G'})  &&  exists($bases{'A'})))) {
            $transverPositions[$i] = 1;
        } elsif (keys(%bases) == 2) {
            $tiPositions[$i] = 1;
        }
    }
}

# print "#fixed : ", scalar(grep {$_} @fixedPositions), "\n";
# print "#transv: ", scalar(grep {$_} @transverPositions), "\n";
# print "#indel : ", scalar(grep {$_} @indelPositions), "\n";
# print "#ti    : ", scalar(grep {$_} @tiPositions), "\n";
# print "#other : ", @fixedPositions-(scalar(grep {$_} @indelPositions)+scalar(grep {$_} @transverPositions)+scalar(grep {$_} @fixedPositions)+scalar(grep {$_} @tiPositions)), "\n";


my @alnMap; # for mapping between alignment & consensus coordinates


## iterate through all reads and compare to each potential contaminate genome
my $inferringConsensus = 0; #whether we are inferring the consensus id
my ($errDenom, $errNumer) = (0,0); #for simple estimate of errors
open INBAM, ($bamFN =~ /.sam$/) ? $bamFN : "samtools view $bamFN |" or
    die "$! attempting to open '$bamFN'\n";
while (<INBAM>) {
    my @F = split /\t/;
    next if ($F[MAPQ] < $mapQthreshold  ||  $F[FLAG] & 0x4); #"unmapped" flag
    --$F[POS]; #convert to 0-based coords;

    # on first line, identify consensus and initialize mapping
    if (@alnMap == 0) {
        if ($warnId  &&  defined($consId)  &&  $F[REFID] ne $consId) {
            print STDERR "Warning: BAM reference id ('$F[REFID]') does not match supplied consensus id ('$consId')\n";
        } else {
            $consId = $F[REFID];
            $inferringConsensus = 1;
        }
        
        my $consI;
        for ($consI = 0;  $consI < @genomes  &&  $genomes[$consI][0] ne $consId;
             ++$consI) {};
        if ($consI == @genomes) {
            if ($inferringConsensus) {
                die "Error: SAM/BAM reference id ('$consId') does not exist in the multiple alignment.  Override with 'consId' command line option if necessary -- but ensure both identifiers do truly correspond to the consensus sequence!\n";
            } else {
                die "Error: Couldn't find user-supplied consensus '$consId' in multiple alignment!\n";
            }
        }
        $alignedConsensus = $genomes[$consI][1];
        # ensure consensus at first position in array
        unshift(@genomes, $genomes[$consI]);
        splice(@genomes, $consI+1, 1);

        # make mapping between consensus and alignment coordinates
        my $seqI = 0;
        for (my $alnI = 0;  $alnI < $alignedLength;  ++$alnI) {
            if ($alignedConsensus->[$alnI] ne '-') {
                $alnMap[$seqI++] = $alnI;
            }
        }
    }


    # Expand read into @genomes alignment space using CIGAR string
    # (which gives read::consensus mapping) and @alnMap (which gives
    # consensus::alignment mapping)
    #print "cigar: ", $F[CIGAR], "\n";
    my @bases;
    my $discard = 0;
    my ($rPos, $sPos) = (0, $F[POS]);
    while ($F[CIGAR] =~ /(\d+)([SMDI])/g) {
        if ($2 eq 'S') { #soft clipping
            $rPos += $1;
        } elsif ($2 eq 'M') { #match
            die "error processing read '$F[ID]' -- mapped position (".($sPos+1).") is off end of reference sequence (length ".scalar(@alnMap).")!\n" if ($sPos >= @alnMap);
            my $aPos = $alnMap[$sPos];
            my @match = split //, substr $F[SEQ], $rPos, $1;
            if (defined($baseQthreshold)) { #"N"-out bases below threshold
                if ($F[QUAL] eq '*') {
                    print STDERR "Warning: cannot use base Q threshold since sam/bam does not include base qualities for read '$F[ID]'\n";
                } else {
                    my @quals = split //, substr $F[QUAL], $rPos, $1;
                    for (my $k = 0;  $k < @quals;  ++$k) {
                        if (ord($quals[$k])-33 < $baseQthreshold) {
                            $match[$k] = 'N';
                        }
                    }
                }
            }
            
            if ($aPos + $1 > @$alignedConsensus) {
                warn "warning: read '$F[ID]' mapping is partially off end of reference sequence (discarding read)!\n" ;
                $discard = 1;
                last;
            } else {
                for (my $i = 0;  $i < $1;  ++$i) {
                    while ($alignedConsensus->[$aPos++] eq '-') {
                        push @bases, '-';
                    }
                    push @bases, $match[$i];
                }
                $rPos += $1;
                $sPos += $1;
            }
        } else { #indel -- discard due to potential flaky alignment
            $discard = 1;
            last;
        }
    }

    #print $alnMap[$F[POS]], "\t", $alnMap[$F[POS]]+@bases-1, "\n";
    #print @bases, "\n";
    #print @{$alignedConsensus}[($alnMap[$F[POS]])..($alnMap[$F[POS]]+@bases-1)], "\n";
    #print "----------\n";
    #if requested, trim bases from each end of read
    my $alnStart = $alnMap[$F[POS]];
    my $aLeft = 0;
    for (my $sI=0;  $aLeft < @bases  &&  $sI < $trimBases;  ++$aLeft) {
        ++$sI if $bases[$aLeft] ne '-';
    }
    my $aRight = @bases-1;
    for (my $sI=0;  $aRight >= 0  &&  $sI < $trimBases;  --$aRight) {
        ++$sI if $bases[$aRight] ne '-';
    }
    if ($aLeft == @bases  ||  $aRight < 0) {
        $discard = 1;
    } else {
        $alnStart += $aLeft;
        @bases = splice @bases, $aLeft, $aRight-$aLeft+1;
    }

    # finally summarize read!
    my @readSummary;
    if (!$discard) {
        for (my $j = 0;  $j < @genomes;  ++$j) {
            my ($M,$N,$I) = (0,0,0);
            for (my ($k, $alnPos) = (0,$alnStart);  $k < @bases;
                 ++$k, ++$alnPos) {
                if ($indelPositions[$alnPos]) {
                    $discard = 1;
                    last;
                }
                my $genomeB = $genomes[$j][1][$alnPos];
                next if (defined $transverOnly  &&
                         !$transverPositions[$alnPos]);
                if ($bases[$k] eq $genomeB) {
                    ++$M if ($genomeB ne 'N') #only count match if not N-to-N
                } else {
                    ++$N if ($bases[$k] ne 'N'  &&  $genomeB ne 'N');
                }
            }
            push @readSummary, [$M, $N];
        }
    }
    if ($discard) {
        #print STDERR "fyi: discarding read $F[ID]\n";
        next;
    }

    # simple estimate of error rate
    my $allFixed = 1;
    for (my $k = 0;  $k < @bases;  ++$k) {
        my $alnPos = $k + $alnStart;
        if ($fixedPositions[$alnPos]) {
            my $genomeB = $genomes[0][1][$alnPos];
            if ($bases[$k] ne 'N'  &&  $genomeB ne 'N') {
                ++$errDenom;
                ++$errNumer if ($bases[$k] ne $genomeB)
            }
        } else {
            $allFixed = 0;
        }
    }
    
    # report read only if informative
    if (!$allFixed  &&  $readSummary[0][0] + $readSummary[0][1] > 0) {
        print $F[ID], "\t";
        for (my $j = 0;  $j < @genomes;  ++$j) {
            print join(' ', @{$readSummary[$j]}), ' ';
        }
        print "\n"
    } elsif (($readSummary[0][0] + $readSummary[0][1]) > 30  &&
        $readSummary[0][1]/($readSummary[0][0] + $readSummary[0][1]) > 0.50) {
        warn("implausibly high number of mismatches in read '$F[ID]'\n");
    }
}
close INBAM;


print "#error:\t";
if ($errDenom > 0) {
    print sprintf('%.4f', $errNumer/$errDenom), "\t$errDenom";
} else { 
    print "0\t0";
}
print "\n";
