#!/usr/bin/perl

# John M. Gaspar (jsh58@unh.edu)
# June 2013

# This program takes a list of mapping files
#   intended for use with QIIME and converts
#   them for use with FlowClus.

# Primers are named in order (0, 1, 2, ...)
#   -- there are no names in the QIIME file!

use strict;
use warnings;

die "Error -- please list file names on command line\n"
  if (!@ARGV); 

# make output file
my $out = "master.csv";
my @spl = split(/\./, $out);
my $z = 1;
while (-f $out) {
  $out = $spl[0] . $z . "." . $spl[1];
  $z++;
}

# variables for saving data
my $num = 0;
my @primer = ();
my @reverse = ();
my @mids = ();
my @samples = ();

# analyze files
while (my $file = shift @ARGV) {

  if (! open(IN, $file)) {
    print "$file: Warning -- cannot open: $!\n";
    next;
  }

  # initialize order
  my @order;
  for (my $x = 0; $x < 4; $x++) {
    $order[$x] = -1;
  }

  my $line = <IN>;
  chomp $line;
  my @spl = split("\t", $line);

  # load order
  for (my $x = 0; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "#SampleID") {
      $order[0] = $x;
    } elsif ($spl[$x] eq "BarcodeSequence") {
      $order[1] = $x; 
    } elsif ($spl[$x] eq "LinkerPrimerSequence") {
      $order[2] = $x;
    } elsif ($spl[$x] eq "ReversePrimer") {
      $order[3] = $x;
    }
  }

  # check for missing information
  for (my $x = 0; $x < 4; $x++) {
    if ($order[$x] == -1) {
      if ($x == 3) {
        print "$file: Warning -- no reverse primer specified\n",
          "\tDo not specify search for reverse primer\n";
      } else {
        die "$file: Error -- information missing\n";
      }
    }
  }

  # analyze lines
  while ($line = <IN>) {
    chomp $line;
    next if (substr($line, 0, 1) eq "#");
    my @spl = split("\t", $line);
    my $prim = $spl[$order[2]];

    # check existing primers
    my $x;
    for ($x = 0; $x < $num; $x++) {
      if ($primer[$x] eq $prim) {
        push @{$mids[$x]}, $spl[$order[1]];
        push @{$samples[$x]}, $spl[$order[0]];
        last;
      }
    }

    if ($x == $num) {
      $primer[$num] = $prim;
      if ($order[3] != -1) {
        $reverse[$num] = $spl[$order[3]];
      } else {
        $reverse[$num] = 0;
      }
      push @{$mids[$num]}, $spl[$order[1]];
      push @{$samples[$num]}, $spl[$order[0]];
      $num++;
    }

  }

}

close IN;

# print output
die "Error -- no files processed\n" if (! scalar @primer);
open(OUT, ">$out");
for (my $x = 0; $x < scalar @primer; $x++) {
  print OUT "primer,$x,$primer[$x]\n";
  print OUT "reverse,$reverse[$x]\n" if ($reverse[$x]);
  for (my $y = 0; $y < scalar @{$mids[$x]}; $y++) {
    print OUT "midtag,$samples[$x][$y],$mids[$x][$y]\n";
  }
}
close OUT;

print "New mapping file: $out\n";
