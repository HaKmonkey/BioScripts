#!/usr/bin/env perl

use strict;
use Parallel::ForkManager;

my %data;
my @files = glob("*.pdf");
my $pm = Parallel::ForkManager->new(3);

DATA_LOOP:
foreach my $file (@files) {
  $pm->start and next DATA_LOOP;
    # system
  $pm->finish;
}
$pm->wait_all_children;
