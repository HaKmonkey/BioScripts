#!/usr/bin/env perl

use strict;
use Parallel::ForkManager;

# Input the max number of child precess you would like to run at once
my $max_cpu = "";
# Input the file type you would like to use ("*.gtk")
my $type = "";

my %data;
my @files = glob("$type");
my $pm = Parallel::ForkManager->new($max_cpu);

DATA_LOOP:
foreach my $file (@files) {
  $pm->start and next DATA_LOOP;
    # Write the command you would like to execute here!
    system "";
  $pm->finish;
}
$pm->wait_all_children;
