# README

These are some scripts that I use to perform biological analysis.

## Table Of Contents
1. [GC Calculator](#gc.py)
2. [Forking](#forking.pl)

## GC Calculator <a name = "gc.py"></a>

`gc.py` gives you node name, length, and gc content (in \%) separated by spaces.

To run this script simply pass `gc.py` your `.fasta` or `.fastq` file and a
`gc_content.csv` will be written to the working directory.

```{bash}
python3 gc.py scaffolds.fasta
```

## General Forking Script <a name = "forking.pl"></a>

This is a very simple forking script for use with a compute cluster (not that
you can't use it locally). To use the script you will have to edit some of the
values and then it should execute fine.

- `$max_cpu`
  - This value helps to set the maximum number of child process that will be
    run at any given time
- `$type`
  - This value tells the program what file type you want to use
  - As this is a simpler forking script it only currently accepts one file
    type, but can easily be expanded to multiple types
- `system ""`
  - Enter the command you would normally type into the bash shell in between
    the quotes ("")

```{perl}
my $mac_cpu = "3";
my $type = "*.gtk"
system "python3 ~/pseudofinder.py annotate -g $file -op pseudo_$file
  -db nr -t 4"
```
