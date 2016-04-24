Experiments to understand how assemblers handle repeat regions.

Examples:

1. `sim.pl data/AB-AIB.fa`

Simulate reads from the `AB-AIB` FASTA file and assemble them with velvet and SPAdes.

2. `sim.pl -e 0 -c 30 -a data/AB-AIB.fa`

Simulate error free reads at 30X coverage and plot all-to-all matches in the 'reference vs assembly' dot plot.

```
Usage: sim.pl [options] ref.fa

       -h          - print this help statement
       -a          - show all to all matches in reference vs assembly dot plot
       -c cov      - fold coverage used in simulating reads (D: 20)
       -e fract    - fraction of standard Illumina HiSeq 2000 error rate (D: 1.0)
       -l len      - simulated read length (D: 100)
       -m mean     - mean of paired end read fragment size (D: 300)
       -s stdev    - standard deviation of read fragment size (D: 10)

```

