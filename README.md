# count-mutations

**Gregory W. Schwartz**

Takes two sequences or a  fasta file in the format of ">>Germline
header\nGermline sequence\n>Mutant header\nMutant sequence...etc" and generates
the mutation counts of the clones from each listed germline in the file per
position in a dataframe type format. Contains flags to bias for or against
silent or replacement mutations, along with only including codons with certain
number of mutations.

To install:
```
cabal update
cabal install
```

```
count-mutations, Gregory W. Schwartz

Usage: count-mutations [-t|--input-mut-type [Silent]|Replacement]
                       [-b|--input-bias [Silent]|Replacement]
                       [-c|--input-codon-mut [0]|1|2|3]
                       [-m|--input-mut-count [1]|2|3|...]
                       [-l|--input-label [count]|STRING] [-N|--remove-N]
                       [-f|--four-fold-redundant] [-I|--input STRING STRING]
                       [-i|--input-fasta FILE] [-o|--output FILE]
  Return the mutation counts with certain biases from the germline to the
  mutants within clones

Available options:
  -h,--help                Show this help text
  -t,--input-mut-type [Silent]|Replacement
                           The type of mutation to be counted
  -b,--input-bias [Silent]|Replacement
                           The type of mutation to bias for when calculating
                           intermediate codons
  -c,--input-codon-mut [0]|1|2|3
                           Only count mutations from codons with this many
                           mutations (0 is the same as include all codons)
  -m,--input-mut-count [1]|2|3|...
                           Only count a unique mutation if it appears this many
                           or more times
  -l,--input-label [count]|STRING
                           The label for the data (usually the dataset)
  -N,--remove-N            Whether to remove N or n in the sequence
  -f,--four-fold-redundant Whether to only count mutations if they are four fold
                           redundant mutations
  -I,--input STRING STRING Two sequences separated by spaces in order to find
                           the mutations between them
  -i,--input-fasta FILE    The name of the input fasta file
  -o,--output FILE         The name of the output file with the counts
```
