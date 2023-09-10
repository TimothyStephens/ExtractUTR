# ExtractUTR
Predict and extract UTRs from transcripts.

## Installation

All required tools can be installed using conda. Simply run the following commands to setup and activate the ExtractUTR environment. 
```bash
conda env create -f env.yaml
conda activate ExtractUTRs
```

## Extract UTRs

The `ExtractUTRs.sh` script needs just a fasta file with transcripts and it will produce bed and fasta files with all the 5- and 3-prime UTR sequences. These result files will have `*.ExUTRs.*` in their names. 
Optionally, a [DIAMOND](https://github.com/bbuchfink/diamond) formatted database can be provided that has a set of proteins that will be used to guide ORF prediction. Examples of what to use as the database include SwissProt or NCBI's nr protein database.
An example run might look like:
```bash
./ExtractUTRs.sh -t test.fa
```
Output files:
> test.fa.ExUTRs.3UTR.bed
>
> test.fa.ExUTRs.3UTR.fasta
>
> test.fa.ExUTRs.5UTR.bed
>
> test.fa.ExUTRs.5UTR.fasta

OR, if you want to use a set of reference proteins to guide ORF prediction:
```bash
./ExtractUTRs.sh -t test.fa -d /path/to/diamond/datbase.dmnd
```
The output files will have the same names as the previous example.


## Assemble samples from SRA

If you want to assemble transcripts using samples from SRA, this can be done with the `Assemble_Transcripts.sh` script.
This script needs, as a minimum, a list (comma separated) of 1 or more SRR IDs, and the name of the output files.

An example run might look like:
```bash
./Assemble_Transcripts.sh -s DRR330246 -o test
```

If you want to assemble multiple SRA samples together then you need to specify them as a comma separated list:
```bash
./Assemble_Transcripts.sh -s DRR330244,DRR330245,DRR330246 -o test
```

In both cases, the final transcripts will be in a file called `test.transcripts.fasta`, and the assembly stats will be in `test.transcripts.fasta.stats`.
If you are happy with the results then you can delete the directory that was created with all the temp data in it (for the above examples, this directory will be called `test`).
The `test.transcripts.fasta` file can then be used directly with the `ExtractUTRs.sh` script to extract the UTR regions.


