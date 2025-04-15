Alignment Experiment
====================

This assignment is designed to be completed in pairs. Choose a lab partner or
wait for the instructors to pair you up.

## The Question ##

When you examine a sequence alignment, for example, from a BLAST search, how do
you know if the alignment is any _good_? One way is to determine if the
alignment could have occurred by chance. Alignments have a _score_, so another
way to state the question is: "how often would a score of _X_ occur at random?"
This depends on several factors:

- The lengths of the sequences being aligned
- The scoring scheme used to align the sequences
	- The scores for matching and mismatching letters
	- The score for gap initiation and extension
- What exactly _random_ sequence means
	- Randomly generated with 5% each amino acid?
	- Randomly generated with typical probabilities of amino acids?
	- Randomly shuffled sequences that exist naturally?
	- Naturally occuring sequences that are unrelated?

Your goal is to design and conduct experiments that determine the relationships
between alignment scores and their probabilities of occurring by chance. To
begin, you will need some real proteins as well as programs that can generate
random sequences and perform pairwise alignments.

## Real Proteins ##

Create a directory for this project in your home directory.

```bash
mkdir ~/alignment
cd ~/alignment
```

Retrive the E. coli proteome as a compressed FASTA file. Uncompress it and name
it something simpler.

```
wget https://raw.githubusercontent.com/iankorf/E.coli/main/GCF_000005845.2_ASM584v2_protein.faa.gz
gunzip -c GCF_000005845.2_ASM584v2_protein.faa.gz > E.coli.fa
```

Examine the file with `less`. Note that the first sequence is very short and
mostly Threonine (weird).

```bash
less E.coli.fa
```

Exercise your command line skills to count how many sequences are in the file
and also count the total number of amino acids.


## Install Emboss and Blast ##

The various flavors of `conda` are a convenient way to install software. On
your personal computer, you would download and install mini-forge. Conda is
already istalled on Hive in the `module` system, so you can enable conda with
the following command:

```bash
module load conda
```

Using `nano` or some other text editor, create the following file and save it
as `alignment.yml`. YAML files are used for a variety of purposes, such as
configuration files for `conda` environments.

```
name: alignment
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - emboss
  - blast-legacy
```

Create the environment. You only need to do this once.

```
conda env create -f alignment.yml
```

Now activate the enviroment to make the programs available.

```
conda activate alignment
```

You now have access to a bunch of new programs in the EMBOSS and BLAST suites.
Note that each time you open a new terminal, you have to load the module and
activate the environment.

```
module load conda
conda activate alignment
```

The EMBOSS suite includes many CLI programs for bioinformatics tasks. The
following programs will be useful for your experiments.

- `makeprotseq` - create random protein sequences
- `pepstats` - calcuate statistics of protein sequences
- `shuffleseq` - shuffle a set of sequences, maintaining composition
- `water` - Smith-Waterman local alignment of sequences

## Random Sequences ##

Try running `makeprotseq -h`. This shows the usage statement. EMBOSS programs
are also interactive, so if you type `makeproteseq` by itself you will be
prompted for various parameters.

We'll start by making one very large protein, the reason for which will become
obvious later.

```bash
makeprotseq -amount 1 -length 99999 -outseq onebigseq.fa
```

When you are prompted for a pepstats file, hit return (we are not using one
yet). Next, run `pepstats`.

```bash
pepstats -sequence onebigseq.fa -outfile onebigseq.pepstats
```

The output of `pepstats` shows that each amino acid is produced with the same
5% probability. This is not realistic of real proteins. In order for
`makeprotseq` to generate realistic-ish sequences, it needs to be provided
typical amino acid probabilities. We can get these from the E. coli proteome.

Try running `pepstats` on the `E.coli.fa` file and examine the output file.

```bash
pepstats -sequence E.coli.fa -outfile E.coli.pepstats
```

Notice that `pepstats` reports the composition of each protein individually.
For example, the first protein is mostly rich in Threonine. We need the
composition of the whole proteome, not individual proteins. Use `grep -v` to
remove the FASTA headers from to make one HUGE E.coli protein. Then run
`pepstats` on that. At the end, save this as `E.coli.pepstats` (over-writing
the old one).

Now run `makeprotseq` to generate more realisitic proteins. Let's make 1000
proteins of average protein length (~300 aa).

```bash
makeprotseq -amount 1000 -length 300 -pepstatsfile E.coli.pepstats -outseq random-N1000-L300.fa
```

We now have 1000 sequences, each 300 aa long, whose composition is similar to
typical proteins. Verify this with your commandline skills.

## Random Alignments ##

The `water` program aligns sequences. Type `water -h` to get its usage
statement.

- `-asequence` is the "query" sequence (just one sequence)
- `-bsequence` is the "database" of multiple sequences
- `-gapopen` defaults to 10, which is okay
- `-gapextend` defaults to 0.5 which is a bit low (1 or 2 is better)
- `-datafile` defaults to BLOSUM62 (other options are BLOSUM45, BLOSUM80)

Try aligning `random-N1000-L300.fa` to itself. Note that `-asequence` only
reads the first sequence in the FASTA file, so only the first sequence from
`random-N1000-L300.fa` will get aligned to all of the other sequences in
`random-N1000-L300.fa`

```bash
water -asequence random-N1000-L300.fa -bsequence random-N1000-L300.fa -gapopen 10 -gapextend 1 -outfile test.water
```

Examine the output file with `less test.water`. Note that the first alignment
is an alignment of the first sequence to itself. It will have an absurdly high,
non-random score. You will probably want to ignore this for future
calculations.

Now `grep` the Score from the output, isolate the number, and then make a dirty
histogram with `sort` and `uniq`. Note again, the one very high score.

```bash
grep Score test.water | cut -f3 -d " " | sort -n | uniq -c
```

- What shape is this distribution?
- Specifically, does it appear normal?
- What is the modal alignment score?
- What is the average alignment score?

You now have some idea what random alignment looks like. But critical questions
remain:

- What happens if you change the length of `-asequence`?
- What happens if you change the `-gapopen` and `-gapextend` penalties?
- What happens if you change the datafile to BLOSUM45 or BLOSUM80?

Answer these questions by doing experiments with `makeprotseq` and `water`.

-----------------------------------------------------------------------------

Using R, graph average alignment score as a function of sequence length (keep
scoring parameters constant). You will want to make a shell script that does
the following:

- Makes sequences of some length
- Aligns sequences with `water`
- Extracts alignment scores
- Saves to a file

Hint: use a for loop

## Biological Alignment ##

Retrieve the second protein from the E. coli proteome using `head` and `tail`
(the first one is too short to be useful).

```bash
head -14 E.coli.fa | tail -12 > NP_414543.1.fa
```

Search this against the whole E. coli proteome. This will take a little longer
(~20 sec) than the previous search because NP_414543.1.fa is 820 aa and the E.
coli proteome is larger than our random sequences.

```bash
water -asequence NP_414543.1.fa -bsequence E.coli.fa -gapopen 10 -gapextend 1 -outfile NP_414543.1.align
```

Use `less` to examine the alignments by eye. Note that the second alignment is
the sequence to itself, so the percent identity is 100% and the score is huge
(4141.0).

```bash
grep Score NP_414543.1.align | cut -f3 -d " " | sort -n | uniq -c
```

Unlike the search against random proteins, some of these alignments have very
high scores. The highest scoring alignment (not to itself) is 915.0. Which
sequence is this? Use `less` again, and use the search function (the slash key)
to find "Score: 915".

Now you can get the name of the matching sequence. Open the E. coli fasta file
and search for that name. Then copy-paste the sequence into a new file.

To verify you did things correctly, try running `water` again on just these two
sequences.

```bash
water -asequence NP_414543.1.fa -bsequence NP_418375.1.fa -gapopen 10 -gapextend 1 -outfile test.pair
```

You should see a score of 915.0 again. Is a score of 915 _good_? It seems it
would be hard to achieve that by chance. But how do we know? When in doubt, a
shuffling experiment is generally a good way to separate signal from random
noise.

The `shuffleseq` program scrambles the letters of a sequence, keeping all of
the amino acid counts the same. Let's use that to create 1000 scrambled
sequences and then examine the score distribution.

```bash
shuffleseq -sequence NP_418375.1.fa -shuffle 1000 -outseq shuff.fa
water -asequence NP_414543.1.fa -bsequence shuff.fa -gapopen 10 -gapextend 1 -outfile shuff.align
grep Score shuff.align | cut -f3 -d " " | sort -n | uniq -c
```

-----------------------------------------------------------------------------

- Compute the mean and standard deviation of the alignment scores
- How many z-scores away from the mean is 915?
- Do you think a score of 915 is significant?
- What score threshold is the typical 5% cutoff?
- Do you think that cutoff is useful here? Explain.
- Should you correct of for multiple hypotheses?
- If so, what is the corrected P-value threshold?
- What threshold score would **you** use for "biologically significant?"

Graph the original alignment scores and shuffled alignment scores in R.

## BLAST ##

BLAST has several advantages over the Smith-Waterman algorithm in `water`.

- Faster
- Provides statistical significance
- Finds suboptimal alignments

To run BLAST, we must first create a database from the E. coli fasta file. The
usage statement shows with `formatdb --help`.

```bash
formatdb -i E.coli.fa
```

Several files are created by the command. In addition to files with `.phr`,
`.pin`, and `.psq` extensions, there is a `formatdb.log`. The log file can be
deleted, but the others comprise the blast database.

To align proteins to each other, we use BLASTP. To see the incredibly long
usage statement for blast, type `blastall --help`. Our command is pretty
simple.

```bash
blastall -p blastp -d E.coli.fa -i NP_414543.1.fa > NP_414543.1.blastp
```

Examine the output file `less NP_414543.1.blastp` and note that the first
alignment is the sequence to itself. Unlike the `water` output, BLAST output is
sorted by score. So we only need to scroll a little way to get to the 2nd best
alignment.

The Score here is reported as 333 bits (853). The score of 853 corresponds to
the previous score of 915. Why is the score not 915? Even though BLASTP and
`water` use the same scoring matrix by default, BLOSUM62, BLASTP uses slightly
different gapping penalties (11, 1 instead of 10, 1). Also, it does some low
complexity masking, which you can see as XXXXXXX in the query sequence.

The "normalized score" of 333 bits is what is used for the calculation of the
E-value: e-101. Small E-values are equivalent to P-values, so the P-value of
e-101 is much below the "usual" threshold 0.05. Clearly, an alignmet this good
cannot happen by chance (according to the BLAST probability model).

Unlike `water`, BLAST will preform alignments on all of the input sequences,
not just the first one. You can align all of the random sequences to E. coli as
follows:

```bash
blastall -p blastp -d E.coli.fa -i random-N1000-L300.fa > E.coli-vs-random.blastp
```

Use `grep`, `cut`, `sort`, and `uniq` to observe the scores in the BLAST
report.

-----------------------------------------------------------------------------

Using the `time` command, compare the speed of `water` and BLASTP in several
different scenarios. How much faster is BLAST?

Make graphs comparing `water` and BLASTP for random and shuffled sequences.

## Posters ##

Summarize all of your knowledge and findings about pairwise alignment in 2
posters. One poster should have a traditional format (lots of detail) while the
other should be more modern (more eye-catching and less detail).
