# Gene Identification Tool

A command-line Python script that takes an unknown DNA sequence, sends it to the NCBI BLAST web service and prints back the closest matches from its database. BLAST is the standard tool for working out what a sequence is by comparing it against millions of known ones.

I first came across BLAST during my biomedical science degree. I built this afterwards to learn how to call a web service from code and turn the response into something readable.

## What it does

1. Reads a DNA sequence from a FASTA file (the standard text format for sequences).
2. Submits it to NCBI BLAST using Biopython.
3. Parses the XML response into a list of matches.
4. Prints each match with its name, length, score, E-value (how likely the match is to be chance, lower is better), identical positions, gaps and where it aligns.

If the file is missing or the request fails, it exits with a clear error message rather than a stack trace.

## Example

`Hippo.fasta` holds a 600-letter section of the hippopotamus COX1 gene, a gene commonly used to identify species. It's a published record (PQ902075.1), which makes it a good test that the script is working.

```bash
pip install -r requirements.txt
python blast_search.py Hippo.fasta
```

Needs an internet connection. NCBI can take a few minutes to respond depending on load.

## Limitations

- Only accepts files with a single sequence
- Results print to the screen and aren't saved
- The database and search type are fixed in the code

## What I'd do next

- Support files with several sequences
- Save results to CSV so they can be loaded elsewhere
- Add command-line options for the database and search type
