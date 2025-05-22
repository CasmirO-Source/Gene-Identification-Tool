#  DNA BLAST Search Tool (Biopython)

This is a simple Python script that performs a BLAST (Basic Local Alignment Search Tool) search using a DNA sequence from a FASTA file. It uses Biopython(https://biopython.org/) to send a query to the NCBI BLAST API and parse the results.

---

##  Project Background

I first learned about bioinformatics during my biomedical science degree and became fascinated with how tools like BLAST help identify genetic sequences. I’ve continued developing my coding skills post-graduation, and this tool is part of that journey.

---

##  What This Script Does

- Reads a DNA sequence from a `.fasta` file.
- Submits the sequence to NCBI’s online BLAST service.
- Retrieves and prints matching sequences, alignment scores, identity matches, and more.

---

##  Example Files

This repo includes:
- `blast_search.py`: The main script.
- `hippo.fasta`: Example DNA sequence for testing.

---

##  Requirements

- Python 3.x
- Biopython

Install Biopython using pip:

```bash
pip install biopython
```

## How to Use

Make sure your FASTA file (like hippo.fasta) is in the same folder as the script.

Run the script in your terminal:

python blast_search.py hippo.fasta
You’ll see BLAST results printed to your terminal, including:

Match title
Score
E-value
Start/End positions of alignments
  
 ## Notes

Only one record per FASTA file is supported.
Internet connection is required (it queries NCBI online).
BLAST queries may take a few minutes depending on the sequence length.
## Future Improvements

1. Add support for multiple sequences in one file
2. Export results to a text or CSV file
3. Option to select different BLAST databases or parameters

## Resources

NCBI BLAST
Biopython Tutorial
## Author

Created by someone passionate about bioinformatics and learning by building.
Feel free to fork, contribute, or ask questions!