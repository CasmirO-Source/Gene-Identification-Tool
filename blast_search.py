# This is a script to do a BLAST search for DNA sequences
# I first learnt about bioinformatics in my degree in biomedical science and hope to pursue it as a career. I continued my practice over studying and took a year after graduating to improve my coding skills.
# I’m using Biopython cause it makes this easier

from Bio import SeqIO # for reading fasta files
from Bio.Blast import NCBIWWW # to connect to BLAST online
from Bio.Blast import NCBIXML # to read the results
import sys # for command line stuff
import time # to wait for BLAST

def search_with_blast(my_dna_sequence, database="nt", search_type="blastn"):
    """
    this function sends a DNA sequence to BLAST to find matches
    args:
        my_dna_sequence: the dna sequence like ATCGATCG
        database: the database to search in, default is nt (big nucleotide one)
        search_type: what kind of blast, default is blastn
    returns:
        a list of matches with all the details
    """
    print("Starting BLAST search now...")
    print("This might take a while cause it’s online!") # added this to let user know
    
    # doing the BLAST search
    result = NCBIWWW.qblast(search_type, database, my_dna_sequence)
    
    # reading results
    blast_data = NCBIXML.parse(result)
    
    # making a list to store matches
    matches = [] # empty list to start
    for record in blast_data:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # saving all the details
                match_info = {
                    'name': alignment.title, # name of sequence
                    'length': alignment.length, # how long it is
                    'score': hsp.score, # how good is the match
                    'e_value': hsp.expect, # how likely by chance
                    'matches': hsp.identities, # exact matches
                    'gaps': hsp.gaps, # any gaps
                    'query_start': hsp.query_start, # start of my sequence
                    'query_end': hsp.query_end, # end of my sequence
                    'subject_start': hsp.sbjct_start, # start in database
                    'subject_end': hsp.sbjct_end # end in database
                }
                matches.append(match_info) # add to list
                print(f"Found a match: {alignment.title}") # extra print to see progress
    
    return matches # return all matches

def show_results(matches_found):
    """
    shows the blast results nicely
    args:
        matches_found: list of matches from blast
    """
    if len(matches_found) == 0: # check if no matches
        print("No matches found, sorry!")
        return
    
    print("\nHere are your BLAST results!!")
    print("=" * 50) # cool line to separate
    
    for i in range(len(matches_found)): # loop through matches
        match = matches_found[i]
        print(f"\nMatch number {i + 1}:") # number each match
        print(f"Name: {match['name']}")
        print(f"Length: {match['length']} base pairs") # added base pairs
        print(f"Score: {match['score']} (bigger is better)") # explain score
        print(f"E-value: {match['e_value']} (smaller is better)") # explain e-value
        print(f"Exact matches: {match['matches']} letters match")
        print(f"Gaps: {match['gaps']} gaps")
        print(f"Your sequence: {match['query_start']} to {match['query_end']}")
        print(f"Database sequence: {match['subject_start']} to {match['subject_end']}")
        print("-" * 50) # another line

def main():
    # check if user gave a file
    if len(sys.argv) != 2:
        print("You need to give a fasta file!")
        print("Like this: python blast_search.py my_dna.fasta")
        sys.exit(1) # exit if wrong
    
    my_sequence_file = sys.argv[1] # get file name
    
    try:
        # read the fasta file
        with open(my_sequence_file, 'r') as file:
            record = SeqIO.read(file, "fasta") # read fasta
            dna_sequence = str(record.seq) # get the sequence
            print(f"Read sequence: {dna_sequence[:10]}...") # show first 10 letters
        
        # do blast search
        print("Starting search...")
        matches = search_with_blast(dna_sequence)
        
        # show results
        show_results(matches)
        
    except FileNotFoundError:
        print(f"Error: Can’t find the file '{my_sequence_file}'!")
        sys.exit(1)
    except Exception as error:
        print(f"Something broke: {error}") # general error
        sys.exit(1)

if __name__ == "__main__":
    main() # run the program
