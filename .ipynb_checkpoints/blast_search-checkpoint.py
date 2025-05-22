
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import time

def perform_blast_search(sequence, database="nt", program="blastn"):
    """
    Perform a BLAST search on the given DNA sequence.
    
    Args:
        sequence (str): DNA sequence to search
        database (str): BLAST database to search against (default: "nt" for nucleotide)
        program (str): BLAST program to use (default: "blastn")
    
    Returns:
        list: List of BLAST hits with their details
    """
    print(f"Performing {program} search against {database} database...")
    print("This may take a few minutes...")
    
    # Perform the BLAST search
    result_handle = NCBIWWW.qblast(program, database, sequence)
    
    # Parse the results
    blast_records = NCBIXML.parse(result_handle)
    
    # Process and return the results
    hits = []
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hit = {
                    'title': alignment.title,
                    'length': alignment.length,
                    'score': hsp.score,
                    'e_value': hsp.expect,
                    'identities': hsp.identities,
                    'gaps': hsp.gaps,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'sbjct_start': hsp.sbjct_start,
                    'sbjct_end': hsp.sbjct_end
                }
                hits.append(hit)
    
    return hits

def print_results(hits):
    """
    Print the BLAST search results in a formatted way.
    
    Args:
        hits (list): List of BLAST hits
    """
    if not hits:
        print("No significant matches found.")
        return
    
    print("\nBLAST Search Results:")
    print("=" * 80)
    
    for i, hit in enumerate(hits, 1):
        print(f"\nMatch {i}:")
        print(f"Title: {hit['title']}")
        print(f"Length: {hit['length']} bp")
        print(f"Score: {hit['score']}")
        print(f"E-value: {hit['e_value']}")
        print(f"Identities: {hit['identities']}")
        print(f"Gaps: {hit['gaps']}")
        print(f"Query range: {hit['query_start']}-{hit['query_end']}")
        print(f"Subject range: {hit['sbjct_start']}-{hit['sbjct_end']}")
        print("-" * 80)

def main():
    if len(sys.argv) != 2:
        print("Usage: python blast_search.py <sequence_file>")
        print("The sequence file should be in FASTA format.")
        sys.exit(1)
    
    sequence_file = sys.argv[1]
    
    try:
        # Read the sequence from the file
        with open(sequence_file, 'r') as handle:
            record = SeqIO.read(handle, "fasta")
            sequence = str(record.seq)
        
        # Perform the BLAST search
        hits = perform_blast_search(sequence)
        
        # Print the results
        print_results(hits)
        
    except FileNotFoundError:
        print(f"Error: File '{sequence_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 