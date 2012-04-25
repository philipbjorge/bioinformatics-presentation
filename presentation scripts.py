### SEARCHING
from Bio import Entrez
Entrez.email = "bjorgep@students.wwu.edu"
search_handle = Entrez.esearch(db="nucleotide", term="Opuntia[orgn] and rpl16", usehistory="y")
opuntia_rpl16_search_results = Entrez.read(search_handle)
search_handle.close()
# Look at opuntia_rpl16_search_results.keys()


### RETRIEVING
file = open("opuntia_rpl16.fasta", "w")
opuntia_rpl16_fasta_data = Entrez.efetch(db="nucleotide",rettype="fasta", retmode="text", webenv=opuntia_rpl16_search_results["WebEnv"], query_key=opuntia_rpl16_search_results["QueryKey"]).read()
file.write(opuntia_rpl16_fasta_data)
# Look at data file



### FUNCTIONS
from Bio import Entrez
Entrez.email = "bjorgep@students.wwu.edu"

def search_and_retrieve_fasta(database, searchterm, filename, batch_size=3):
  # SEARCH
  search_results = Entrez.read( Entrez.esearch(db=database, term=searchterm, usehistory="y") )
  # RETRIEVE
  file = open(filename, "w")
  count = int(search_results['Count'])
  for start in range(0,count,batch_size):
    end = min(count, start+batch_size)
    print "Going to download record %i to %i" % (start+1, end)
    data = Entrez.efetch(db=database,rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=search_results["WebEnv"], query_key=search_results["QueryKey"]).read()
    file.write(data)

# You can now do all these steps in just one line!
search_and_retrieve_fasta("nucleotide", "Blossfeldia[orgn] and rpl16", "blossfeldia_rpl16.fasta")
# Look at output file

### PARSING
from Bio import SeqIO
blossfeldia_rpl16_sequences = list(SeqIO.parse("blossfeldia_rpl16.fasta", "fasta"))
# Look at sequence list and blossfeldia_rpl16_sequences[0]

### SEQUENCE OBJECTS
_first_blossfeldia_rpl16_sequence = blossfeldia_rpl16_sequences[0].seq
first_blossfeldia_rpl16_sequence = blossfeldia_rpl16_sequences[0].seq
# GC %
from Bio.SeqUtils import GC
GC(first_blossfeldia_rpl16_sequence)

# DNA --> RNA --> DNA
first_blossfeldia_rpl16_sequence.transcribe()
first_blossfeldia_rpl16_sequence.back_transcribe()

# DNA Coding Strand --> Protein
first_blossfeldia_rpl16_sequence.translate()


### BLASTING
from Bio.Blast import NCBIWWW
from Bio import SeqIO
result_handle = NCBIWWW.qblast("blastn", "nt", _first_blossfeldia_rpl16_sequence)
save_file = open("blast_search_on_first_blossfeldia_rpl16_sequence.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()


### Clustal
import os
from Bio.Align.Applications import ClustalwCommandline

clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
clustalw_cmd_line = ClustalwCommandline(clustalw_exe, infile="opuntia_rpl16.fasta")
stdout, stderr = clustalw_cmd_line()   #outputs two files opuntia_rpl16.aln, opuntia_rpl16.dnd

# Read Multiple Alignment
from Bio import AlignIO
opuntia_rpl16_alignment = AlignIO.read("opuntia_rpl16.aln", "clustal")
print opuntia_rpl16_alignment

# Draw Phylo Tree
from Bio import Phylo
opuntia_rpl16_tree = Phylo.read("opuntia_rpl16.dnd", "newick")
Phylo.draw(opuntia_rpl16_tree)