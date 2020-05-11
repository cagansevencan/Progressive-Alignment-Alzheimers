# Alzheimer's Genomes Used:
APOE, PICALM, ABCA7, and TREM2

The problem we are facing is there have been many genes found to be associated with late-onset Alzheimer’s disease (AD) but not a lot of studies in all genes have been explored. Even with some identified sequences, there is still no preventative care like vaccines for this disease. Examples of genes that are associated with late-onset AD are APOE, PICALM, ABCA7, and TREM2. Using Python with a library called ‘Scikit-bio’, we have implemented the Needleman-Wunsch algorithm using dynamic programming and extrapolating to multiple sequence alignment implementing progressive pairwise alignment. Our goal is to check for distances, identities, and differences throughout the sequences that store the positions of misalignments and to determine the prominence of these genes in late-onset AD.

We have implemented ‘scikit-bio’ library tools to support our guide tree process. As the match/mismatch and gap penalty parameters were derived from NCBI BLAST, we could successfully run the algorithm with different sequences, and alignments of their products to have a prominent result at the end. Overall the steps follows as:
Take a set of query sequences-> build a guide tree-> perform progressive multiple sequence alignment, -> return the guide tree and the alignment.
• Distance matrix 
– Each pairwise alignment O(n^2)
 – Number of pairwise alignments O(k^2)
 • Iterative construction of MSA 
– Number of merge steps O(k) 
– Each pairwise alignment O(k^2n^2) 
Entire method O(k^2n^2)
Summary: 
It is also important to point out several caveats with progressive alignment heuristics:
It is not guaranteed to give the optimal MSA
Bad choice of gaps propagates
Complexity - only worse compared to few sequences:  Progressive: O(k^2 n^2) vs Dynamic Programming: O(n^k 2^k k^2)
Typically, better to merge the most closely related sequences first.

