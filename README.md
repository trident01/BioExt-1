Some tips for working with the BioExt package. This is meant to describe a new batch of updates, NOT the package as a whole.
Email: tristan.pollner@gmail.com

IMPORTANT - A note about running the scripts - I changed their names to end with '.py'. I had to do this to import functions from a few of them, and wanted the naming to remain consistent. Run these from the command line by typing their name (ex. 'bealign.py') with necessary and optional arguments. Use the -h argument to show a list of these. Just typing 'bealign' will not work (or will not run the latest version). After typing 'bealign.py -h' you should see an option to set the extend gap penalty. 

NOTE - There are still some unobtrusive print statements. 

To setup, from the command line run "python setup.py build; python setup.py install". Any changes made to the code will be present after this setup.

**Notes about alignFinal in general:**

-alignFinal first runs bealign to align the sequences in a fasta file to a reference (if no reference is specifically selected, it will use the first sequence in the input fasta file). 
-Using this alignment, alignFinal will generate a consensus using msaConsensus (more about this below).
-Then, it bealign will run again, aligning the sequences in the input file to this consensus as reference. 

-As output, you will see two .bam files and one .fasta file. The .bam files represent the two alignments (denoted by _FIRST.bam and _FINAL.bam, where _ is replaced by the name of the input file). The .fasta file is the consensus generated. 

**Notes about msaConsensus in general:**

-This was developed because aligning to consensus often produces a nicer alignment than just aligning to reference. Simply choosing an 'A', 'C', 'G', 'T', or '-' majority at each nucleotide site does not keep frame because gaps must be stripped. 

-In an attempt to remedy this, I implemented an option to create a "multiple sequence alignment", where we can keep insertions (they were previously lost). For each sequence, we store what is matched to each position in the reference with an entry in a matrix of strings. From this, we take a consensus of each column (where each column contains what matches to a certain position in the reference). 

-In our test example, doing this and then gap stripping kept frame and created a very good alignment in the next step. I am unconvinced that this works in the general case, and more work may be needed to check this.

-Ideas about using a threshold when looking at majority gap sites (to help with detecting minority variants) and keeping insertions in groups of 3 (which are accessible from the command line) were proposed but did not correctly keep frame. They are still included, but should not be used if not updated.  


**Notes about the the aligner (bealign) in general:**

-All sequences in the input fasta file as well as the reference are gap-stripped before being aligned. 
-If no reference is selected from the option in the commandline, the first sequence in the input fasta file will be used. 
-The reference must be a multiple of 3 because we are doing a codon-based alignment which keeps frame - if this is not the case, the best option would be to add 1 or 2 N's to the end based on the error from the command line (which tells you the length of the reference mod 3). 

**Some notes about the options I added to the aligner:**

*GlobalStartingPoint:* This should be used if we know that all sequences should start at the same point at the beginning of the reference. For the default option (of false), for each pairwise alignment the aligner does not penalize for opening gaps in either the query or the reference. This led to some strange behaviors in my test data set for a few of the sequences. For more information, please refer to the help statement for the function:
                        *Sequences are penalized for not starting at the*
                        starting point of the reference (the first row and
                        column of the scoring matrix are initialized with
                        penalties). Sequences are not penalized for ending
                        early, with the caveat that at least one sequence be
                        used fully (the backtrack starts with the max score in
                        the bottom row or rightmost column of the dynamic programming
                        matrix). This option is therefore different from either a
                        full global or local alignment. *
There is as of now no option for a full global alignment. This might be of interest to add for completeness, but it will not have a large effect on the alignment. 

*The extend gap penalty:* the default of 2.5 % is based what we thought produced the nicest picture of the alignment based on a simple visualizer (which is included in this package). However, this is arbitrary, and is therefore subject to change from the command line. It would be an interesting project to find an empirical way to set the extend gap penalty, but I am not sure that it will make much of a difference to the alignment. 

*The empirical codon matrix:* 64x64 scoring matrix, with entires in the order of 'AAA' to 'TTT' (the link to the paper describing it in much more detail can be found in a comment where the matrix is initialized). This was originally included as a separate option in the command line, but now is the default. A different scoring matrix can be selected by using the --score-matrix option from the command line.
I added this in the simplest way possible by just creating a ScoreMatrix object, not using the parsing options available in the package. As far as I know of, it has all the functionality of any ScoreMatrix in the aligner. One final note - I changed the min function to not accept a value of -50 for the min if we are using a CodonScoringMatrix, as this is the score for matching a codon to a stop codon, as is not representative of the actual minimum. 




