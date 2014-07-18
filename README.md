gibsonsimulator
===============
This program is meant as a maximum ignorance gibson assembly predictor. Given a list of sequences, the simulator will figure out how they should be assembled.

To this end it uses a local BLAST installation (get it here ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and an implementation of the gillespie stochastic reaction simulation algorithm to actually simulate the assembly process. Eventually, you will get a sequence output (though so far that is not implemented)