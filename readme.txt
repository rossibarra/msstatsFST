A modification of  K. Thornton’s msstats program (http://molpopgen.org/software/lseqsoftware.html) that calculates Weir & Cockerham’s Fst estimator from ms simulations.

You must have libsequence installed for this program to function.  

To use msstatsFST with Hudson's ms, try this example:

ms 10 1 -t 2 -I 2 6 4 0 -ej 4 2 1 | msstats -I 2 6 4

which will calculate Fst between two simulated populations of size 6 and 4, respectively.  
It should output two lines of stats, one for each pop.  It outputs values of "nan" for all but the final population, which shows the Fst between the two populations.  It also works for more than 2 populations, but I cannot remember if  the correctness of the calculations for >2 pops has been checked.  Please email me for questions or bug reports.  
If you use this code, please cite the below papers:

Eckert, A.J., J. van Heerwaarden, J.L. Wegrzyn, C.D. Nelson, J. Ross-Ibarra, S.C. González-Martínez, and D.B. Neale. 2010. Patterns of population structure and environmental associations with aridity across the range of loblolly pine (Pinus taeda L, Pinaceae). Genetics 185: 969-982
Thornton K (2003) libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19: 2325-2327.

Note that this is NOT an up-to-date version of msstats.  Please check the most recent version for msstats for additional changes/bug fixes at the website above.
