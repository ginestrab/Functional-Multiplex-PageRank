%++++++++ Functional Multiplex PageRank ++++++++++++++++++++++++++++

This folder contains 4 MATLAB codes for calculating the 
Functional Multiplex PageRank of duplex networks and of 
networks with arbitrary number of layers:

1) functionalPageRank_duplex.m

	calculating the Functional Multiplex PageRank of a duplex 
	network given the influence vector z=[z^(1,0),z^(0,1),z^(1,1)].

2) fPR.m

	calculating the Functional Multiplex PageRank of a duplex network 
	for all values of the influence vector.
	This code makes use of the code functionalPageRank_duplex.m
	
3) functionalPageRank_multiplicity.m

	calculating the Functional Multiplex PageRank of a multiplex network 
	with arbitrary number of layers and with a specific influence vector 
	depending only of the multiplicity of the overlap of the links.

4) fPRm.m

	calculating the Functional Multiplex PageRank of a multiplex network 
	with arbitrary number of layers with varying influence parameters.
	This code makes use of the code functionalPageRank_multiplicity.m


These programs are distributed ny the authors in the hope that it will be 
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
 If you use any of these codes please cite 

 [1]  J. Iacovacci, C. Rahmede, A. Arenas and G. Bianconi, "Functional Multiplex PageRank." 
        arxiv:1608.06328 (2016)  

 (c) Jacopo Iacovacci (mriacovacci@hotmail.it) 
     Christoph Rahmede (c.rahmede@kit.edu)
     Ginestra Bianconi (g.bianconi@qmul.ac.uk)  
