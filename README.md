# Ripple-MPBNN-for-FPGAs

Part of paper: Genetically Optimized Massively Parallel Binary Neural Networks for Intrusion Detection Systems, Murovič T., Trost A.

--> build_layers.m is the main script.
To change the size of hidden layer and to grab the precomputed trained network parameters the LT variable must be changed to either 64,128,192 or 256.

SEL_T = 0, creates the straight-forward combinational verilog BNN model.
SEL_T = 1, creates the ripple architecture combinational verilog BNN model.

Net_Params\ holds the traiend networks parameters.
Verilog Models\ holds the prebuild HDL models.


Networks were trained with help of: "T. Murovič, A. Trost, Massively Parallel Combinational Binary Neural Networks for Edge Processing, Elektrotehniški vestnik, vol. 86, no. 1-2, pp. 47-53, 2019".
Researchgate link: https://www.researchgate.net/publication/333563328_Massively_parallel_combinational_binary_neural_networks_for_edge_processing
Paper link: https://ev.fe.uni-lj.si/1-2-2019/Murovic.pdf
