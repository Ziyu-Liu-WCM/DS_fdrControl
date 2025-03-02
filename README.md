## ScienceData_benchmark

This folder holds the simulation framework using the data from the Science paper. You can use the Load_Simulator.R to generate simulation data and use Simulation_Main.R for analysis, by manually setting parameters of interest.


## maaslin3_benchmark

This folder holds the simulation framework from the MaAsLin3 paper https://pmc.ncbi.nlm.nih.gov/articles/PMC11661281/#SD1. Though it's still in preprint, but it's said to have been accepted by Nature Methods, so I think their framework must make some sense.

They used the SparseDOSSA2 simulator to generate log-normal microbiome simulation datasets. Please refer to general_evaluations/data_generation/SD2.R for generating data(which I've already done) and general_evaluations/run_tools/run_DSBin.R for analysis. After analysis, the metric evaluation tool is provided at general_evaluations/evaluate_tools/evaluate_tools.R.


## The two simulators both produced FDR(1 - precision) around 0.5