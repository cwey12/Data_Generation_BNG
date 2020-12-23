This is a pipeline used to:
1) Run BioNetGen simulations (on both Linux and Window's Systems). Note: Has not been tested on MacOS yet.
2) Format the data to feed into Sincerities

Running Simulations:
Run simulations using the 'Run_Simulation.py' script. It takes 1 argument, the number of cells to simulation.
Ex. To run a 1000 cell simulation, do-
python Run_Simulation 1000

The current code is written to run in Clark's HPC3 and Windows 10 directories. 
These hardcoded directories contain the location to save simulation data.
Hardcoded directories are on lines: 
-22,27 (Where to save sims)

In Progress Work:
Extract simulation data into one dataframe, and format for Junhao's data analysis.