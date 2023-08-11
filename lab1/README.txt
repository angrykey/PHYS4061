----------------lab3 README-------------------------
I use conda to install my julia compiler. 
I am using Windows Subsystem for Linux (WSL), version is "Ubuntu 20.04 LTS".

To install conda, one could look at
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Then, create new environments and install julia in that environment.

Once julia is installed, install package in julia, run the following comment
	julia -i #enter REPL

	import Pkg
	Pkg.add("Plots")
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	Pkg.add("StatsBase")
	exit() #leave REPL
------------------------------------------------------------

To run the script, use
	julia main.jl

To change the cell size and period, look at function main inside main.jl. 
There is two variables, cell_period (number of cells in each of xyz direction) and cell_size (a, size of conventional cell).
The file is output to ./output
	
