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

For part a, the volumn of unit cell and reciprocal vector is printed to the terminal.
As expected, the primitive vector and reciprocal vector is realted as 
			sc <=> sc
			fcc <=> bcc.

For part b, if you need to change the distance vector (x) or the unit vector (cell_vector)
go to main2b and change the corresponding variable.

For part c, I used a bc lattice with conventional cell length, a = 1, cut-off distance=1.1.
Each should give 6 nightbour. The output text file is in ./output/sc.nearest_neighbour.txt.


