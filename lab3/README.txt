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

For part a), the coefficient for LJ potential are epsilon = 0.0104, sigma = 3.4Å, which
is from [1]. The experimental value for total energy is -0.080 eV per atom [2]. The calculated is -0.088eV, with 10% difference compared to experiment.

For part b, the coefficient for Buckingham potential is from [3]. The experimental value is calculated from lattice energy, −3795 kJ/mol,
equal to 39.3324 eV per basis (one atom of Mg and one atom of O), The calculated value is -40.5eV, with 3% difference.


[1] https://www.researchgate.net/publication/319412425_Phase_Change_Characteristics_of_Ultra-Thin_Liquid_Argon_Film_over_different_Flat_Substrates_at_High_Wall_Superheat_for_HydrophilicHydrophobic_Wetting_Condition_A_Non-Equilibrium_Molecular_Dynamics_Stu
[2] C. Kittel, Introduction to Solid State Physics, P.50
[3] https://journals.aps.org/prb/pdf/10.1103/PhysRevB.72.115437
