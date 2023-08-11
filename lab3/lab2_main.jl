using Random
using Statistics

using CSV
using DataFrames
using Plots
using StatsBase

include("./phy_constant.jl")
using .phy_constant
include("./vector_algebra.jl")
using .vector_algebra
include("./VMD.jl")
using .VMD
include("./lattice.jl")
using .lattice
include("./mk_crystal.jl")
using .mk_crystal

#auto output the lattice to xyz
function make_crystal(crystal, fout, cell_period, cell_size)
	position = crystal(cell_period..., cell_size=cell_size)
	atoms = lazy_Atoms(value=position)
	Atoms_to_VMD(fout, atoms)
	return atoms
end

#for part a, check if reciprocal lattice is working as intended using FCC and BCC
function check_reciprocal(primitive, name_primitive)
	result = reciprocal(primitive, scale=1)
	println("For ", name_primitive, ",")
	println("Its primitive vector: ", primitive)
	println("Its reciprocal vector: ", result)
	println("Its unit cell volume: ", get_volume(primitive))
	println("Its reciprocal volume: ", get_volume(result))
end

function main2a()
	bc_primitive = rs_vector([1.0,0,0], [0.0,1,0], [0.0,0,1])
	fcc_primitive = rs_vector([0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5])
	bcc_primitive = rs_vector([0.5,0.5,-0.5], [-0.5,0.5,0.5], [0.5,-0.5,0.5])
	check_reciprocal(bc_primitive, "bc")
	check_reciprocal(bcc_primitive, "bcc")
	check_reciprocal(fcc_primitive, "fcc")
end

function main2b()
	cell_vector = rs_vector([5.0,0,0], [0,4.0,0], [0.0, 0.0, 3])
	x = [1.0,2,2.9]

	bs = reciprocal(cell_vector, scale=1)
	fc_before_pbc = get_fractional_coordinate(bs,x)
	fc_after_pbc = apply_pbc(fc_before_pbc)
	x_after_pbc = fc_to_real_space(cell_vector, fc_after_pbc)

	println("Original vector, x: ", x)
	println("Fractional coordinate before PBC: ", fc_before_pbc)
	println("Fractional coordinate after PBC: ", fc_after_pbc)
	println("x after PBC", x_after_pbc)
end

function export_neightbourlist(
		n_list,
		fout)
	open(fout, "w") do io
		println(io, "atom1 atom2 distance")
		for i in n_list
			println(io, i.atom1, " ", i.atom2, " ", i.distance)
		end
	end
end

function main2c()
	period = 9 
	cell_period = (period, period, period)
	cell_size = 1
	fout_dir = "./output"
	sc_lattice = make_crystal(sc, "$(fout_dir)/sc.xyz", cell_period, cell_size)
	bcc_lattic = make_crystal(bcc, "$(fout_dir)/bcc.xyz", cell_period, cell_size)
	fcc_lattice = make_crystal(fcc, "$(fout_dir)/fcc.xyz", cell_period, cell_size)

	unit_vector = rs_vector([period,0.0,0], [0.0, period,0], [0.0,0, period])
	n_list = neighbour_list(unit_vector, sc_lattice, cell_size*1.1)
	export_neightbourlist(n_list, "$(fout_dir)/sc.nearest_neighbour.txt")

	n_list = neighbour_list(unit_vector, fcc_lattice, cell_size*1.1)
	export_neightbourlist(n_list, "$(fout_dir)/fcc.nearest_neighbour.txt")

end

println("For part a, ")
@time main2a()
println()
println("For part b, ")
@time main2b()
println()
@time main2c()


