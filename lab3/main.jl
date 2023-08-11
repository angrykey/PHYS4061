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
include("./potential.jl")
using .potential

#auto output the lattice to xyz
function make_crystal(crystal, fout, cell_period, cell_size)
	position = crystal(cell_period..., cell_size=cell_size)
	atoms = lazy_Atoms(value=position)
	Atoms_to_VMD(fout, atoms)
	return atoms
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

function cal_energy(unit_vector, lattice, cut_off, pot)
	total_energy, nl_list = total_energy_2body(unit_vector, lattice, cut_off, pot)
	println("Cut_off: ", cut_off)
	println("total_energy: ", total_energy)
	println("Energy per atom: ", total_energy/size(lattice.value, 2))
	println()
end

function main3a()
	period = 10
	cell_period = (period, period, period)
	cell_size = lattice_const_Ar
	fout_dir = "./output"
	fcc_lattice = make_crystal(fcc, "$(fout_dir)/fcc.xyz", cell_period, cell_size)

	unit_const = period*cell_size
	unit_vector = rs_vector([unit_const,0.0,0], [0.0, unit_const,0], [0.0,0, unit_const])

	println("For LJ_potential")
	cal_energy(unit_vector, fcc_lattice, 2*lattice_const_Ar, LJ_potential)
	cal_energy(unit_vector, fcc_lattice, 3*lattice_const_Ar, LJ_potential)
	cal_energy(unit_vector, fcc_lattice, 4*lattice_const_Ar, LJ_potential)
	cal_energy(unit_vector, fcc_lattice, 5*lattice_const_Ar, LJ_potential)
	println()
end

function make_MgO(fcc_lattice, l_constant)
	mg_value = fcc_lattice.value
	o_value = fcc_lattice.value .+ l_constant*[0.5,0,0]
	all_value = hcat(mg_value, o_value)

	mg_type = fill("Li", size(mg_value,2))
	o_type = fill("O", size(mg_value,2))
	all_type = vcat(mg_type, o_type)

	result_atom = lazy_Atoms(value=all_value, type = all_type)

	mg_intype = fill(1, size(mg_value,2))
	o_intype = fill(2, size(mg_value,2))
	all_intype = vcat(mg_intype, o_intype)

	return result_atom, all_intype
end

function cal_buck(unit_vector, bin, cut_off)
	total_energy, nl_list = total_Buck_pot(unit_vector, bin, cut_off)
	println("total_energy: ", total_energy)
	println("Energy per basis: ", total_energy*(2/size(bin.atom_lat.value, 2)))
	println()
end

function main3b()
	period = 10
	cell_period = (period, period, period)
	cell_size = lattice_const_buck
	fout_dir = "./output"
	fcc_lattice = make_crystal(fcc, "$(fout_dir)/fcc.xyz", cell_period, cell_size)

	unit_const = period*cell_size
	unit_vector = rs_vector([unit_const,0.0,0], [0.0, unit_const,0], [0.0,0, unit_const])

	mgo_lat, mgo_type = make_MgO(fcc_lattice, cell_size)
	mgo_bin = Buck_info(mgo_lat, mgo_type, type_charge, buck_pot_A, buck_pot_B, buck_pot_C)

	println("For Buckingham potential")
	cal_buck(unit_vector, mgo_bin, 200*cell_size)
	println()
	

end

println()
@time main3a()
@time main3b()


