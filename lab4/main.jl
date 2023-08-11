using Random
using Statistics

using CSV
using DataFrames
using Plots
using StatsBase
using Accessors

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
include("./optimization.jl")
using .optimization

function get_name(func)
	return String(Symbol(func))
end

#auto output the lattice to xyz
function make_crystal(crystal, fout_dir, cell_period, cell_size)
	position = crystal(cell_period..., cell_size=cell_size)
	atoms = lazy_Atoms(value=position)
	fout = "$fout_dir/$(get_name(crystal)).xyz"
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

function get_f1_LJ(cv::rs_vector,
		in_atom_lat::Atoms,
		cut_off::Real;
		) where {T<:Real}
	f1(x) = gradient_LJ(cv, in_atom_lat, cut_off;new_pos=x)[1]
	return f1
end


function main3a()
	period = 5
	cell_period = (period, period, period)
	cell_size = lattice_const_Ar
	cell_size = sigma*2^(1/6)*1.000001
	fout_dir = "./output"
	lattice = make_crystal(sc, fout_dir, cell_period, cell_size)

	unit_const = period*cell_size
	unit_vector = rs_vector([unit_const,0.0,0], [0.0, unit_const,0], [0.0,0, unit_const])

	changed_lattice = deepcopy(lattice)
	changed_lattice.value[3,2] +=0.00001
	#changed_lattice.value[2,2] +=0.001
	#changed_lattice.value[3,2] +=0.001

	f1_LJ = get_f1_LJ(unit_vector, lattice, 5*lattice_const_Ar)

	println("For LJ_potential")

	cal_energy(unit_vector, lattice, 100*lattice_const_Ar, LJ_potential)
	old_opt_lattice = @set lattice.value = SG(f1_LJ, lattice.value, 0.0000001)
	cal_energy(unit_vector, old_opt_lattice, 100*lattice_const_Ar, LJ_potential)

	cal_energy(unit_vector, changed_lattice, 100*lattice_const_Ar, LJ_potential)
	opt_lattice = @set changed_lattice.value = SG(f1_LJ, changed_lattice.value, 0.0000001)
	cal_energy(unit_vector, opt_lattice, 100*lattice_const_Ar, LJ_potential)


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
#@time main3b()


