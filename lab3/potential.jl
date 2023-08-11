module potential
using ..phy_constant
using ..lattice
using ..VMD

export total_energy_2body
export LJ_potential, Buckingham_potential
export total_Buck_pot, Buck_info

function LJ_potential(
		ir::Real; 
		iepsilon::Real = epsilon, 
		isigma::Real = sigma,
		)
	#normalize r
	r = ir/sigma
	r3 = r*r*r
	r6 = r3*r3
	r12 = r6*r6
	PE = 4*epsilon*(1/r12 - 1/r6)
	return PE
end

function Buckingham_potential(
		r::Real;
		A::Real=buck_pot_A,
		B::Real=buck_pot_B,
		C::Real=buck_pot_C,
		q1::Real,
		q2::Real,
		)
	temp = 4*pi*e_permittivity
	energy = A*exp(-B*r)-C/r^6 + q1*q2/(temp*r)
	return energy
end

#info need for cal total energy using Buckingham potential
#atom_lat is position vector of lattice point
#atom_type is the type of atom on the lattice point, using index
#type_charge is the charge of atom, if for lattice 2, atom_type[2] = 1, then its charge is type_charge[atom_type[2]] =
#type-charge[1]
#cof_A is the coefficient for Buckingham potential, is a 2d array where for atom_type = 1 and 2, cof_A[1][2] = A between type 1 and 2
struct Buck_info
	atom_lat::Atoms
	atom_type::Vector{Int32}
	type_charge::Vector{Float64}
	cof_A::Matrix{Float64}
	cof_B::Matrix{Float64}
	cof_C::Matrix{Float64}
end

#cal total energy using buckingham potential
function total_Buck_pot(
		cv::rs_vector,
		bin::Buck_info,
		cut_off::Real,
		) where {T<:Real}
	nl_list = neighbour_list(cv, bin.atom_lat, cut_off)
	energy = 0
	for ineight in nl_list
		at1_type = bin.atom_type[ineight.atom1]
		at2_type = bin.atom_type[ineight.atom2]
		energy += Buckingham_potential(ineight.distance,
									   A=bin.cof_A[at1_type, at2_type],
									   B=bin.cof_B[at1_type, at2_type],
									   C=bin.cof_C[at1_type, at2_type],
									   q1=bin.type_charge[at1_type],
									   q2=bin.type_charge[at2_type],
									   )
	end
	return energy, nl_list
end

#pot is the two body potential, f(r), r is distance between the atoms
function total_energy_2body(
		cv::rs_vector,
		atom_lat::Atoms,
		cut_off::Real,
		pot::Function,
		) where {T<:Real}
	nl_list = neighbour_list(cv, atom_lat, cut_off)
	energy = 0
	for ineight in nl_list
		energy += pot(ineight.distance)
	end
	return energy, nl_list
end

end
