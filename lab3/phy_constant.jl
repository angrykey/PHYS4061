module phy_constant

export epsilon, sigma, lattice_const_Ar
export buck_pot_A, buck_pot_B, buck_pot_C
export type_charge, lattice_const_buck
export e_permittivity

const eV_2_joules = 1.602176634e-19

#https://www.researchgate.net/publication/319412425_Phase_Change_Characteristics_of_Ultra-Thin_Liquid_Argon_Film_over_different_Flat_Substrates_at_High_Wall_Superheat_for_HydrophilicHydrophobic_Wetting_Condition_A_Non-Equilibrium_Molecular_Dynamics_Stu
const epsilon = 0.0104
const sigma = 3.4 #Å

#https://aip.scitation.org/doi/10.1063/1.1726009
const lattice_const_Ar = 5.311

#https://aip.scitation.org/doi/10.1063/1.3425844
#type1: Li
#type2: 0

#=
const buck_pot_A = [0 653.84
					653.84 0]
const buck_pot_B = [1/1 1/0.285723
				   1/0.285723 1/1]
const buck_pot_C = [0 0
					0 76.651]
const type_charge = [0.944, -0.944]
const lattice_const_buck = 4.581
=#

#=
#type1: Mg
#type2: O
const buck_pot_A = [0 929.69
					929.69 4870]
const buck_pot_B = [1/1 1/0.29909
				   1/0.29909 1/0.2670]
const buck_pot_C = [0 0
					0 77.0]
const type_charge = [1.7, -1.7]
const lattice_const_buck = 4.1
=#
const buck_pot_A = [0 821.6
					821.6 22764]
const buck_pot_B = [1/1 1/0.3242
				   1/0.3242 1/0.149]
const buck_pot_C = [0 0
					0 27.88]
const type_charge = [2, -2]
const lattice_const_buck = 4.26

const e_permittivity = 55.26349406e-4 #e2⋅eV^−1⋅Å^−1
end
