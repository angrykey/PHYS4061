module mk_crystal

export sc, bcc, fcc, dc

#simple cubic
sc(period_x::Int, period_y::Int, period_z::Int; 
			 cell_size::Real=1) = square_crystal(period_x, 
											   period_y, period_z; cell_size=cell_size, cell_atoms=[[0.0, 0.0, 0.0]])

#body centered cubic
bcc(period_x::Int, period_y::Int, period_z::Int; 
			 cell_size::Real=1) = square_crystal(period_x, 
											   period_y, period_z; cell_size=cell_size, cell_atoms=[
																									[0.0, 0.0, 0.0],
																									[0.5, 0.5, 0.5]
																									])

#face centered cubic
fcc(period_x::Int, period_y::Int, period_z::Int; 
			 cell_size::Real=1) = square_crystal(period_x, 
											   period_y, period_z; cell_size=cell_size, cell_atoms=[
																									[0.0, 0.0, 0.0],
																									[0.0, 0.5, 0.5],
																									[0.5, 0.0, 0.5],
																									[0.5, 0.5, 0.0],
																									])

#diamond cubic
dc(period_x::Int, period_y::Int, period_z::Int;
		cell_size::Real=1) = square_crystal(period_x,
													period_y, period_z; cell_size=cell_size, cell_atoms=[
																									[0.0, 0.0, 0.0],
																									[0.25, 0.25, 0.25],
																									[0.0, 0.5, 0.5],
																									[0.25, 0.75, 0.75],
																									[0.5, 0.0, 0.5],
																									[0.75, 0.25, 0.75],
																									[0.5, 0.5, 0.0],
																									[0.75, 0.75, 0.25],
																									])

#copy the square cell multiple times to get the position of atoms
function square_crystal(
		period_x::Int, 
		period_y::Int, 
		period_z::Int;
		cell_size::Real,
		cell_atoms::Vector{Vector{T}},
		) where {T<:Real}
	nb_cell_atoms = size(cell_atoms)[1]
	nb_atoms = period_x*period_y*period_z*nb_cell_atoms
	position_arr = zeros(Float64, 3, nb_atoms)
	
	for i in 1:period_x
		i_pos = (i-1)*period_y*period_z*nb_cell_atoms
		for j in 1:period_y
			ij_pos = i_pos + (j-1)*period_z*nb_cell_atoms
			for k in 1:period_z
				ijk_pos = ij_pos + (k-1)*nb_cell_atoms
				for atom_index in 1:nb_cell_atoms
					p_index = ijk_pos+atom_index
					position_arr[1,p_index] = ((i-1) + cell_atoms[atom_index][1])*cell_size
					position_arr[2,p_index] = ((j-1) + cell_atoms[atom_index][2])*cell_size
					position_arr[3,p_index] = ((k-1) + cell_atoms[atom_index][3])*cell_size
				end
			end
		end
	end
	return position_arr
end

end
