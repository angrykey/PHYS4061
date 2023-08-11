module lattice
using ..vector_algebra
using ..VMD

export get_volume, reciprocal, get_fractional_coordinate!, get_fractional_coordinate
export apply_pbc, apply_pbc!
export fc_to_real_space, fc_to_real_space!
export neighbour_list

#global array to reduce overhead
const temp1 = [0.0, 0, 0]
const temp2 = [0.0, 0, 0]
const default_scale = 2*pi

get_volume(a1, a2, a3) = dot_3d(a1, cross_product!(a2, a3, temp1))

#calculate reciprocal vector from primitive vector
function reciprocal(
		a1::Vector{T},
		a2::Vector{T},
		a3::Vector{T},
		;scale::Real=default_scale,
		) where {T<:Real}
	volume = get_volume(a1, a2, a3)
	b1 = (scale/volume)*cross_product(a2, a3)
	b2 = (scale/volume)*cross_product(a3, a1)
	b3 = (scale/volume)*cross_product(a1, a2)
	return (b1, b2, b3)
end

#convert vector in normal unit into fractional coordinate
#for x is a vector
function get_fractional_coordinate!(b1::Vector{T},
		b2::Vector{T},
		b3::Vector{T},
		x::AbstractVector{T},
		result::AbstractVector{T}; 
		scale=default_scale
		) where {T<:Real}
	result[1] = dot_3d(b1, x)/scale
	result[2] = dot_3d(b2, x)/scale
	result[3] = dot_3d(b3, x)/scale
	return result
end

#same as above
#for x is a matrix
function get_fractional_coordinate!(b1::AbstractVector{T},
        b2::AbstractVector{T},
        b3::AbstractVector{T},
		x::AbstractMatrix{T},
		result::AbstractMatrix{T};
		scale=default_scale
		) where {T<:Real}
	for j in 1:size(x,2)
		@views get_fractional_coordinate!(b1, b2, b3, x[1:3, j], result[1:3,j], scale = scale)
	end
	return result
end

#lazy function, for x is a vector
function get_fractional_coordinate(b1, b2, b3, x; scale=default_scale)
	result = [0.0,0.0,0.0]
	get_fractional_coordinate!(b1,b2,b3,x,result,scale=scale)
	return result
end

#convert fractional coordinate to real space
function fc_to_real_space!(b1::AbstractVector{T},
        b2::AbstractVector{T},
        b3::AbstractVector{T},
		fc_vector::AbstractVector{T},
		rs_vector::AbstractVector{T};
		scale=default_scale
		) where {T<:Real}
	for i in 1:3
		rs_vector[i] = fc_vector[1]*b1[i] + fc_vector[2]*b2[i] + fc_vector[3]*b3[i]
	end
	rs_vector *= scale
	return rs_vector
end

function fc_to_real_space(b1, b2, b3, x; scale=default_scale)
    result = [0.0,0.0,0.0]
    fc_to_real_space!(b1,b2,b3,x,result,scale=scale)
    return result
end

#fc is fractional coordinate
#u_vector is vectors specific for the whole domain
function get_fc_driver(
		a1::Vector{T},
		a2::Vector{T},
		a3::Vector{T},
		lattice::Atoms,
		) where {T<:Real}
	position = lattice.value
	fc_lattice = zeros(size(position))
	rec_vec = reciprocal(a1, a2, a3, scale=1)
	get_fractional_coordinate!(rec_vec..., position, fc_lattice, scale=1)
	return fc_lattice, rec_vec
end

struct neighbour
	atom1::Int32
	atom2::Int32
	distance::Float64
end

#get all neighbour within cut_off from the position of target
#target is the index of target atom in fc_lattice
#fc_lattice is the position vector of the whole lattice
function get_neighbour(
		cv1::Vector{T},
		cv2::Vector{T},
		cv3::Vector{T},
		fc_lattice::Matrix{T},
		target_index::Int,
		cut_off::Real
		) where {T<:Real}
	result = Vector{neighbour}()

	for i in 1:size(fc_lattice,2)
		if i == target_index
			continue
		end		
		@views diff_3d!(fc_lattice[:,i], fc_lattice[:,target_index], temp1)
		apply_pbc!(temp1, temp2)
		fc_to_real_space!(cv1, cv2, cv3, temp2, temp1, scale = 1)
		distance = dot_product(temp1, temp1)
		if distance <= cut_off
			record = neighbour(target_index, i, distance)
			push!(result, record)
		end
	end
	return result
end

#get neighbour_list from cell vector (unit cell, whole domain), lattice(position of atoms)
#and cut_off_distance
function neighbour_list(
		cv1::Vector{T},
		cv2::Vector{T},
		cv3::Vector{T},
		lattice::Atoms,
		cut_off_distance::Real,
		) where {T<:Real}
	fc_lattice, rec_vec = get_fc_driver(cv1, cv2, cv3, lattice)
	record_neighbour = []
	for i in 1:size(fc_lattice,2)
		irecord = get_neighbour(cv1, cv2, cv3, fc_lattice, i ,cut_off_distance)
		push!(record_neighbour, irecord)
	end
	result = reduce(vcat, record_neighbour)
	return result
end


function apply_1d_pbc(x)
	t1 = mod(x,1)
	if  t1 >=0.5
		result = t1-1
	elseif t1 <= -0.5
		result = t1+1
	else
		result = t1
	end
	return result
end

#x is in FRACTIONAL COORDINATE
#remove integer part and move it by 0.5
function apply_pbc!(x::Vector{T}, result::Vector{T}) where {T<:Real}
	result[1] = apply_1d_pbc(x[1])
	result[2] = apply_1d_pbc(x[2])
	result[3] = apply_1d_pbc(x[3])
end

function apply_pbc!(x::Matrix{T}, result::Matrix{T}) where {T<:Real}
	for j in 1:size(x,2)
		@views apply_pbc!(x[1:3, j], result[1:3,j])
	end
	return result
end

#lazy function 
function apply_pbc(x::Vector{T}) where{T<:Real}
	result = [0.0,0.0,0.0]
	apply_pbc!(x, result)
	return result
end

end
