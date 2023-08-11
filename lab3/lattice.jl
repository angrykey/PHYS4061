module lattice
using ..vector_algebra
using ..VMD

export get_volume, reciprocal, get_fractional_coordinate!, get_fractional_coordinate
export apply_pbc, apply_pbc!
export fc_to_real_space, fc_to_real_space!
export neighbour_list
export rs_vector, rec_vector

#global array to reduce overhead
const temp1 = [0.0, 0, 0]
const temp2 = [0.0, 0, 0]
const default_scale = 2*pi

#unit cell vector
struct rs_vector{T<:Real}
a1::Vector{T}
a2::Vector{T}
a3::Vector{T}
end

#reciprocal vector
struct rec_vector{T<:Real}
b1::Vector{T}
b2::Vector{T}
b3::Vector{T}
scale::T
end

get_volume(a1, a2, a3) = dot_3d(a1, cross_product!(a2, a3, temp1))
get_volume(rc::rec_vector) = dot_3d(rc.b1, cross_product!(rc.b2, rc.b3, temp1))
get_volume(rs::rs_vector) = dot_3d(rs.a1, cross_product!(rs.a2, rs.a3, temp1))

#calculate reciprocal vector from primitive vector
function reciprocal(
		rs::rs_vector{T},
		;scale::Real=default_scale,
		) where {T<:Real}
	volume = get_volume(rs.a1, rs.a2, rs.a3)
	b1 = (scale/volume)*cross_product(rs.a2, rs.a3)
	b2 = (scale/volume)*cross_product(rs.a3, rs.a1)
	b3 = (scale/volume)*cross_product(rs.a1, rs.a2)
	rec_vec = rec_vector(b1, b2, b3, T(scale))
	return rec_vec
end

#convert vector in normal unit into fractional coordinate
#for x is a vector
function get_fractional_coordinate!(
		rc::rec_vector{T},
		x::AbstractVector{T},
		result::AbstractVector{T}
		) where {T<:Real}
	result[1] = dot_3d(rc.b1, x)/rc.scale
	result[2] = dot_3d(rc.b2, x)/rc.scale
	result[3] = dot_3d(rc.b3, x)/rc.scale
	return result
end

#same as above
#for x is a matrix
function get_fractional_coordinate!(
		rc::rec_vector{T},
		x::AbstractMatrix{T},
		result::AbstractMatrix{T}
		) where {T<:Real}
	for j in 1:size(x,2)
		@views get_fractional_coordinate!(rc, x[1:3, j], result[1:3,j])
	end
	return result
end

#lazy function, for x is a vector
function get_fractional_coordinate(rc, x)
	result = [0.0,0.0,0.0]
	get_fractional_coordinate!(rc,x,result)
	return result
end

#convert fractional coordinate to real space
function fc_to_real_space!(
		rs::rs_vector{T},
		fc_vector::AbstractVector{T},
		converted_vector::AbstractVector{T},
		) where {T<:Real}
	for i in 1:3
		converted_vector[i] = fc_vector[1]*rs.a1[i] + fc_vector[2]*rs.a2[i] + fc_vector[3]*rs.a3[i]
	end
	return rs_vector
end

function fc_to_real_space(rc, x)
    result = [0.0,0.0,0.0]
    fc_to_real_space!(rc,x,result)
    return result
end

#fc is fractional coordinate
#u_vector is vectors specific for the whole domain
function get_fc_driver(
		rs::rs_vector{T},
		lattice::Atoms;
		scale=default_scale,
		) where {T<:Real}
	position = lattice.value
	fc_lattice = zeros(size(position))
	rec_vec = reciprocal(rs)
	get_fractional_coordinate!(rec_vec, position, fc_lattice)
	return fc_lattice, rec_vec
end

struct Neighbour
	atom1::Int32
	atom2::Int32
	distance::Float64
end

#get all neighbour within cut_off from the position of target
#fc_lattice is the position vector of the whole lattice
function get_neighbour!(
		cv::rs_vector{T},
		fc_lattice::Matrix{T},
		cut_off::Real,
		neighbour_list::Vector{Neighbour}
		) where {T<:Real}
	cut_off_square = cut_off*cut_off
	for i in 1:size(fc_lattice,2)
		for j in (i+1):size(fc_lattice,2)
			@views diff_3d!(fc_lattice[:,i], fc_lattice[:,j], temp1)
			apply_pbc!(temp1, temp2)
			fc_to_real_space!(cv, temp2, temp1)
			distance_square = dot_product(temp1, temp1)
			if distance_square <= cut_off_square
				distance = sqrt(distance_square)
				record = Neighbour(i, j, distance)
				push!(neighbour_list, record)
			end
		end
	end
	return neighbour_list
end

#get neighbour_list from cell vector (unit cell, whole domain), lattice(position of atoms)
#and cut_off_distance
function neighbour_list(
		cv::rs_vector{T},
		lattice::Atoms,
		cut_off_distance::Real,
		) where {T<:Real}
	fc_lattice, rec_vec = get_fc_driver(cv, lattice)
	record_neighbour::Vector{Neighbour} = []
	get_neighbour!(cv, fc_lattice ,cut_off_distance, record_neighbour)
	return record_neighbour
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
