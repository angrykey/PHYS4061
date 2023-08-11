#dot product and cross product
module vector_algebra

export dot_product, cross_product, cross_product!, dot_3d, diff_3d!

#get the norm of x-y
function difference_norm(
	x::AbstractVector{T},
	y::AbstractVector{T},
	) where {T<:Real}
	result = 0
	for i in 1:length(x)
		result += (x[i] - y[i])^2
	end
	return result
end


function dot_product(
	x::AbstractVector{T},
	y::AbstractVector{T},
	) where {T<:Real}
	return sum(x.*y)
end

function cross_product!(x::AbstractVector{T},
		y::AbstractVector{T},
		result::AbstractVector{T}
		) where {T<:Real}
	length(x) == 3 || error("x need to be legnth 3")
	length(y) == 3 || error("y need to be legnth 3")
	length(result) == 3 || error("result need to be legnth 3")

	result[1] = x[2]*y[3] - x[3]*y[2]
	result[2] = x[3]*y[1] - x[1]*y[3]
	result[3] = x[1]*y[2] - x[2]*y[1]
	return result
end
function cross_product(x::Vector{T},
        y::Vector{T},
        ) where {T<:Real}
	result = Vector{Float64}(undef, 3)
	cross_product!(x,y,result)
	return result
end

#the following 3d is used to reduce allocation thus better performance
#assum x,y only have 3 element and make no check
function dot_3d(
		x::AbstractVector{T},
		y::AbstractVector{T},
		) where {T<:Real}
	result = 0.0
	for i in 1:3
		result += x[i]*y[i]
	end
	return result
end

#assum x,y only have 3 element
function diff_3d!(
		x::AbstractVector{T},
		y::AbstractVector{T},
		result::AbstractVector{T},
		) where {T<:Real}
	for i in 1:3
		result[i] = x[i] - y[i]
	end
end
end #vector_algebra

