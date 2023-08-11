module optimization
using ..lattice

const secant_allow_negative = true

export SG, dot
export secant
export mk_tmp

struct Temp_sheet{T<:Real}
	x_i::Matrix{T}
	x_next::Matrix{T}
	g_i::Matrix{T}
	x_diff::Matrix{T}
end

function mk_tmp(x)
	x_i = deepcopy(x)
	x_next = deepcopy(x)
	g_i = deepcopy(x)
	x_diff = deepcopy(x)
	tmp = Temp_sheet(x_i, x_next, g_i, x_diff)
	return tmp
end

function SG(
		f1::Function,
		ini_val::Matrix{T},
		threshold::Real;
		max_i=1000,
		) where {T<:Real}
	n = length(ini_val)
	x_i = deepcopy(ini_val)
	for i in 1:1000
		g_i = -f1(x_i)
		if dot(g_i,g_i)/n <= threshold*1e-15
			println("SG: g_i get to threshold")
			break
		end
		x_i = secant(f1, x_i, 0.1, g_i, 0.001)
	end
	return x_i
end

function dot(x::Matrix, y::Matrix)
	result = 0.0
	for j in 1:size(x,2)
		for i in 1:size(x,1)
			result += x[i,j]*y[i,j]
		end
	end
	return result
end

function secant(f1::Function,
		x1::Matrix{T},
		dx0::Real,
		g::Matrix{T},
		threshold::Real;
		max_i = 1000
		) where {T<:Real}
	alpha_i_n1 = 0.0
	alpha_i = dx0
	x_i_n1 = deepcopy(x1)
	local dx
	pref = 0
	for i in 1:max_i
		dx = alpha_i - alpha_i_n1
		x_i = x1+alpha_i*g
		if (!secant_allow_negative) & (alpha_i<0 || isnan(dx))
			println("i:", i)
			#println("x_i",x_i)
			#println("x_i_n1", x_i_n1)
			println("alpha_i:", alpha_i)
			#println("g:",g)
			if  alpha_i<0
				println("alpha_i<0 which indicate the function is exiting a local minimum, there maybe no minmum")
			elseif isnan(dx)
				println("dx is nan")
			end
			exit()
		end
		if abs(dx) <= threshold
			break
		end
		f1_i_n1 = f1(x_i_n1)
		f1_i = f1(x_i)
		pref = dot(f1_i, g)/dot(f1_i-f1_i_n1, g)
		alpha_i_p1 = alpha_i - pref*dx

		alpha_i_n1 = alpha_i
		alpha_i = alpha_i_p1
		x_i_n1 = x_i
	end
	return x_i_n1
end


end
