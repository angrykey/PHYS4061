#for VMD io
module VMD
export Atoms, VMD_to_atoms, empty_Atoms, lazy_Atoms, Atoms_to_VMD, Atoms_to_log
export VMD_get_first_record

#value correspond to VMD position, but since Atoms will also used for other
#quantity, value is a suitable name
struct Atoms
	value::Array{Float64,2}
	type::Vector{String}
	time::Float64

	Atoms(value, type, time) = begin
		if (size(value,1) == 3)
			size(value,2) == length(type) ?
			new(value, type, time) :
			throw(
			BoundsError("dimention 2 of value need to be the same as type")
			)
		else
			throw(
			BoundsError("dimention 1 of value need to have length 3")
			)
		end
	end
end

#allow the use of size(x), x is Atoms class
#return the number of atoms in x
function Base.size(atoms::Atoms)
	return size(atoms.type)
end

#allow the use of length(x), x is Atoms class
#same as size(x)
function Base.length(atoms::Atoms)
	return size(atoms)[1]
end

#initial an Atoms class with all default value
function empty_Atoms(number_of_atoms::Int, 
		time::Real; 
		default_type::String="H"
		) where {T<:Real}
	value = zeros(Float64, 3, number_of_atoms)
	type = fill(default_type, number_of_atoms)
	result = Atoms(value, type, time)
	return result
end

#initialize Atoms, with auto fill
function lazy_Atoms(;value::Union{Nothing,Array{T,2}} = nothing,
		type::Union{Nothing, Vector{String}} = nothing,
		number_of_atoms::Union{Nothing, Int} = nothing,
		time::Real=0
		) where {T<:Real}
	if number_of_atoms isa Nothing
		if !(value isa Nothing) 
			number_of_atoms = size(value)[2]
		elseif !(type isa Nothing)
			number_of_atoms = length(type)
		else
			throw(
				  error("Atleast one of value, type or number_of_atoms need be specify")
				  )
		end
	end
	if value isa Nothing
		value = zeros(Float64, 3, number_of_atoms)
	end
	if type isa Nothing
		type = fill("H", number_of_atoms)
	end
	return Atoms(value, type, time)
end

#output animation as .xyz
function Atoms_to_VMD(fout::String, models::Vector{Atoms})
	open(fout, "w") do io
		for model in models
			Atoms_to_log(io, model)
		end
	end
end

#output one slide to .xyz
function Atoms_to_VMD(fout::String, atoms::Atoms)
	Atoms_to_VMD(fout, [atoms])
end

#for inner loop when output to .xyz
function Atoms_to_log(io::IOStream, model::Atoms)
	println(io, size(model.value,2))
	println(io, "time=", model.time)
	for (atom, type) in zip(eachcol(model.value), model.type)
		println(io, type, " ", atom[1], " ", atom[2], " ", atom[3])
	end
end

#for reading a .xyz line
function read_VMD_line!(s::String, atoms::Atoms, i::Int64)
	cpf(x) = parse(Float64, x)

	tp = split(s, " ")
	type = tp[1]
	p1 = cpf(tp[2]) 
	p2 = cpf(tp[3]) 
	p3 = cpf(tp[4]) 

	atoms.value[:,i] = [p1, p2, p3]
	atoms.type[i] = type
end

#for reading .xyz, output as atoms
function VMD_to_atoms(io::IOStream)
	number_of_atom::Int64 = parse(Int64, readline(io))
	ts = split(readline(io), "=")[2]
	time::Float64 = parse(Float64, ts)

	atoms = empty_Atoms(number_of_atom, time)
	for i in 1:number_of_atom
		read_VMD_line!(readline(io), atoms, i)
	end
	return atoms
end

#read the first atom of a animation
function VMD_get_first_record(fin::String)
	local atoms
	open(fin, "r") do io
		atoms = VMD_to_atoms(io)
	end
	return atoms
end


end




