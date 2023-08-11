using Random
using Statistics

include("./phy_constant.jl")
using .phy_constant
include("./VMD.jl")
using .VMD
include("./mk_crystal.jl")
using .mk_crystal

#auto output the lattice to xyz
function make_crystal(crystal, fout, cell_period, cell_size)
	position = crystal(cell_period..., cell_size=cell_size)
	atoms = lazy_Atoms(value=position)
	Atoms_to_VMD(fout, atoms)
end

function main()
	cell_period = (4,4,4)
	cell_size = 1
	fout_dir = "./output"

	make_crystal(sc, "$(fout_dir)/sc.xyz", cell_period, cell_size)
	make_crystal(bcc, "$(fout_dir)/bcc.xyz", cell_period, cell_size)
	make_crystal(fcc, "$(fout_dir)/fcc.xyz", cell_period, cell_size)
	make_crystal(dc, "$(fout_dir)/dc.xyz", cell_period, cell_size)
end

main()

