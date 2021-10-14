using DCAlign
using DCAbuild
using Printf


function call_dcabuild_infer_ins(fileseed::String, L::Int64, filename_ins::String)
	ctype=:amino
	seed = DCAlign.readfull(fileseed, ctype=ctype, pos = true)
	println("length(seed)=", length(seed))	
	l_o, l_e = DCAbuild.infer_ins_pen(seed, L)
	println("### Printing insertions penalties in ", filename_ins, " ###")
	file_ins = open(filename_ins, "w")
	for i in 1:L
	@printf(file_ins, "%.4f\t%.4f\n", l_o[i], l_e[i])
	end
	flush(file_ins)
end
