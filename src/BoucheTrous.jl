module BoucheTrous

using SimpleAlgebra, AbstractFFTs,OptimPackNextGen, Zygote

export BoucheTrou!

function BoucheTrou!(img::AbstractArray{T,N}, badpix::BitArray{N}, pupil_radius) where {T,N}

	sz = size(img)
	pupil = sqrt.(rfftfreq(sz[1]).^2 .+ (fftfreq(sz[2]).^2)') .< pupil_radius


	F = LinOpDFT(T,sz[1:2])
	#pupil[(1:n for n ∈ outputsize(F))...]
	D = LinOpDiag(T,.!pupil)
	C = CostL2(T,size(pupil),0)

	Threads.@threads for n=1:size(img,3)
		data = img[:,:,n]
		S = LinOpSelect(T,badpix[:,:,n])
		H = C * D *F * (S' + data)
		cost(x) = H*x
		p = vmlmb( cost, S*data; lower=0., autodiff=true)
		img[:,:,n] .+=  S'*p
	end
end
end # module BoucheTrous
