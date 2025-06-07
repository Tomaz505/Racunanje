function Lagrange(xi::Vector{Number},yi::Vectot{Number})

	n = length(xi)
	
	if n< length(yi)
		error("yi prekratek")
	end


	f(x::Number) = sum([ yi[i]/prod(xi[setdiff(1:n,i)] .- xi[i] )*prod(xi[setdiff(1:n,i)].-x) for i = 1:n])
	
	f(x::Vector{Number}) = [ f(x::Number) for x =x ]

	return f
end
