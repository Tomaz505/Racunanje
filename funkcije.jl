using LinearAlgebra
"""Lagrange(xi::vector{<:Number},yi::Vector{<:Number})

	f = Lagrange(xi,yi)
	f(xi) = yi
	
	f(x::Number) & f(x::Vector{<:Number})

	
"""
function Lagrange(xi::Vector{<:Number},yi::Vector{<:Number})

	n = length(xi)
	
	if n< length(yi)
		error("yi prekratek")
	end

	ci = [ yi[i]/prod(xi[setdiff(1:n,i)].-xi[i]) for i=1:n]


	f(x::Number) = sum([ ci[i] * prod(xi[setdiff(1:n,i)].-x) for i = 1:n])
	
	f(x::Vector{<:Number}) = [ f(x::Number) for x =x ]

	return f
end



"""LeastSquares(xi::vector{<:Number},yi::Vector{<:Number};n = 1)

	f = LeastSquares(xi,yi) = sum( ai * x^i )

	f(x::Number) & f(x::Vector{<:Number})
"""
function LeastSquares(xi::Vector{<:Number},yi::Vector{<:Number};n = 1)
	A = hcat([[ sum(xi.^(i+j)) for i=0:n ] for j=0:n] ...)
	b = [ sum(yi.*(xi.^i)) for i =0:n]
	a = A\b


	f(x::Number) = a' * x.^(0:n)
	f(x::Vector{<:Number}) = [f(x::Number) for x=x]

	return f
end


function GaussQuadrature(n::Int64)
	b = [ (i+1)/((2*i+1)*(2*i+3))^0.5 for i=0:n-2]
	K = diagm(1=>b,-1=>b)
	E = eigen(K)

	wg = 2*E.vectors[1,:].^2
	xg = E.values
	return xg,wg
end
