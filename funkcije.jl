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




"""GaussQuadrature(n::Int64)

	x,w = GaussQuadrature(n)

	f(x)*w ~> Int(f,-1,1)
"""
function GaussQuadrature(n::Int64)
	b = [ (i+1)/((2*i+1)*(2*i+3))^0.5 for i=0:n-2]
	K = diagm(1=>b,-1=>b)
	E = eigen(K)

	wg = 2*E.vectors[1,:].^2
	xg = E.values
	return xg,wg
end

"""GaussIntegrate(f,a1,a2)

I = GaussIntegrate(f,a1,a2) = wi'*f(xi)
"""
function GaussIntegrate(f::Function,a1,a2;n=30)
	x,w = GaussQuadrature(n::Int64)

	I =w'*f.(x*(a2-a1)/2 .+(a2+a1)/2)*(a2-a1)/2
	I = round(I,digits = 14)
end
	

function TrigExpans(f::Function,a1::Number,a2::Number;n::Int64 = 60)
	x,w = GaussQuadrature(2*n)
	l1 = Lagrange([-1,1],[a1,a2])
	l2 = Lagrange([-1,1],[0,2*pi])
	l3 = Lagrange([a1,a2],[0,2*pi])


	s = [w'*(f.(l1.(x)).*sin.(i*l2(x))) for i = 1:n ]
	c = [w'*(f.(l1.(x)).*cos.(i*l2.(x))) for i = 1:n ] 
	a = GaussIntegrate(f,a1,a2;n=n)/(a2-a1)
	
	g(h) = a + s'*sin.((1:n)*l3(h))+c'*cos.((1:n)*l3(h))
	return g
end

	
