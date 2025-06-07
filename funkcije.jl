function Lagrange(xi,yi)
	n = length(xi)
	f(x) = sum([ yi[i]/prod(xi[setdiff(1:n,i)] .- xi[i] )*prod(xi[setdiff(1:n,i)].-x) for i = 1:n])


	return f
end
