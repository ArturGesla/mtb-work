function p=prob(N,x)
p=factorial(N)/factorial((N-x)/2)/factorial((N+x)/2)/2^N;
end