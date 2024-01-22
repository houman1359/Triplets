function Om=Omega_Student(n,r,d)



Om=( log( ((beta(n/2,1/2)^d)*gamma(d/2))/((pi^(d/2))*beta(n/2,d/2)) ) -(n*(d-1))/2*psi(n/2) +(d*(n+1))/2*psi((n+1)/2) -(n+d)/2*psi((n+d)/2) )/log(2)-0.5*log2(1-r^2);
