function a = coefficients(p,x)
[c,t] = coeffs(p,x);
degrees = double(diff(log(t),x)*x);
n = max(degrees);
a = sym(zeros(1,n+1));
a(1+degrees)=c;