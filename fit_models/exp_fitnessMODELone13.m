function ee=exp_fitnessMODEL13(k,alpha, x, s, ppno,pphi,ifile) 
 lambda1=k(1); 
 lambda2=k(2); 
 lambda3=k(3); 
 lambda4=k(4); 
 lambda5=k(5); 
 A1=k(6); 
 A2=k(7); 
 A3=k(8); 
 A4=k(9); 
 A5=1-A1-A2-A3-A4; 
 k1=k(10); 
 k2=k(11); 
 k3=k(12); 
 k4=k(13); 
 k5=k(14); 
 k6=k(15); 
 k7=k(16); 
 k8=k(17); 
 k9=k(18); 
 k10=k(19); 
 L1 = lambda1 + lambda2 + lambda3 + lambda4 + lambda5; 
 L2 = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5; 
 L3 = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5; 
 L4 = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5; 
 L5 = lambda1*lambda2*lambda3*lambda4*lambda5; 
 S1 = A1*lambda1  + A2*lambda2   + A3*lambda3   + A4*lambda4   + A5*lambda5; 
 S2 = A1*lambda1^2+ A2*lambda2^2 + A3*lambda3^2 + A4*lambda4^2 + A5*lambda5^2; 
 S3 = A1*lambda1^3+ A2*lambda2^3 + A3*lambda3^3 + A4*lambda4^3 + A5*lambda5^3; 
 S4 = A1*lambda1^4+ A2*lambda2^4 + A3*lambda3^4 + A4*lambda4^4 + A5*lambda5^4; 
 LL1 = - k1 - k2 - k3 - k4 - k5 - k6 - k7 - k8 - k9 - k10; 
 LL2 = k1*k3 + k1*k4 + k1*k5 + k2*k4 + k1*k6 + k2*k5 + k1*k7 + k2*k6 + k3*k5 + k1*k8 + k2*k7 + k3*k6 + k1*k9 + k2*k8 + k3*k7 + k4*k6 + k1*k10 + k2*k9 + k3*k8 + k4*k7 + k5*k6 + k2*k10 + k3*k9 + k4*k8 + k5*k7 + k3*k10 + k4*k9 + k4*k10 + k5*k9 + k5*k10 + k6*k9 + k7*k8 + k7*k10 + k8*k9 + k9*k10; 
 LL3 = - k1*k3*k5 - k1*k3*k6 - k1*k3*k7 - k1*k4*k6 - k1*k3*k8 - k1*k4*k7 - k1*k5*k6 - k2*k4*k6 - k1*k3*k9 - k1*k4*k8 - k1*k5*k7 - k2*k4*k7 - k2*k5*k6 - k1*k3*k10 - k1*k4*k9 - k2*k4*k8 - k2*k5*k7 - k3*k5*k6 - k1*k4*k10 - k1*k5*k9 - k2*k4*k9 - k3*k5*k7 - k1*k5*k10 - k1*k6*k9 - k1*k7*k8 - k2*k4*k10 - k2*k5*k9 - k2*k5*k10 - k2*k6*k9 - k2*k7*k8 - k3*k5*k9 - k1*k7*k10 - k1*k8*k9 - k3*k5*k10 - k3*k6*k9 - k3*k7*k8 - k2*k7*k10 - k2*k8*k9 - k4*k6*k9 - k4*k7*k8 - k1*k9*k10 - k3*k7*k10 - k3*k8*k9 - k2*k9*k10 - k4*k7*k10 - k4*k8*k9 - k3*k9*k10 - k5*k7*k10 - k4*k9*k10 - k5*k9*k10; 
 LL4 = k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k5*k9 + k1*k3*k5*k10 + k1*k3*k6*k9 + k1*k3*k7*k8 + k1*k4*k6*k9 + k1*k4*k7*k8 + k1*k3*k7*k10 + k1*k3*k8*k9 + k2*k4*k6*k9 + k2*k4*k7*k8 + k1*k4*k7*k10 + k1*k4*k8*k9 + k1*k3*k9*k10 + k1*k5*k7*k10 + k2*k4*k7*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k5*k7*k10 + k1*k5*k9*k10 + k2*k4*k9*k10 + k3*k5*k7*k10 + k2*k5*k9*k10 + k3*k5*k9*k10; 
 LL5 = -k1*k3*k5*k10*(k7 + k9); 
 SS2 = -k5*k10; 
 SS3 = -k5*k10*(L1 + k1 + k2 + k3 + k7 + k9); 
 SS4 = -k10*(L1*k5*(k1 + k2 + k3 + k7 + k9) - k5*(L2 - L1^2) + k5*(k7 + k9)*(k1 + k2 + k3) + k1*k3*k5); 
 pt=k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k3*k7*k8 + k1*k4*k6*k9 + k1*k4*k7*k8 + k1*k3*k7*k10 + k1*k3*k8*k9 + k2*k4*k6*k9 + k2*k4*k7*k8 + k1*k4*k7*k10 + k1*k4*k8*k9 + k1*k3*k9*k10 + k2*k4*k7*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10; 
 p1= k2*k4*(k6*k9 + k7*k8 + k7*k10 + k8*k9 + k9*k10)/pt; 
 p2= k1*k4*(k6*k9 + k7*k8 + k7*k10 + k8*k9 + k9*k10)/pt; 
 p3= k1*k3*(k6*k9 + k7*k8 + k7*k10 + k8*k9 + k9*k10)/pt; 
 p4= k1*k3*k5*k6/pt; 
 p5= k1*k3*k5*(k7 + k9)/pt; 
 SS1 = abs(A1*lambda1)+ abs(A2*lambda2) + abs(A3*lambda3) + abs(A4*lambda4) + abs(A5*lambda5); 
 N = length(x); 
 sN= sqrt(N); 
 fact1=sqrt(1-alpha); 
 fact2=sqrt(alpha); 
 if ifile==2 
     pp=pphi; 
 else 
     pp=ppno; 
 end 
 constr=[(L1-LL1)/L1;(L2-LL2)/L2;(L3-LL3)/L3;(L4-LL4)/L4;(L5-LL5)/L5;S1/SS1;(S2-SS2)/S2;(S3-SS3)/S3;(S4-SS4)/S4;(p1-pp(1))/pp(1);(p2-pp(2))/pp(2);(p3-pp(3))/pp(3);(p4-pp(4))/pp(4);(p5-pp(5))/pp(5)]; 
 ps = A1*exp(lambda1*x) + A2*exp(lambda2*x)+A3*exp(lambda3*x)+A4*exp(lambda4*x)+A5*exp(lambda5*x); 
 if alpha == 1 
     ee =  [(ps-s)/sN;constr]/sqrt(15); 
 else 
     pss = ps./s;pss(pss<=0)=100; 
     ee = [(log(pss))/sN*fact1;(ps-s)/sN*fact2;constr]/sqrt(15); 
 end 
 