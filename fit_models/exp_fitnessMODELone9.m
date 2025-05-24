function ee=exp_fitnessMODEL9(k,alpha, x, s, ppno,pphi,ifile) 
 lambda1=k(1); 
 lambda2=k(2); 
 lambda3=k(3); 
 lambda4=k(4); 
 A1=k(6); 
 A2=k(7); 
 A3=k(8); 
 A4=1-A1-A2-A3; 
 k1=k(8); 
 k2=k(9); 
 k3=k(10); 
 k4=k(11); 
 k5=k(12); 
 k9=k(13); 
 k10=k(14); 
 L1 = lambda1 + lambda2 + lambda3 + lambda4; 
 L2 = lambda1*lambda2 + lambda1*lambda3 + lambda1*lambda4 + lambda2*lambda3 + lambda2*lambda4 + lambda3*lambda4; 
 L3 = lambda1*lambda2*lambda3 + lambda1*lambda2*lambda4 + lambda1*lambda3*lambda4   +   lambda2*lambda3*lambda4; 
 L4 = lambda1 * lambda2 * lambda3 * lambda4; 
 S1 = A1*lambda1  + A2*lambda2   + A3*lambda3   + A4*lambda4; 
 S2 = A1*lambda1^2+ A2*lambda2^2 + A3*lambda3^2 + A4*lambda4^2; 
 S3 = A1*lambda1^3+ A2*lambda2^3 + A3*lambda3^3 + A4*lambda4^3; 
 LL1 = - k1 - k2 - k3 - k4 - k5 - k9 - k10; 
 LL2 = k1*k3 + k1*k4 + k1*k5 + k2*k4 + k2*k5 + k3*k5 + k1*k9 + k1*k10 + k2*k9 + k2*k10 + k3*k9 + k3*k10 + k4*k9 + k4*k10 + k5*k9 + k5*k10; 
 LL3 = - k1*k3*k5 - k1*k3*k9 - k1*k3*k10 - k1*k4*k9 - k1*k4*k10 - k1*k5*k9 - k2*k4*k9 - k1*k5*k10 - k2*k4*k10 - k2*k5*k9 - k2*k5*k10 - k3*k5*k10; 
 LL4 = k1*k3*k5*k10; 
 SS2 = -k5*k10; 
 SS3 = -k5*k10*(L1 + k1 + k2 + k3); 
 pt=k1*k3*k5 + k1*k3*k9 + k1*k3*k10 + k1*k4*k9 + k1*k4*k10 + k1*k5*k9 + k2*k4*k9 + k2*k4*k10 + k2*k5*k9; 
 p1= k2*(k4*k9 + k4*k10 + k5*k9)/pt; 
 p2= k1*(k4*k9 + k4*k10 + k5*k9)/pt; 
 p3= k1*k3*(k9 + k10)/pt; 
 p4= k1*k3*k5/pt; 
 SS1 = abs(A1*lambda1)+ abs(A2*lambda2) + abs(A3*lambda3) + abs(A4*lambda4); 
 N = length(x); 
 sN= sqrt(N); 
 fact1=sqrt(1-alpha); 
 fact2=sqrt(alpha); 
 if ifile==2 
     pp=pphi; 
 else 
     pp=ppno; 
 end 
 constr=[(L1-LL1)/L1;(L2-LL2)/L2;(L3-LL3)/L3;(L4-LL4)/L4;S1/SS1;(S2-SS2)/S2;(S3-SS3)/S3;(p1-pp(1))/pp(1);(p2-pp(2))/pp(2);(p3-pp(3))/pp(3);(p4-pp(4))/pp(4)]; 
 ps = A1*exp(lambda1*x) + A2*exp(lambda2*x)+A3*exp(lambda3*x)+A4*exp(lambda4*x); 
 if alpha == 1 
     ee =  [(ps-s)/sN;constr]/sqrt(12); 
 else 
     pss = ps./s;pss(pss<=0)=100; 
     psl = pl./sl;psl(psl<=0)=100; 
     ee = [(log(pss))/sN*fact1;(ps-s)/sN*fact2;constr]/sqrt(12); 
 end 
 