LL1 = - k1 - k2 - k3 - k4 - k5 - k6 - k8 - k9 - k10; 
 LL2 = k1*k3 + k1*k4 + k1*k5 + k2*k4 + k1*k6 + k2*k5 + k2*k6 + k3*k5 + k1*k8 + k3*k6 + k1*k9 + k2*k8 + k4*k6 + k1*k10 + k2*k9 + k3*k8 + k5*k6 + k2*k10 + k3*k9 + k4*k8 + k3*k10 + k4*k9 + k4*k10 + k5*k9 + k5*k10 + k6*k9 + k8*k9 + k9*k10; 
 LL3 = - k1*k3*k5 - k1*k3*k6 - k1*k4*k6 - k1*k3*k8 - k1*k5*k6 - k2*k4*k6 - k1*k3*k9 - k1*k4*k8 - k2*k5*k6 - k1*k3*k10 - k1*k4*k9 - k2*k4*k8 - k3*k5*k6 - k1*k4*k10 - k1*k5*k9 - k2*k4*k9 - k1*k5*k10 - k1*k6*k9 - k2*k4*k10 - k2*k5*k9 - k2*k5*k10 - k2*k6*k9 - k3*k5*k9 - k1*k8*k9 - k3*k5*k10 - k3*k6*k9 - k2*k8*k9 - k4*k6*k9 - k1*k9*k10 - k3*k8*k9 - k5*k6*k9 - k2*k9*k10 - k4*k8*k9 - k3*k9*k10 - k4*k9*k10 - k5*k9*k10; 
 LL4 = k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k5*k10 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k1*k5*k9*k10 + k2*k4*k9*k10 + k2*k5*k9*k10 + k3*k5*k9*k10; 
 LL5 = -k1*k3*k5*k9*k10; 
 SS1 = 0; 
 SS2 = -k5*k10; 
 SS3 = -k5*k10*(L1 + k1 + k2 + k3 + k9); 
 SS4 = -k10*(L1*(k5*k9 + k5*(k1 + k2 + k3)) - k5*(L2 - L1^2) + k1*k3*k5 + k5*k9*(k1 + k2 + k3)); 
 p1= (k2*k9*(k4*k6 + k5*k6 + k4*k8 + k4*k10))/(k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10) 
 p2= (k1*k9*(k4*k6 + k5*k6 + k4*k8 + k4*k10))/(k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10) 
 p3= (k1*k3*k9*(k6 + k8 + k10))/(k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10) 
 p4= (k1*k3*k5*k6)/(k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10) 
 p5= (k1*k3*k5*k9)/(k1*k3*k5*k6 + k1*k3*k5*k9 + k1*k3*k6*k9 + k1*k4*k6*k9 + k1*k3*k8*k9 + k1*k5*k6*k9 + k2*k4*k6*k9 + k1*k4*k8*k9 + k2*k5*k6*k9 + k1*k3*k9*k10 + k2*k4*k8*k9 + k1*k4*k9*k10 + k2*k4*k9*k10) 
 