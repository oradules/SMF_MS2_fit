LL1 = - k1 - k2 - k3 - k4 - k5 - k6 - k7 - k8; 
 LL2 = k1*k3 + k1*k4 + k1*k5 + k2*k4 + k1*k6 + k2*k5 + k1*k7 + k2*k6 + k3*k5 + k1*k8 + k2*k7 + k3*k6 + k2*k8 + k3*k7 + k4*k6 + k3*k8 + k4*k7 + k5*k6 + k4*k8 + k5*k7 + k5*k8 + k7*k8; 
 LL3 = - k1*k3*k5 - k1*k3*k6 - k1*k3*k7 - k1*k4*k6 - k1*k3*k8 - k1*k4*k7 - k1*k5*k6 - k2*k4*k6 - k1*k4*k8 - k1*k5*k7 - k2*k4*k7 - k2*k5*k6 - k1*k5*k8 - k2*k4*k8 - k2*k5*k7 - k3*k5*k6 - k2*k5*k8 - k3*k5*k7 - k1*k7*k8 - k3*k5*k8 - k2*k7*k8 - k3*k7*k8 - k4*k7*k8 - k5*k7*k8; 
 LL4 = k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k5*k8 + k1*k3*k7*k8 + k1*k4*k7*k8 + k1*k5*k7*k8 + k2*k4*k7*k8 + k2*k5*k7*k8 + k3*k5*k7*k8; 
 LL5 = -k1*k3*k5*k7*k8; 
 SS1 = 0; 
 SS2 = -k5*k8; 
 SS3 = -k5*k8*(L1 + k1 + k2 + k3 + k7); 
 SS4 = -k5*k8*(L1*k1 - L2 + L1*k2 + L1*k3 + L1*k7 + k1*k3 + k1*k7 + k2*k7 + k3*k7 + L1^2); 
 p1= (k2*k4*k7*k8)/(k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k7*k8 + k1*k4*k7*k8 + k2*k4*k7*k8) 
 p2= (k1*k4*k7*k8)/(k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k7*k8 + k1*k4*k7*k8 + k2*k4*k7*k8) 
 p3= (k1*k3*k7*k8)/(k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k7*k8 + k1*k4*k7*k8 + k2*k4*k7*k8) 
 p4= (k1*k3*k5*k6)/(k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k7*k8 + k1*k4*k7*k8 + k2*k4*k7*k8) 
 p5= (k1*k3*k5*k7)/(k1*k3*k5*k6 + k1*k3*k5*k7 + k1*k3*k7*k8 + k1*k4*k7*k8 + k2*k4*k7*k8) 
 