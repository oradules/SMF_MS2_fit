%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, August
%%%% 2024
%%%%%%% compute the constraints symbolically for 5
%%%%%%% states models, generates code with objective functions for each model and ML problem 
%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MODELS: 1-8, 13,14,17,18
%%%%%%% The objective functions are named exp_fitnessMODELXmergeonekini_mergep4p5ponkinishort.m
%%%%%%%  data consists of: survival functions for no and hi Tat, occupation probabilities, ponkini 
%%%%%%%  p4 and p5 are merged and imposed as a sum rather than individually









%%%% remember to add generation of tikz code for transition diagrams  


clear all



syms k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 lambda lambda1 lambda2 lambda3 lambda4 lambda5 C1 C2 C3 C4 C5 u1 u2 u3 u4 u5 A1 A2 A3 A4 A5 S1 S2 S3 S4 S5 L1 L2 L3 L4 p1 p2 p3 p4
syms L1 L2 L3 L4 L5 S1 S2 S3 S4 h1 h2 h3 
syms p1 p2 p3 p4 p5
syms l1 l2 l3   

%%%% generate matrix Q

%Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),0,0;0,0,0,-k7,k6;0,0,k5,k7,-k8-k6];

%Q = [-k1, k2, 0, 0, 0; ...
%      k1,-(k2+k3),k4,k10,0;...
%      0,k3,-(k4+k5),0,k9;...
%      0,0,0,-k7-k10,k6;...
%      0,0,k5,k7,-k8-k6-k9];



for iii=1:12
iii

switch iii
%%%%% Model definition
    case 1
%%%%% MODEL 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),0,0;0,0,0,-k7,k6;0,0,k5,k7,-k8-k6];
kini=k8;nini=8;
npars =8; nparinv=5;
idel =[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
%%%%% MODEL 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,k9,0;0,k3,-(k4+k5),0,k8;0,0,0,-(k7+k9),k6;0,0,k5,k7,-(k8+k6+k10)];
kini=k10;nini=10;
npars=10; nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
%%%%% MODEL 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, k9, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),0,k8;0,0,0,-(k7+k9),k6;0,0,k5,k7,-(k8+k6+k10)];
kini=k10;nini=10;
npars=10; nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
%%%%% MODEL 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3; %%% 9 pars
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,k9,0;0,k3,-(k4+k5),0,k8;0,0,0,-k9,k6;0,0,k5,0,-(k8+k6+k10)];
kini=k10;nini=10;
npars =10; nparinv=5;
idel = [7]; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
%%%%% MODEL 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3; %%%% 9 pars
Q = [-k1, k2, 0, k9, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),0,k8;0,0,0,-k9,k6;0,0,k5,0,-(k8+k6+k10)];
kini=k10;nini=10;
npars = 10; nparinv=5;
idel = [7]; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
%%%%% MODEL 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;  %%% 8 pars
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,k9,0;0,k3,-(k4+k5),0,0;0,0,0,-k9,k6;0,0,k5,0,-(k6+k10)];
kini=k10;nini=10;
npars = 10; nparinv=5;
idel = [7,8]; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
%%%%% MODEL 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3; %%%% 8 pars
Q = [-k1, k2, 0, k9, 0; k1,-(k2+k3),k4,0,0; 0,k3,-(k4+k5),0,0; 0,0,0,-k9,k6;0,0,k5,0,-(k6+k10)];
kini=k10;nini=10;
npars = 10; nparinv=5;
idel = [7,8]; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
%%%%% MODEL 8 no common par %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;  %%% 8 pars
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0; 0,k3,-(k4+k5),k9,0; 0,0,0,-k9,k6;0,0,k5,0,-(k6+k10)];
kini=k10;nini=10;
npars = 10; nparinv=5;
idel = [7,8]; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9
%%%%% MODEL 13 no common par %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),k9,k8;0,0,0,-(k7+k9),k6;0,0,k5,k7,-(k8+k6+k10)];
kini=k10;nini=10;
npars=10; nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 10
%%%%% MODEL 14 no common par %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, k12, 0; k1,-(k2+k3),k4,k11,0; 0,k3,-(k4+k5),k9,k8; 0,0,0,-(k7+k9+k11+k12),k6; 0,0,k5,k7,-(k8+k6+k10)];
kini=k10;nini=10;
npars=12;nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 11
%%%%% MODEL 17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0;0,k3,-(k4+k5),0,k8;0,0,0,-k7,k6;0,0,k5,k7,-(k8+k6+k10)];
kini=k10;nini=10;
npars=10; nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 12
%%%%% MODEL 18 no common par %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; s=3;  %%% 8 pars
Q = [-k1, k2, 0, 0, 0; k1,-(k2+k3),k4,0,0; 0,k3,-(k4+k5),k9,k8; 0,0,0,-k9,k6;0,0,k5,0,-(k6+k8+k10)];
kini=k10;nini=10;
npars = 10; nparinv=5;
idel = []; %%%%% parameters to be set to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%%%%% Qred
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qp = Q;
Qp(s,N)=Qp(s,N)+kini;



%%%% compute occupation probabilities
res5=solve ([transpose((Qp*[p1;p2;p3;p4;p5]==[0;0;0;0;0])),...
    p1+p2+p3+p4+p5==1],p1,p2,p3,p4,p5);



%%%% system for eigenvectors start from state s
u = [];
for i=1:N
    u = [u;sym(['u',num2str(i)])];
end
uu = transpose(u);
u(s)=1;
uu(s)=[];

M=lambda*eye(N)-Q;
%%%% characteristic polynomial
plambda = collect(det(M),'lambda');
d = coefficients(plambda,lambda); %%%% [a0,a1,...,aN-1] coefficients of characteristic polynomial in order of increasing degrees: 0, 1, 2, 3, 4


%%%%%% compute eigenvectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syst=M*u;
Syst = [syst(1) == 0, syst(2) == 0, syst(3) == 0, syst(4) == 0, syst(5)==0];
Syst(s)=[];
res = solve (Syst, uu(1),uu(2),uu(3),uu(4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% polynomials Dsi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sym(zeros(N,1));
for i=1:N
    MM=M;MM(3,:)=[];MM(:,i)=[];
    D(i)=collect(det(MM),lambda);   
end
c=coefficients(D(N),lambda); %%%% DN =[c0,...,cN-2]


%%% test eigenvectors
for i=1:N-1
   fn=char(uu(i));
   ii=str2num(strrep(fn,'u',''));
   uii=getfield(res,fn);
   simplify(uii - D(ii)/D(s)*(-1)^(abs(s-ii)))
end



%%% compute Si
SS=sym(zeros(N-1,1));
sig = (-1)^(abs(s-N)+1);
if s < N
   SS(1)=0;
   SS(2)= sig*kini*c(4);
   SS(3)= sig*kini*(c(4)*L1+c(3));
   SS(4)= sig*kini*(c(4)*(L1^2-L2)+c(3)*L1+c(2));
else
   SS(1) = -kini;
   SS(2) = -kini*(c(5)*L1 + c(4));
   SS(3) = -kini*(c(5)*(L1^2-L2) + c(4)*L1 + c(3));
   SS(4) = -kini*(c(5)*(L1^3-2*L1*L2+L3)+c(4)*(L1^2-L2) + c(3)*L1 + c(2));
end
    







%%% test C
%%%  
%res2 =solve( ...
%C1*subs(res.u1,lambda,lambda1) + C2*subs(res.u1,lambda,lambda2) + C3*subs(res.u1,lambda,lambda3) + C4*subs(res.u1,lambda,lambda4) + C5*subs(res.u1,lambda,lambda5) == 0, ... 
%C1*subs(res.u2,lambda,lambda1) + C2*subs(res.u2,lambda,lambda2) + C3*subs(res.u2,lambda,lambda3) + C4*subs(res.u2,lambda,lambda4) + C5*subs(res.u2,lambda,lambda5) == 0, ... 
%C1 + C2 + C3 + C4 + C5 == 1,...
%C1*subs(res.u4,lambda,lambda1) + C2*subs(res.u4,lambda,lambda2) + C3*subs(res.u4,lambda,lambda3) + C4*subs(res.u4,lambda,lambda4) + C5*subs(res.u4,lambda,lambda5) == 0, ...
%C1*subs(res.u5,lambda,lambda1) + C2*subs(res.u5,lambda,lambda2) + C3*subs(res.u5,lambda,lambda3) + C4*subs(res.u5,lambda,lambda4) + C5*subs(res.u5,lambda,lambda5) == 0, ...    
%C1, C2, C3, C4, C5);

%simplify(res2.C1 - subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5)))
%simplify(res2.C2 - subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5)))
%simplify(res2.C3 - subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5)))
%simplify(res2.C4 - subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5)))
%simplify(res2.C5 - subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4)))

%simplify(subs(subs(subs(subs(subs(C1*subs(res.u1,lambda,lambda1) + C2*subs(res.u1,lambda,lambda2) + C3*subs(res.u1,lambda,lambda3) + C4*subs(res.u1,lambda,lambda4) + C5*subs(res.u1,lambda,lambda5),...
%C1, subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))),C2,...
%subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))),C3,...
%subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))),C4,...
%subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))),C5,...
%subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))))


%simplify(subs(subs(subs(subs(subs(C1*subs(res.u4,lambda,lambda1) + C2*subs(res.u4,lambda,lambda2) + C3*subs(res.u4,lambda,lambda3) + C4*subs(res.u4,lambda,lambda4) + C5*subs(res.u4,lambda,lambda5),...
%C1, subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))),C2,...
%subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))),C3,...
%subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))),C4,...
%subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))),C5,...
%subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))))


%simplify(subs(subs(subs(subs(subs(C1*subs(res.u2,lambda,lambda1) + C2*subs(res.u2,lambda,lambda2) + C3*subs(res.u2,lambda,lambda3) + C4*subs(res.u2,lambda,lambda4) + C5*subs(res.u2,lambda,lambda5),...
%C1, subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))),C2,...
%subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))),C3,...
%subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))),C4,...
%subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))),C5,...
%subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))))

%simplify(subs(subs(subs(subs(subs(C1*subs(res.u5,lambda,lambda1) + C2*subs(res.u5,lambda,lambda2) + C3*subs(res.u5,lambda,lambda3) + C4*subs(res.u5,lambda,lambda4) + C5*subs(res.u5,lambda,lambda5),...
%C1, subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))),C2,...
%subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))),C3,...
%subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))),C4,...
%subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))),C5,...
%subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))))

%simplify(subs(subs(subs(subs(subs(C1+C2+C3+C4+C5-1,...
%C1, subs(D3,lambda,lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))),C2,...
%subs(D3,lambda,lambda2)/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))),C3,...
%subs(D3,lambda,lambda3)/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))),C4,...
%subs(D3,lambda,lambda4)/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))),C5,...
%subs(D3,lambda,lambda5)/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))))






%collect(simplify((Q- lambda*eye(N))*[D1;-D2;D3;-D4;D5]/D3),lambda)

%%%  
%res2 =solve( ...
%C1*subs(res.u1,lambda,lambda1) + C2*subs(res.u1,lambda,lambda2) + C3*subs(res.u1,lambda,lambda3) + C4*subs(res.u1,lambda,lambda4) + C5*subs(res.u1,lambda,lambda5) == 0, ... 
%C1*subs(res.u2,lambda,lambda1) + C2*subs(res.u2,lambda,lambda2) + C3*subs(res.u2,lambda,lambda3) + C4*subs(res.u2,lambda,lambda4) + C5*subs(res.u2,lambda,lambda5) == 0, ... 
%C1 + C2 + C3 + C4 + C5 == 1,...
%C1*subs(res.u4,lambda,lambda1) + C2*subs(res.u4,lambda,lambda2) + C3*subs(res.u4,lambda,lambda3) + C4*subs(res.u4,lambda,lambda4) + C5*subs(res.u4,lambda,lambda5) == 0, ...
%C1*subs(res.u5,lambda,lambda1) + C2*subs(res.u5,lambda,lambda2) + C3*subs(res.u5,lambda,lambda3) + C4*subs(res.u5,lambda,lambda4) + C5*subs(res.u5,lambda,lambda5) == 0, ...    
%C1, C2, C3, C4, C5);
%r=rank([subs(res.u1,lambda,lambda1),subs(res.u1,lambda,lambda2),subs(res.u1,lambda,lambda3),subs(res.u1,lambda,lambda4),subs(res.u1,lambda,lambda5);...
%    subs(res.u2,lambda,lambda1),subs(res.u2,lambda,lambda2),subs(res.u2,lambda,lambda3),subs(res.u2,lambda,lambda4),subs(res.u2,lambda,lambda5);...
%    1,1,1,1,1;...
%subs(res.u4,lambda,lambda1),subs(res.u4,lambda,lambda2),subs(res.u4,lambda,lambda3),subs(res.u4,lambda,lambda4),subs(res.u4,lambda,lambda5);...
%subs(res.u5,lambda,lambda1),subs(res.u5,lambda,lambda2),subs(res.u5,lambda,lambda3),subs(res.u5,lambda,lambda4),subs(res.u5,lambda,lambda5)]);  



%d=coefficients(det(M),lambda);


%res=solve(-L1 == d(5),L2 == d(4),...
%           S2 == -k8*c(4),S3==-k8*(c(4)*L1+c(3)),...
%           S4 == -k8*(c(4)*(L1^2-L2)+c(3)*L1+c(2)),k3,k4,k5,k6,k7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iii < 9
fname = ['proba5MODEL',num2str(iii),'.m']
else
    if iii < 11
fname = ['proba5MODEL',num2str(iii+4),'.m']
    else
fname = ['proba5MODEL',num2str(iii+6),'.m']        
    end
end
fid=fopen(fname,'w');
if iii < 9
	str=['function [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL',num2str(iii),'(k);'];
else
    if iii < 11
    str=['function [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL',num2str(iii+4),'(k);'];
    else
    str=['function [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL',num2str(iii+6),'(k);'];    
    end
end
fprintf(fid,'%s \n ',str);

j=19;%%% index of k
if nparinv > 0
for i=1:nparinv
   if ~ismember(i,idel)   
        str=['k',num2str(i),'=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end
end
end
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=['k',num2str(i),'no=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end 
end
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=['k',num2str(i),'hi=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end 
end


[N,D] = numden(res5.p1);
str=['pt=',char(D),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p1no= (',char(simplify(res5.p1*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p2no= (',char(simplify(res5.p2*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p3no= (',char(simplify(res5.p3*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p4no= (',char(simplify(res5.p4*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p5no= (',char(simplify(res5.p5*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));

[N,D] = numden(res5.p1);
str=['pt=',char(D),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p1hi= (',char(simplify(res5.p1*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p2hi= (',char(simplify(res5.p2*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p3hi= (',char(simplify(res5.p3*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p4hi= (',char(simplify(res5.p4*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p5hi= (',char(simplify(res5.p5*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));



fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iii < 9
    filename =['matrices_MODEL',num2str(iii),'.txt']  ;
else
    if iii < 11
        filename =['matrices_MODEL',num2str(iii+6),'.txt']  ;
    else
       filename =['matrices_MODEL',num2str(iii+6),'.txt']  ;     
    end
end
fid=fopen(filename,'w');

str= 'Q=';
fprintf(fid,'%s \n ',str);
for i=1:5
str= char(Q(i,:));
fprintf(fid,'%s \n ',str);
end
str= 'Qp=';
fprintf(fid,'%s \n ',str);
for i=1:5
str= char(Qp(i,:));
fprintf(fid,'%s \n ',str);
end


fclose(fid);

%%%%%%%%%%%%%changed!
if iii < 9
fname = ['exp_fitnessMODEL',num2str(iii),'mergeonekini_mergep4p5ponkinishort.m']
else
    if iii < 11
fname = ['exp_fitnessMODEL',num2str(iii+4),'mergeonekini_mergep4p5ponkinishort.m']
    else
fname = ['exp_fitnessMODEL',num2str(iii+6),'mergeonekini_mergep4p5ponkinishort.m']        
    end
end
fid=fopen(fname,'w');

%%%%%%%%% changed!
if iii < 9
	str=['function ee=exp_fitnessMODEL',num2str(iii),'mergeonekini_mergep4p5ponkinishort(k, alpha, xno, sno, xhi, shi,  ppno,pphi, ponkinino, ponkinihi);'];
else
    if iii < 11
        str=['function ee=exp_fitnessMODEL',num2str(iii+4),'mergeonekini_mergep4p5ponkinishort(k, alpha, xno, sno, xhi, shi,  ppno,pphi, ponkinino, ponkinihi);'];
    else
        str=['function ee=exp_fitnessMODEL',num2str(iii+6),'mergeonekini_mergep4p5ponkinishort(k, alpha, xno, sno, xhi, shi,  ppno,pphi, ponkinino, ponkinihi);'];   
    end
end

fprintf(fid,'%s \n ',str);





for i=1:5
    str=['lambda',num2str(i),'=k(',num2str(i),');'];
    fprintf(fid,'%s \n ',str);
end
for i=1:4
    str=['A',num2str(i),'=k(',num2str(i+5),');'];
    fprintf(fid,'%s \n ',str);
end
str='A5=1-A1-A2-A3-A4;';
fprintf(fid,'%s \n ',str);

str='L1no = lambda1 + lambda2 + lambda3 + lambda4 + lambda5;';
fprintf(fid,'%s \n ',str);
str='L2no = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L3no = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L4no = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L5no = lambda1*lambda2*lambda3*lambda4*lambda5;';
fprintf(fid,'%s \n ',str);

%%%%%%%%% changed!
str='pno = A1*exp(lambda1*xno) + A2*exp(lambda2*xno)+A3*exp(lambda3*xno)+A4*exp(lambda4*xno)+A5*exp(lambda5*xno);';
fprintf(fid,'%s \n ',str);


str='S1no = A1*lambda1  + A2*lambda2   + A3*lambda3   + A4*lambda4   + A5*lambda5;'; 
fprintf(fid,'%s \n ',str);
str='S2no = A1*lambda1^2+ A2*lambda2^2 + A3*lambda3^2 + A4*lambda4^2 + A5*lambda5^2;';
fprintf(fid,'%s \n ',str);
str='S3no = A1*lambda1^3+ A2*lambda2^3 + A3*lambda3^3 + A4*lambda4^3 + A5*lambda5^3;';
fprintf(fid,'%s \n ',str);
str='S4no = A1*lambda1^4+ A2*lambda2^4 + A3*lambda3^4 + A4*lambda4^4 + A5*lambda5^4;';
fprintf(fid,'%s \n ',str);
str='SS1no = abs(A1*lambda1)+ abs(A2*lambda2) + abs(A3*lambda3) + abs(A4*lambda4) + abs(A5*lambda5);';
fprintf(fid,'%s \n ',str);

for i=1:5
    str=['lambda',num2str(i),'=k(',num2str(i+9),');'];
    fprintf(fid,'%s \n ',str);
end
for i=1:4
    str=['A',num2str(i),'=k(',num2str(i+14),');'];
    fprintf(fid,'%s \n ',str);
end
str='A5=1-A1-A2-A3-A4;';
fprintf(fid,'%s \n ',str);


str='L1hi = lambda1 + lambda2 + lambda3 + lambda4 + lambda5;';
fprintf(fid,'%s \n ',str);
str='L2hi = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L3hi = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L4hi = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5;';
fprintf(fid,'%s \n ',str);
str='L5hi = lambda1*lambda2*lambda3*lambda4*lambda5;';
fprintf(fid,'%s \n ',str);

%%%% changed !
str='phi = A1*exp(lambda1*xhi) + A2*exp(lambda2*xhi)+A3*exp(lambda3*xhi)+A4*exp(lambda4*xhi)+A5*exp(lambda5*xhi);';
fprintf(fid,'%s \n ',str);


str='S1hi = A1*lambda1  + A2*lambda2   + A3*lambda3   + A4*lambda4   + A5*lambda5;'; 
fprintf(fid,'%s \n ',str);
str='S2hi = A1*lambda1^2+ A2*lambda2^2 + A3*lambda3^2 + A4*lambda4^2 + A5*lambda5^2;';
fprintf(fid,'%s \n ',str);
str='S3hi = A1*lambda1^3+ A2*lambda2^3 + A3*lambda3^3 + A4*lambda4^3 + A5*lambda5^3;';
fprintf(fid,'%s \n ',str);
str='S4hi = A1*lambda1^4+ A2*lambda2^4 + A3*lambda3^4 + A4*lambda4^4 + A5*lambda5^4;';
fprintf(fid,'%s \n ',str);
str='SS1hi = abs(A1*lambda1)+ abs(A2*lambda2) + abs(A3*lambda3) + abs(A4*lambda4) + abs(A5*lambda5);';
fprintf(fid,'%s \n ',str);

j=19;%%% index of k
if nparinv > 0
for i=1:nparinv
   if ~ismember(i,idel)   
        str=['k',num2str(i),'=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end
end
end
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=['k',num2str(i),'no=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end 
end
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=['k',num2str(i),'hi=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end 
end
%%% add 2 more parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%str = 'ponkinino = kinino;';
%fprintf(fid,'%s \n ',str);
%str = 'ponkinihi = kinihi;';
%fprintf(fid,'%s \n ',str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%str = ['kinino = kinino / k(',num2str(j),');'];
%fprintf(fid,'%s \n ',str);
%j=j+1;
%str = ['kinihi = kinihi / k(',num2str(j),');'];
%fprintf(fid,'%s \n ',str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


str=['LL1 = ', char(simplify(-d(5))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['LL2 = ',char(simplify(d(4))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['LL3 = ',char(simplify(-d(3))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['LL4 = ',char(simplify(d(2))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['LL5 = ',char(simplify(-d(1))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['SS2no = ',char(simplify(SS(2))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['SS3no = ',char(simplify(SS(3))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['SS4no = ',char(simplify(SS(4))),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
[N,D] = numden(res5.p1);
str=['pt=',char(D),';'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p1no= (',char(simplify(res5.p1*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p2no= (',char(simplify(res5.p2*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p3no= (',char(simplify(res5.p3*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p4no= (',char(simplify(res5.p4*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));
str=['p5no= (',char(simplify(res5.p5*D)),')/pt;'];
fprintf(fid,'%s \n ',strrepno(str,idel,npars,nparinv));


str=['LL1 = ', char(simplify(-d(5))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['LL2 = ',char(simplify(d(4))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['LL3 = ',char(simplify(-d(3))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['LL4 = ',char(simplify(d(2))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['LL5 = ',char(simplify(-d(1))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['SS2hi = ',char(simplify(SS(2))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['SS3hi = ',char(simplify(SS(3))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['SS4hi = ',char(simplify(SS(4))),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
[N,D] = numden(res5.p1);
str=['pt=',char(D),';'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p1hi= (',char(simplify(res5.p1*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p2hi= (',char(simplify(res5.p2*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p3hi= (',char(simplify(res5.p3*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p4hi= (',char(simplify(res5.p4*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));
str=['p5hi= (',char(simplify(res5.p5*D)),')/pt;'];
fprintf(fid,'%s \n ',strrephi(str,idel,npars,nparinv));

    str=['ponkininopred = p5no*k',num2str(nini),'no;'];
    fprintf(fid,'%s \n ',str);
    str=['ponkinihipred = p5hi*k',num2str(nini),'hi;'];
    fprintf(fid,'%s \n ',str);
%%%%%% fit model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% changed!
str='NS = length(xno);';fprintf(fid,'%s \n ',str);
str='NSno= sqrt(NS);';fprintf(fid,'%s \n ',str); 
str='fact1no=sqrt(1-alpha);';fprintf(fid,'%s \n ',str);
str='fact2no=sqrt(alpha);';fprintf(fid,'%s \n ',str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str='NS = length(xhi);';fprintf(fid,'%s \n ',str);
str='NShi= sqrt(NS);';fprintf(fid,'%s \n ',str); 
str='fact1hi=sqrt(1-alpha);';fprintf(fid,'%s \n ',str);
str='fact2hi=sqrt(alpha);';fprintf(fid,'%s \n ',str);

%%%% changed!
%str='constr=[(L1no-LL1no)/L1no;(L2no-LL2no)/L2no;(L3no-LL3no)/L3no;(L4no-LL4no)/L4no;(L5no-LL5no)/L5no;S1no/SS1no;(S2no-SS2no)/S2no;(S3no-SS3no)/S3no;(S4no-SS4no)/S4no;(p1no-ppno(1))/ppno(1);(p2no-ppno(2))/ppno(2);(p3no-ppno(3))/ppno(3);(p4no-ppno(4))/ppno(4);(p5no-ppno(5))/ppno(5);(L1hi-LL1hi)/L1hi;(L2hi-LL2hi)/L2hi;(L3hi-LL3hi)/L3hi;(L4hi-LL4hi)/L4hi;(L5hi-LL5hi)/L5hi;S1hi/SS1hi;(S2hi-SS2hi)/S2hi;(S3hi-SS3hi)/S3hi;(S4hi-SS4hi)/S4hi;(p1hi-pphi(1))/pphi(1);(p2hi-pphi(2))/pphi(2);(p3hi-pphi(3))/pphi(3);(p4hi-pphi(4))/pphi(4);(p5hi-pphi(5))/pphi(5)];';
str='constr=[(L1no-LL1no)/L1no;(L2no-LL2no)/L2no;(L3no-LL3no)/L3no;(L4no-LL4no)/L4no;(L5no-LL5no)/L5no;S1no/SS1no;(S2no-SS2no)/S2no;(S3no-SS3no)/S3no;(S4no-SS4no)/S4no;(p1no-ppno(1))/ppno(1);(p2no-ppno(2))/ppno(2);(p3no-ppno(3))/ppno(3);(p4no+p5no-ppno(4)-ppno(5))/(ppno(4)+ppno(5));(ponkininopred-ponkinino)/ponkinino;(L1hi-LL1hi)/L1hi;(L2hi-LL2hi)/L2hi;(L3hi-LL3hi)/L3hi;(L4hi-LL4hi)/L4hi;(L5hi-LL5hi)/L5hi;S1hi/SS1hi;(S2hi-SS2hi)/S2hi;(S3hi-SS3hi)/S3hi;(S4hi-SS4hi)/S4hi;(p1hi-pphi(1))/pphi(1);(p2hi-pphi(2))/pphi(2);(p3hi-pphi(3))/pphi(3);(p4hi+p5hi-pphi(4)-pphi(5))/(pphi(4)+pphi(5));(ponkinihipred-ponkinihi)/ponkinihi];';

%%% 28 elements
fprintf(fid,'%s \n ',str);

%%%% add the two kini constraints
%str=['constrnew=[(k',num2str(nini),'no-kinino)/kinino;(k',num2str(nini),'hi-kinihi)/kinihi];'];
%%% 2 elements
%fprintf(fid,'%s \n ',str);


str='if alpha == 1'; fprintf(fid,'%s \n ',str);
%str='    ee =  [(pno-sno)/NSno;(phi-shi)/NShi;constr;constrnew]/sqrt(32);';fprintf(fid,'%s \n ',str);
str='    ee =  [(pno-sno)/NSno;(phi-shi)/NShi;constr]/sqrt(30);';fprintf(fid,'%s \n ',str);
str='else';fprintf(fid,'%s \n ',str);
str='    psno = pno./sno;psno(psno<=0)=100;';fprintf(fid,'%s \n ',str);
str='    pshi = phi./shi;pshi(pshi<=0)=100;';fprintf(fid,'%s \n ',str);
%str='    ee = [ (log(psno))/NSno*fact1no;(pno-sno)/NSno*fact2no;(log(pshi))/NShi*fact1hi;(phi-shi)/NShi*fact2hi;constr;constrnew]/sqrt(32);';
str='    ee = [ (log(psno))/NSno*fact1no;(pno-sno)/NSno*fact2no;(log(pshi))/NShi*fact1hi;(phi-shi)/NShi*fact2hi;constr]/sqrt(30);';
fprintf(fid,'%s \n ',str);
str='end';fprintf(fid,'%s \n ',str);
 
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%changed!
if iii < 9
fname = ['exp_fitnessMODELone',num2str(iii),'.m']
fid=fopen(fname,'w');
str=['function ee=exp_fitnessMODEL',num2str(iii),'(k,alpha, x, s, ppno,pphi,ifile)'];
fprintf(fid,'%s \n ',str);
else
    if iii < 11
fname = ['exp_fitnessMODELone',num2str(iii+4),'.m']
fid=fopen(fname,'w');
str=['function ee=exp_fitnessMODEL',num2str(iii+4),'(k,alpha, x, s, ppno,pphi,ifile)'];
fprintf(fid,'%s \n ',str);
    else
fname = ['exp_fitnessMODELone',num2str(iii+6),'.m']
fid=fopen(fname,'w');
str=['function ee=exp_fitnessMODEL',num2str(iii+6),'(k,alpha, x, s, ppno,pphi,ifile)'];
fprintf(fid,'%s \n ',str);        
    end
end


for i=1:5
    str=['lambda',num2str(i),'=k(',num2str(i),');'];
    fprintf(fid,'%s \n ',str);
end
for i=1:4
    str=['A',num2str(i),'=k(',num2str(i+5),');'];
    fprintf(fid,'%s \n ',str);
end
str='A5=1-A1-A2-A3-A4;';
fprintf(fid,'%s \n ',str);
j=10;
for i=1:npars
   if ~ismember(i,idel)   
        str=['k',num2str(i),'=k(',num2str(j),');']; 
        fprintf(fid,'%s \n ',str);
        j=j+1;
   end
end

fprintf(fid,'%s \n ','L1 = lambda1 + lambda2 + lambda3 + lambda4 + lambda5;');
fprintf(fid,'%s \n ','L2 = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5;');
fprintf(fid,'%s \n ','L3 = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5;');
fprintf(fid,'%s \n ','L4 = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5;');
fprintf(fid,'%s \n ','L5 = lambda1*lambda2*lambda3*lambda4*lambda5;');

fprintf(fid,'%s \n ','S1 = A1*lambda1  + A2*lambda2   + A3*lambda3   + A4*lambda4   + A5*lambda5;'); 
fprintf(fid,'%s \n ','S2 = A1*lambda1^2+ A2*lambda2^2 + A3*lambda3^2 + A4*lambda4^2 + A5*lambda5^2;');
fprintf(fid,'%s \n ','S3 = A1*lambda1^3+ A2*lambda2^3 + A3*lambda3^3 + A4*lambda4^3 + A5*lambda5^3;');
fprintf(fid,'%s \n ','S4 = A1*lambda1^4+ A2*lambda2^4 + A3*lambda3^4 + A4*lambda4^4 + A5*lambda5^4;'); 

str=['LL1 = ', char(simplify(-d(5))),';'];
fprintf(fid,'%s \n ',str);
str=['LL2 = ',char(simplify(d(4))),';'];
fprintf(fid,'%s \n ',str);
str=['LL3 = ',char(simplify(-d(3))),';'];
fprintf(fid,'%s \n ',str);
str=['LL4 = ',char(simplify(d(2))),';'];
fprintf(fid,'%s \n ',str);
str=['LL5 = ',char(simplify(-d(1))),';'];
fprintf(fid,'%s \n ',str);
%str=['SS1 = ',char(simplify(SS(1))),';'];
%fprintf(fid,'%s \n ',str);
str=['SS2 = ',char(simplify(SS(2))),';'];
fprintf(fid,'%s \n ',str);
str=['SS3 = ',char(simplify(SS(3))),';'];
fprintf(fid,'%s \n ',str);
str=['SS4 = ',char(simplify(SS(4))),';'];
fprintf(fid,'%s \n ',str);

[N,D] = numden(res5.p1);
str=['pt=',char(D),';'];
fprintf(fid,'%s \n ',str);
str=['p1= ',char(simplify(res5.p1*D)),'/pt;'];
fprintf(fid,'%s \n ',str);
str=['p2= ',char(simplify(res5.p2*D)),'/pt;'];
fprintf(fid,'%s \n ',str);
str=['p3= ',char(simplify(res5.p3*D)),'/pt;'];
fprintf(fid,'%s \n ',str);
str=['p4= ',char(simplify(res5.p4*D)),'/pt;'];
fprintf(fid,'%s \n ',str);
str=['p5= ',char(simplify(res5.p5*D)),'/pt;'];
fprintf(fid,'%s \n ',str);

fprintf(fid,'%s \n ','SS1 = abs(A1*lambda1)+ abs(A2*lambda2) + abs(A3*lambda3) + abs(A4*lambda4) + abs(A5*lambda5);');


%%%%%%% changed!
fprintf(fid,'%s \n ','N = length(x);');
fprintf(fid,'%s \n ','sN= sqrt(N);'); 
fprintf(fid,'%s \n ','fact1=sqrt(1-alpha);');
fprintf(fid,'%s \n ','fact2=sqrt(alpha);');

fprintf(fid,'%s \n ','if ifile==2');
fprintf(fid,'%s \n ','    pp=pphi;');
fprintf(fid,'%s \n ','else');
fprintf(fid,'%s \n ','    pp=ppno;');
fprintf(fid,'%s \n ','end');
fprintf(fid,'%s \n ','constr=[(L1-LL1)/L1;(L2-LL2)/L2;(L3-LL3)/L3;(L4-LL4)/L4;(L5-LL5)/L5;S1/SS1;(S2-SS2)/S2;(S3-SS3)/S3;(S4-SS4)/S4;(p1-pp(1))/pp(1);(p2-pp(2))/pp(2);(p3-pp(3))/pp(3);(p4-pp(4))/pp(4);(p5-pp(5))/pp(5)];');
%%%%% 14 elements
%%%%%% changed!
fprintf(fid,'%s \n ','ps = A1*exp(lambda1*x) + A2*exp(lambda2*x)+A3*exp(lambda3*x)+A4*exp(lambda4*x)+A5*exp(lambda5*x);');
fprintf(fid,'%s \n ','if alpha == 1');
fprintf(fid,'%s \n ','    ee =  [(ps-s)/sN;constr]/sqrt(15);');                
fprintf(fid,'%s \n ','else');
fprintf(fid,'%s \n ','    pss = ps./s;pss(pss<=0)=100;');
fprintf(fid,'%s \n ','    ee = [(log(pss))/sN*fact1;(ps-s)/sN*fact2;constr]/sqrt(15);'); 
fprintf(fid,'%s \n ','end');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iii < 9
fname = ['objective_MODELone',num2str(iii),'.m']
else
    if iii < 11
fname = ['objective_MODELone',num2str(iii+4),'.m']
    else
fname = ['objective_MODELone',num2str(iii+6),'.m']        
    end
end
fid=fopen(fname,'w');



%%% display system
%str=['(L1 + ', char(simplify(d(5))),')/L1; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(-L2 + ',char(simplify(d(4))),')/L2; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(L3 + ',char(simplify(d(3))),')/L3; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(-L4 + ',char(simplify(d(2))),')/L4; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(L5 + ',char(simplify(d(1))),')/L5; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(S1 + ',char(simplify(-SS(1))),')/S1; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(S2 + ',char(simplify(-SS(2))),')/S2; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(S3 + ',char(simplify(-SS(3))),')/S3; ...'];
%fprintf(fid,'%s \n ',str);
%str=['(S4 + ',char(simplify(-SS(4))),')/S4; ...'];
%fprintf(fid,'%s \n ',str);


%%% display system
str=['LL1 = ', char(simplify(-d(5))),';'];
fprintf(fid,'%s \n ',str);
str=['LL2 = ',char(simplify(d(4))),';'];
fprintf(fid,'%s \n ',str);
str=['LL3 = ',char(simplify(-d(3))),';'];
fprintf(fid,'%s \n ',str);
str=['LL4 = ',char(simplify(d(2))),';'];
fprintf(fid,'%s \n ',str);
str=['LL5 = ',char(simplify(-d(1))),';'];
fprintf(fid,'%s \n ',str);
str=['SS1 = ',char(simplify(SS(1))),';'];
fprintf(fid,'%s \n ',str);
str=['SS2 = ',char(simplify(SS(2))),';'];
fprintf(fid,'%s \n ',str);
str=['SS3 = ',char(simplify(SS(3))),';'];
fprintf(fid,'%s \n ',str);
str=['SS4 = ',char(simplify(SS(4))),';'];
fprintf(fid,'%s \n ',str);


str=['p1= ',char(simplify(res5.p1))];
fprintf(fid,'%s \n ',str);
str=['p2= ',char(simplify(res5.p2))];
fprintf(fid,'%s \n ',str);
str=['p3= ',char(simplify(res5.p3))];
fprintf(fid,'%s \n ',str);
str=['p4= ',char(simplify(res5.p4))];
fprintf(fid,'%s \n ',str);
str=['p5= ',char(simplify(res5.p5))];
fprintf(fid,'%s \n ',str);
fclose(fid);




if iii < 9
fname = ['PON_Model',num2str(iii),'.m']
else
    if iii < 11
fname = ['PON_Model',num2str(iii+4),'.m']
    else
fname = ['PON_Model',num2str(iii+6),'.m']        
    end
end
fid=fopen(fname,'w');




ans=char(simplify(res5.p5));
ans=strrep(ans,'k10','kxo');
ans=strrep(ans,'k11','kxl');
ans=strrep(ans,'k12','kxxl');
for i=1:5
   ans=strrep(ans,['k',num2str(i)],['M.k',num2str(i),'(1)']);
end
ans=strrep(ans,'kxo','k10');
ans=strrep(ans,'kxl','k11');
ans=strrep(ans,'kxxl','k12');
anshi=ans;
ansno=ans;
for i=6:12
    anshi=strrep(anshi,['k',num2str(i)],['M.k',num2str(i),'hi(1)']);
    ansno=strrep(ansno,['k',num2str(i)],['M.k',num2str(i),'no(1)']);
end

str=['p5no= ',ansno,';'];
fprintf(fid,'%s \n ',str);
str=['p5hi= ',anshi,';'];
fprintf(fid,'%s \n ',str);

fclose(fid);

end %%% for iii

return




sign=-(-1)^abs(s-N);
if s < N
    str='S1 ; ...';
    fprintf(fid,'%s \n ',str);
    if sign==-1
        str=['(S2 + ',char(simplify(kini*c(4))),')/S2; ...'];
        fprintf(fid,'%s \n ',str);
        str=['(S3 + ',char(simplify(kini*(c(4)*L1+c(3)))),')/S3; ...'];
        fprintf(fid,'%s \n ',str);
        str=['(S4 + ',char(simplify(kini*(c(4)*(L1^2-L2)+c(3)*L1+c(2)))),')/S4; ...'];
        fprintf(fid,'%s \n ',str);    
    end
    if sign==1
        str=['(-S2 + ',char(simplify(kini*c(4))),')/S2; ...'];
        fprintf(fid,'%s \n ',str);
        str=['(-S3 + ',char(simplify(kini*(c(4)*L1+c(3)))),')/S3; ...'];
        fprintf(fid,'%s \n ',str);
        str=['(-S4 + ',char(simplify(kini*(c(4)*(L1^2-L2)+c(3)*L1+c(2)))),')/S4; ...'];
        fprintf(fid,'%s \n ',str);
    end
else   
    str='(-S1 + kini)/S1; ...';
    fprintf(fid,'%s \n ',str);
    str=['(-S2 + ',char(simplify(kini*(L1 + c(4)))),')/S2; ...'];
    fprintf(fid,'%s \n ',str);
    str=['(-S3 + ',char(simplify(kini*((L1^2-L2) + c(4)*L1+c(3)))),')/S3; ...'];
    fprintf(fid,'%s \n ',str);
    str=['(-S4 + ',char(simplify(kini*(L1^3 - 2*L2*L1 + L3+c(4)*(L1^2-L2)+c(3)*L1+c(2)))),')/S4; ...'];
    fprintf(fid,'%s \n ',str);
end
    



str=['(prob1 - ',char(simplify(res5.p1)),')/prob1; ...'];
fprintf(fid,'%s \n ',str);
str=['(prob2 - ',char(simplify(res5.p2)),')/prob2; ...'];
fprintf(fid,'%s \n ',str);
str=['(prob3 - ',char(simplify(res5.p3)),')/prob3; ...'];
fprintf(fid,'%s \n ',str);
str=['(prob4 - ',char(simplify(res5.p4)),')/prob4; ...'];
fprintf(fid,'%s \n ',str);
str=['(prob5 - ',char(simplify(res5.p5)),')/prob5]'];
fprintf(fid,'%s \n ',str);
fclose(fid);

return
%%%% problem to solve
res=solve([-L1 == d(5),L2 == d(4),-L3 == d(3),...
L4 == d(2), -L5 == d(1),...
S2 == -kini*c(4),S3==-kini*(c(4)*L1+c(3)),...
S4 == -kini*( c(4)*(L1^2-L2) + c(3)*L1+c(2))],k1,k2,k3,k4,k5,k6,k7,k8)

res1=solve([-L1 == d(5),L2 == d(4),-L3 == d(3),...
    S2 ==-kini*c(4),S3==-kini*(c(4)*L1+c(3))],k1,k2,k3,k4,k5)

res2=solve([L4 == d(2),-L5 == d(1),...
    S4 ==-kini*( c(4)*(L1^2-L2) + c(3)*L1+c(2))],k6,k7,k8)



return

res=solve(-L1 == d(5),L2 == d(4),-L3 == d(3),...
    L4 == d(2),-L5 == d(1), ...
    S2 == -k8*c(4),S3==-k8*(c(4)*L1+c(3)),...
    S4 == -k8*(c(4)*(L1^2-L2)+c(3)*L1+c(2)),k1,k2,k3,k4,k5,k6,k7,k8)


return

h1=simplify(expand(...
lambda1^5/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))+...
lambda2^5/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))+...
lambda3^5/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))+...
lambda4^5/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))+...
lambda5^5/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))...
));
h2=simplify(expand(...
lambda1^6/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))+...
lambda2^6/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))+...
lambda3^6/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))+...
lambda4^6/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))+...
lambda5^6/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))...
));
h3=simplify(expand(...
lambda1^7/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))+...
lambda2^7/((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))+...
lambda3^7/((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))+...
lambda4^7/((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))+...
lambda5^7/((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))...
));





return
simplify(expand(...
res2.C1  - subs(D3,lambda,lambda1)/...
((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-lambda4)*(lambda1-lambda5))...    
))
simplify(expand(...
res2.C2  - subs(D3,lambda,lambda2)/...
((lambda2-lambda1)*(lambda2-lambda3)*(lambda2-lambda4)*(lambda2-lambda5))...    
))
simplify(expand(...
res2.C3  - subs(D3,lambda,lambda3)/...
((lambda3-lambda1)*(lambda3-lambda2)*(lambda3-lambda4)*(lambda3-lambda5))...    
))
simplify(expand(...
res2.C4  - subs(D3,lambda,lambda4)/...
((lambda4-lambda1)*(lambda4-lambda2)*(lambda4-lambda3)*(lambda4-lambda5))...    
))
simplify(expand(...
res2.C5  - subs(D3,lambda,lambda5)/...
((lambda5-lambda1)*(lambda5-lambda2)*(lambda5-lambda3)*(lambda5-lambda4))...    
))





P1 = simplify((lambda1 - lambda2)*(lambda1 - lambda3)*(lambda1 - lambda4)*(lambda1 - lambda5)*res2.C1);
c1=coefficients(P1,lambda1);
P2 = simplify((lambda2 - lambda1)*(lambda2 - lambda3)*(lambda2 - lambda4)*res2.C2);
c2=coefficients(P2,lambda2);
P3 = simplify((lambda3 - lambda1)*(lambda3 - lambda2)*(lambda3 - lambda4)*res2.C3);
c3=coefficients(P3,lambda3);
P4 = simplify((lambda4 - lambda1)*(lambda4 - lambda2)*(lambda4 - lambda3)*res2.C4);
c4=coefficients(P4,lambda4);






return
