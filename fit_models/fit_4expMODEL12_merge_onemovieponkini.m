%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, August
%%%%% 2024
%%%%% needs short & long movies survival functions storead as s_functions.mat
%%%%% computes parameters of four states models; 
%%%%% in this version no tat and high tat datasets are
%%%%% merged and the parameters k1-5 are considered not dependent on tat 
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% fits the MODELS 12  using the functions exp_fitnessMODELXmergeoneponkini.m

clear all
close all

global alpha xs xl ss sl prob1 prob2 prob3 prob4 prob5;

%%%% occupation probabilities no tat
prob1=0.7; %%% 
prob2=0.18; %%% 
prob3=0.004; %%% 
prob4=0.11 + 0.006; %%% 

ppvalno=[prob1,prob2,prob3,prob4];

%%%% occupation probabilities high tat
prob1=0.21; %%% 0.2
prob2=0.31; %%% 0.31
prob3=0.06; %%% 0.057
prob4=0.29 + 0.13; %%% 0.13
ppvalhi=[prob1,prob2,prob3,prob4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% kini %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kinino=0.0022138;
kinihi=0.23962;


alpha = 0.3; %%% relative importance linear
thresh=1200; %%% threshold long movie    
%%% parameters
fsz=16; msz=5;lw=1;
outliers_long=1;
outliers_short=0;
visualize = 1; 
%%%%%%%%%%%%%%%%%%%%%%% 



    
for imodel = 12:12 

%%%% MODEL DEPENDENT
DataFilePath0 = ['results_tat_MODEL',num2str(imodel),'_',num2str(100*alpha),'par_merge_onemoviekini']; %%%% where to write results
textfilename1 = [DataFilePath0,'/resultsMODEL',num2str(imodel),'_',num2str(100*alpha),'_all.txt'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir(DataFilePath0);

%xlsfilename1 = [DataFilePath0,'/results5_2',num2str(100*alpha),'_all.xlsx'];



fid=fopen(textfilename1,'w');
str=['type, ','OBJ, ','lambda1no,','lambda2no,','lambda3no,',...
    'lambda4no,','A1no,','A2no,','A3no,','A4no,','S1no,',...
    'lambda1hi,','lambda2hi,','lambda3hi,',...
    'lambda4hi,','A1hi,','A2hi,','A3hi,','A4hi,','S1hi,',...
    'k1,','k2, ','k3, ','k4, ','k5, ','k9no, ','k10no, ',...
    'k11no, ','k12no, ','k9hi, ','k10hi, ','k11hi, ','k12hi'];


fprintf(fid,'%s \n ',str);
%line1=1;
%xlswrite(xlsfilename1,{'Fname','OBJ','lambda1','lambda2','lambda3',...
%    'lambda4','lambda5','A1','A2','A3','A4','A5','S1','S2','S3','S4','k1','k2','k3','k4','k5','k6','k7','k8',...
%    'p1','p2','p3','p4','p5','p6'},1,'A1');



%%%%% changed!
%%%%% one movie survival functions stored in  
%%%%% xhi shi xno sno 
%%%% need to modify exp_fitnessMODEL1merge
load('s_functions.mat');



opts = optimoptions(@lsqnonlin,'Display','none','TolFun', 1e-8,'MaxIter',1e3, 'MaxFunEvals',1e6,'TolX', 1e-10);
%%%%%%%%%%%%%%%  initial guess of 27 pars
%%%%%%%%%%%%%%%  lambda1no,lambda2no,lambda3no,lambda4no,lambda1hi,lambda2hi,lambda3hi,lambda4hi,
%%%%%%%%%%%%%%%  A1no,A2no,A3no,A4no,A1hi,A2hi,A3hi,A4hi,
%                k1,k2,k3,k4,k5,k9no,k10no,k9hi,k10hi
%                k1,k2,k3,k4,k5,k9no,k10no,k11no,k12no,k9hi,k10hi,k11hi,k12hi
%%%%%%%%%%%%%%%   
k00=[-0.1,-0.01,-0.001,-0.001,0.1,0.1,0.1,...
-0.1,-0.01,-0.001,-0.001,0.1,0.1,0.1,...    
0.000203162,5.54388E-05,0.000866101,0.003286996,1.040101851,...
0.001018218,0.023415801,1e-4,1e-4,...
0.000130657,0.219531747,1-5,1e-5];
 

fmin = 1/2; fmax = 2;



iter = 5000;
%iter=40;
%iter =20;
Obj = zeros(1,iter);
Pars = cell(1,iter);

tic
%%%%%% parallel computing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mycluster = parcluster;   %%% create cluster %%%%%%%%%%%%%%%%%%%%%
set(mycluster,'NumWorkers',20) ; %%% change max number of workers
parpool(mycluster,20); %%% open a pool of workers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor MC=1:iter
%for MC=1:iter
    %MC
    k0=k00;
    lf0 =  log(fmin)+log(fmax/fmin)*rand(1,4);
    k0(1:4)=k0(1:4).*exp(lf0); %%% lambdas
    k0(5:7)=rand(1,3); %%% As
    lf0 =  log(fmin)+log(fmax/fmin)*rand(1,4);
    k0(8:11)=k0(8:11).*exp(lf0); %%% lambdas
    k0(12:14)=rand(1,3); %%% As
    lf1 =  log(fmin)+log(fmax/fmin)*rand(1,13);
    k0(15:27)=k0(15:27).*exp(lf1); %%% ks 

 fobj=@(k)exp_fitnessMODEL12mergeoneponkini(k, alpha,   xno, sno, xhi, shi, ppvalno,ppvalhi,kinino,kinihi);



sobj= sum(fobj(k0));
if ~isinf(sobj) && ~isnan(sobj)
    [k, obj] = lsqnonlin(fobj,k0,...
    [-10,-10,-10,-10,-10,-10,-10, ...
    -10,-10,-10,-10,-10,-10,-10,...
    0,0,0,0,0,0,0,0,...
    0,0,0,0,0],...
    [0,0,0,0, ...
    10,10,10, ...
     0,0,0,0,...
     10,10,10,...
     10,10,10,10,10,10,10,10,10,10,10,10,10],opts);
  


%if obj < omin
%    omin = obj
%    kmin = k;
%end



%%%%% store optimization result
    Pars{MC}=k;
    Obj(MC)=obj;
else
    Obj(MC)=0;
end
end  %%% for MC


delete(gcp('nocreate'))
toc
omin = 100;
ind = find(Obj>0);
Obj=Obj(ind);
Pars=Pars(ind);


[omin,ii]=min(Obj);
kmin = Pars{ii};

%%%%%%select pars
isel = find(Obj < omin*5); %%% overflow 5 fold

for iii=1:length(isel)+1

if iii==1    
    k = kmin;o=omin;
else
    k = Pars{isel(iii-1)};o=Obj(isel(iii-1));
end


lambda1=k(1);
lambda2=k(2);
lambda3=k(3);
lambda4=k(4);
%lambda5=k(5);
lambda1no=k(1);
lambda2no=k(2);
lambda3no=k(3);
lambda4no=k(4);
%lambda5no=k(5);
L1no = lambda1 + lambda2 + lambda3 + lambda4 ;
L2no = lambda1*lambda2 + lambda1*lambda3 + lambda1*lambda4 + lambda2*lambda3 + lambda2*lambda4 + lambda3*lambda4;
L3no = lambda1*lambda2*lambda3 + lambda1*lambda2*lambda4 + lambda1*lambda3*lambda4   +   lambda2*lambda3*lambda4;
L4no = lambda1 * lambda2 * lambda3 * lambda4 ;
%L5no = lambda1*lambda2*lambda3*lambda4*lambda5;
A1no=k(5);
A2no=k(6);
A3no=k(7);
%A4no=k(9);
A4no=1-A1no-A2no-A3no;

S1no = A1no*lambda1  + A2no*lambda2   + A3no*lambda3   + A4no*lambda4   ;


tmax = max(xno);
ts=log(tmax)/100;
time=exp(-1:ts:log(tmax));
sfitno = A1no*exp(lambda1*time)+A2no*exp(lambda2*time)+A3no*exp(lambda3*time)+...
+A4no*exp(lambda4*time);
timeno=time;

lambda1=k(8);
lambda2=k(9);
lambda3=k(10);
lambda4=k(11);
%lambda5=k(14);
lambda1hi=k(8);
lambda2hi=k(9);
lambda3hi=k(10);
lambda4hi=k(11);
%lambda5hi=k(14);
L1hi = lambda1 + lambda2 + lambda3 + lambda4;
L2hi = lambda1*lambda2 + lambda1*lambda3 + lambda1*lambda4 + lambda2*lambda3 + lambda2*lambda4 + lambda3*lambda4;
L3hi = lambda1*lambda2*lambda3 + lambda1*lambda2*lambda4 + lambda1*lambda3*lambda4   +   lambda2*lambda3*lambda4;
L4hi = lambda1 * lambda2 * lambda3 * lambda4 ;
%L5hi = lambda1*lambda2*lambda3*lambda4*lambda5;
A1hi=k(12);
A2hi=k(13);
A3hi=k(14);
%A4hi=k(18);
A4hi=1-A1hi-A2hi-A3hi;

S1hi = A1hi*lambda1  + A2hi*lambda2   + A3hi*lambda3   + A4hi*lambda4 ;

tmax = max(xhi);
ts=log(tmax)/100;
time=exp(-1:ts:log(tmax));
sfithi = A1hi*exp(lambda1*time)+A2hi*exp(lambda2*time)+A3hi*exp(lambda3*time)+...
+A4hi*exp(lambda4*time);
timehi=time;

k1=k(15); 
k2=k(16); 
 k3=k(17); 
 k4=k(18); 
 k5=k(19); 
 k9no=k(20); 
 k10no=k(21); 
 k11no=k(22); 
 k12no=k(23); 
 k9hi=k(24); 
 k10hi=k(25); 
 k11hi=k(26); 
 k12hi=k(27); 


if iii==1


%%%% show fit    
h=figure(1000);
subplot(2,1,1)
hold off
loglog(xno,sno,'or')
hold on
loglog(timeno,sfitno,'k','Linewidth',2)
xlabel('Time [s]','Fontsize',fsz)
ylabel('Freq','Fontsize',fsz)
subplot(2,1,2)
hold off
loglog(xhi,shi,'or')
hold on
loglog(timehi,sfithi,'k','Linewidth',2)
xlabel('Time [s]','Fontsize',fsz)
ylabel('Freq','Fontsize',fsz)



figfile=[DataFilePath0,'/figure_optimal.png'];
print(h,'-dpng',figfile)
end

precision=32;

%%%%%% write par values in file
if iii==1
str=['optimal ',', ',num2str(o),', ', ...
num2str(lambda1no,precision),', ',num2str(lambda2no,precision),', ',num2str(lambda3no,precision),', ',...
num2str(lambda4no,precision),', ',num2str(A1no,precision),', ',...
num2str(A2no,precision),', ',num2str(A3no,precision),', ',...
num2str(A4no,precision),', ',num2str(S1no,precision),', ',...
num2str(lambda1hi,precision),', ',num2str(lambda2hi,precision),', ',num2str(lambda3hi,precision),', ',...
num2str(lambda4hi,precision),', ',num2str(A1hi,precision),', ',...
num2str(A2hi,precision),', ',num2str(A3hi,precision),', ',num2str(A4hi,precision),', ',...
num2str(S1hi,precision),', ',...
num2str(k1,precision),', ',num2str(k2,precision),', ',num2str(k3,precision), ', ',...
num2str(k4,precision),', ',num2str(k5,precision),', ',...
num2str(k9no,precision),', ',num2str(k10no,precision),', ', ...
num2str(k11no,precision),', ',num2str(k12no,precision),', ', ...
num2str(k9hi,precision),', ',num2str(k10hi,precision),', ', ...
num2str(k11hi,precision),', ',num2str(k12hi,precision)];
fprintf(fid,'%s \n ',str);
else
str=['suboptimal ',', ',num2str(o),', ', ...
num2str(lambda1no,precision),', ',num2str(lambda2no,precision),', ',num2str(lambda3no,precision),', ',...
num2str(lambda4no,precision),', ',num2str(A1no,precision),', ',...
num2str(A2no,precision),', ',num2str(A3no,precision),', ',num2str(A4no,precision),', ',...
num2str(S1no,precision),', ',...
num2str(lambda1hi,precision),', ',num2str(lambda2hi,precision),', ',num2str(lambda3hi,precision),', ',...
num2str(lambda4hi,precision),', ',num2str(A1hi,precision),', ',...
num2str(A2hi,precision),', ',num2str(A3hi,precision),', ',num2str(A4hi,precision),', ',...
num2str(S1hi,precision),', ',...
num2str(k1,precision),', ',num2str(k2,precision),', ',num2str(k3,precision), ', ',...
num2str(k4,precision),', ',num2str(k5,precision),', ',...
num2str(k9no,precision),', ',num2str(k10no,precision),', ', ...
num2str(k11no,precision),', ',num2str(k12no,precision),', ', ...
num2str(k9hi,precision),', ',num2str(k10hi,precision),', ', ...
num2str(k11hi,precision),', ',num2str(k12hi,precision)];

fprintf(fid,'%s \n ',str);    
end



end %%% iii

fclose(fid);

end %%%% imodel























