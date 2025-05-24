%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, August
%%%%% 2024
%%%%% needs short & long movies survival functions storead as s_functions.mat
%%%%% computes parameters of five states models; 
%%%%% in this version no tat and high tat datasets are
%%%%% merged and the parameters k1-5 are considered not dependent on tat 
%%%%  p4 and p5 are merged and imposed as a sum rather than individually
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% fits the MODELS 2,3,13  using the functions exp_fitnessMODELXmergeonekini_mergep4p5ponkinishort.m

clear all
close all

global alpha xs xl ss sl prob1 prob2 prob3 prob4 prob5;

%%%% occupation probabilities no tat
prob1=0.7; %%% 0.2
prob2=0.18; %%% 0.31
prob3=0.004; %%% 0.057
prob4=0.11; %%% 0.29
prob5=0.006; %%% 0.13
ppvalno=[prob1,prob2,prob3,prob4,prob5];

%%%% occupation probabilities high tat
prob1=0.21; %%% 0.2
prob2=0.31; %%% 0.31
prob3=0.06; %%% 0.057
prob4=0.29; %%% 0.29
prob5=0.13; %%% 0.13
ppvalhi=[prob1,prob2,prob3,prob4,prob5];
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
%fact=0.75; %%% [fact/(1-fact)]^2 is the relative importance of linear scale in the objective function
%alpha=0.9;  %%%% if alpha = 1 linear 
%%%%%%%%%%%%%%%%%%%%%%% Parameters 



    
for imodel = [2,3,13]
%for imodel = 3
%%%% MODEL DEPENDENT
DataFilePath0 = ['results_tat_MODEL',num2str(imodel),'_',num2str(100*alpha),'par_merge_onemoviekini_mergep4p5ponkinishort']; %%%% where to write results
textfilename1 = [DataFilePath0,'/resultsMODEL',num2str(imodel),'_',num2str(100*alpha),'_all.txt'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir(DataFilePath0);

%xlsfilename1 = [DataFilePath0,'/results5_2',num2str(100*alpha),'_all.xlsx'];



fid=fopen(textfilename1,'w');
str=['type, ','OBJ, ','lambda1no,','lambda2no,','lambda3no,',...
    'lambda4no,','lambda5no,','A1no,','A2no,','A3no,','A4no,','A5no,','S1no,',...
    'lambda1hi,','lambda2hi,','lambda3hi,',...
    'lambda4hi,','lambda5hi,','A1hi,','A2hi,','A3hi,','A4hi,','A5hi,','S1hi,',...
    'k1,','k2, ','k3, ','k4, ','k5, ','k6no, ','k7no, ','k8no, ','k9no, ','k10no, ',...
    'k6hi, ','k7hi, ','k8hi, ','k9hi, ','k10hi'];
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
%%%%%%%%%%%%%%%  initial guess of 2*9+10+5 = 33 pars
%%%%%%%%%%%%%%%  lambda1no,lambda2no,lambda3no,lambda4no,lmabda5no,lambda1hi,lambda2hi,lambda3hi,lambda4hi,lambda5hi,
%%%%%%%%%%%%%%%  A1no,A2no,A3no,A4no,A1hi,A2hi,A3hi,A4hi,k1,k2,k3,k4,k5,k6no,k7no,k8no,k9no,k10no,k6hi,k7hi,k8hi,k9hi,k10hi
%%%%%%%%%%%%%%%   
k00=[-0.1,-0.01,-0.001,-0.001,-0.001,0.1,0.1,0.1,0.1,...
-0.1,-0.01,-0.001,-0.001,-0.001,0.1,0.1,0.1,0.1,...    
0.001,0.001,0.01,0.01,0.01,0.01,0.01,0.01,...
0.01,0.1,0.01,0.01,0.01,0.01,0.1]; 
    

fmin = 1/100; fmax = 100;



iter = 5000;
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
    lf0 =  log(fmin)+log(fmax/fmin)*rand(1,5);
    k0(1:5)=k0(1:5).*exp(lf0); %%% lambdas
    k0(6:9)=rand(1,4); %%% As
    lf0 =  log(fmin)+log(fmax/fmin)*rand(1,5);
    k0(10:14)=k0(10:14).*exp(lf0); %%% lambdas
    k0(15:18)=rand(1,4); %%% As
    lf1 =  log(fmin)+log(fmax/fmin)*rand(1,15);
    k0(19:33)=k0(19:33).*exp(lf1); %%% ks 
    k0(34:35)=rand(1,2).*[ppvalno(4)+ppvalno(5),ppvalhi(4)+ppvalhi(5)];
switch imodel    
    case 2,    
    fobj=@(k)exp_fitnessMODEL2mergeonekini_mergep4p5ponkinishort(k, alpha,  xno, sno, xhi, shi, ppvalno,ppvalhi,kinino,kinihi);
    case 3,
    fobj=@(k)exp_fitnessMODEL3mergeonekini_mergep4p5ponkinishort(k, alpha,   xno, sno, xhi, shi, ppvalno,ppvalhi,kinino,kinihi);
    case 13,
    fobj=@(k)exp_fitnessMODEL13mergeonekini_mergep4p5ponkinishort(k, alpha,   xno, sno, xhi, shi, ppvalno,ppvalhi,kinino,kinihi);    
end


sobj= sum(fobj(k0));
if ~isinf(sobj) && ~isnan(sobj)
    [k, obj] = lsqnonlin(fobj,k0,...
    [-10,-10,-10,-10,-10,-10,-10,-10,-10, ...
    -10,-10,-10,-10,-10,-10,-10,-10,-10,...
    0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0],...
    [ 0,0,0,0,0, ...
      10,10,10,10, ...
      0,0,0,0,0,...
      10,10,10,10,...
      10,10,10,10,10,10,10,10,10,10,10,10,10,10,10],opts);
  


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
lambda5=k(5);
lambda1no=k(1);
lambda2no=k(2);
lambda3no=k(3);
lambda4no=k(4);
lambda5no=k(5);
L1no = lambda1 + lambda2 + lambda3 + lambda4 + lambda5;
L2no = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5;
L3no = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + ...
    lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5;
L4no = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5;
L5no = lambda1*lambda2*lambda3*lambda4*lambda5;
A1no=k(6);
A2no=k(7);
A3no=k(8);
A4no=k(9);
A5no=1-A1no-A2no-A3no-A4no;

S1no = A1no*lambda1  + A2no*lambda2   + A3no*lambda3   + A4no*lambda4   + A5no*lambda5;


tmax = max(xno);
ts=log(tmax)/100;
time=exp(-1:ts:log(tmax));
sfitno = A1no*exp(lambda1*time)+A2no*exp(lambda2*time)+A3no*exp(lambda3*time)+...
+A4no*exp(lambda4*time)+A5no*exp(lambda5*time);
timeno=time;

lambda1=k(10);
lambda2=k(11);
lambda3=k(12);
lambda4=k(13);
lambda5=k(14);
lambda1hi=k(10);
lambda2hi=k(11);
lambda3hi=k(12);
lambda4hi=k(13);
lambda5hi=k(14);
L1hi = lambda1 + lambda2 + lambda3 + lambda4 + lambda5;
L2hi = lambda1 * lambda2 + lambda1 * lambda3 + lambda1 * lambda4 + lambda1 * lambda5 + lambda2 * lambda3 + lambda2 * lambda4 + lambda2 * lambda5 + lambda3 * lambda4 + lambda3 * lambda5 + lambda4 * lambda5;
L3hi = lambda1 * lambda2 * lambda3 + lambda1 * lambda2 * lambda4 + lambda1 * lambda2 * lambda5 + lambda1 * lambda3 * lambda4 + lambda1 * lambda3 * lambda5 + lambda1 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 + ...
    lambda2 * lambda3 * lambda5 + lambda2 * lambda4 * lambda5 + lambda3 * lambda4 * lambda5;
L4hi = lambda1 * lambda2 * lambda3 * lambda4 + lambda1 * lambda2 * lambda3 * lambda5 + lambda1 * lambda2 * lambda4 * lambda5 + lambda1 * lambda3 * lambda4 * lambda5 + lambda2 * lambda3 * lambda4 * lambda5;
L5hi = lambda1*lambda2*lambda3*lambda4*lambda5;
A1hi=k(15);
A2hi=k(16);
A3hi=k(17);
A4hi=k(18);
A5hi=1-A1hi-A2hi-A3hi-A4hi;

S1hi = A1hi*lambda1  + A2hi*lambda2   + A3hi*lambda3   + A4hi*lambda4   + A5hi*lambda5;

tmax = max(xhi);
ts=log(tmax)/100;
time=exp(-1:ts:log(tmax));
sfithi = A1hi*exp(lambda1*time)+A2hi*exp(lambda2*time)+A3hi*exp(lambda3*time)+...
+A4hi*exp(lambda4*time)+A5hi*exp(lambda5*time);
timehi=time;

k1=k(19);
k2=k(20);
k3=k(21);
k4=k(22);
k5=k(23);
k6no=k(24);
k7no=k(25);
k8no=k(26);
k9no=k(27);
k10no=k(28);
k6hi=k(29);
k7hi=k(30);
k8hi=k(31);
k9hi=k(32);
k10hi=k(33);
%p5no=k(34);
%p5hi=k(35);

if iii==1


%%%% show fit    
h=figure(1000);
subplot(2,1,1)
hold off
loglog(xno,sno,'og')
hold on
loglog(timeno,sfitno,'k','Linewidth',2)
xlabel('Time [s]','Fontsize',fsz)
ylabel('Freq','Fontsize',fsz)
subplot(2,1,2)
hold off
loglog(xhi,shi,'og')
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
num2str(lambda4no,precision),', ',num2str(lambda5no,precision),', ',num2str(A1no,precision),', ',...
num2str(A2no,precision),', ',num2str(A3no,precision),', ',num2str(A4no,precision),', ',...
num2str(A5no,precision),', ',num2str(S1no,precision),', ',...
num2str(lambda1hi,precision),', ',num2str(lambda2hi,precision),', ',num2str(lambda3hi,precision),', ',...
num2str(lambda4hi,precision),', ',num2str(lambda5hi,precision),', ',num2str(A1hi,precision),', ',...
num2str(A2hi,precision),', ',num2str(A3hi,precision),', ',num2str(A4hi,precision),', ',...
num2str(A5hi,precision),', ',num2str(S1hi,precision),', ',...
num2str(k1,precision),', ',num2str(k2,precision),', ',num2str(k3,precision), ', ',...
num2str(k4,precision),', ',num2str(k5,precision),', ',num2str(k6no,precision),', ',...
num2str(k7no,precision),', ',num2str(k8no,precision),', ',num2str(k9no,precision),', ',num2str(k10no,precision),', ', ...
num2str(k6hi,precision),', ',...
num2str(k7hi,precision),', ',num2str(k8hi,precision),', ',num2str(k9hi,precision),', ',num2str(k10hi,precision)];
fprintf(fid,'%s \n ',str);
else
str=['suboptimal ',', ',num2str(o),', ', ...
num2str(lambda1no,precision),', ',num2str(lambda2no,precision),', ',num2str(lambda3no,precision),', ',...
num2str(lambda4no,precision),', ',num2str(lambda5no,precision),', ',num2str(A1no,precision),', ',...
num2str(A2no,precision),', ',num2str(A3no,precision),', ',num2str(A4no,precision),', ',...
num2str(A5no,precision),', ',num2str(S1no,precision),', ',...
num2str(lambda1hi,precision),', ',num2str(lambda2hi,precision),', ',num2str(lambda3hi,precision),', ',...
num2str(lambda4hi,precision),', ',num2str(lambda5hi,precision),', ',num2str(A1hi,precision),', ',...
num2str(A2hi,precision),', ',num2str(A3hi,precision),', ',num2str(A4hi,precision),', ',...
num2str(A5hi,precision),', ',num2str(S1hi,precision),', ',...
num2str(k1,precision),', ',num2str(k2,precision),', ',num2str(k3,precision), ', ',...
num2str(k4,precision),', ',num2str(k5,precision),', ',num2str(k6no,precision),', ',...
num2str(k7no,precision),', ',num2str(k8no,precision),', ',num2str(k9no,precision),', ',num2str(k10no,precision),', ', ...
num2str(k6hi,precision),', ',...
num2str(k7hi,precision),', ',num2str(k8hi,precision),', ',num2str(k9hi,precision),', ',num2str(k10hi,precision)];
fprintf(fid,'%s \n ',str);    
end



end %%% iii



fclose(fid);

end %%% imodel



