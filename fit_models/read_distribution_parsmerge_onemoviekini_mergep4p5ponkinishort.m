%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, August
%%%%% 2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% read results different fits, generate figures %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% needs data from directories %%%%%%%%%
%directory='results_tat_MODEL1_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL2_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL3_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL4_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL5_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL6_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL7_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL8_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL9_30par_merge_onemovieponkini\';
%directory='results_tat_MODEL10_30par_merge_onemovieponkini\';
%directory='results_tat_MODEL11_30par_merge_onemovieponkini\';
%directory='results_tat_MODEL12_30par_merge_onemoviekini\';
%directory='results_tat_MODEL13_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL14_30par_merge_onemoviekini_mergep4p5ponkinishortbis\';
%directory='results_tat_MODEL15_30par_merge_onemoviekini\';
%directory='results_tat_MODEL16_30par_merge_onemoviekini\';
%directory='results_tat_MODEL17_30par_merge_onemoviekini_mergep4p5ponkinishort\';
%directory='results_tat_MODEL18_30par_merge_onemoviekini_mergep4p5ponkinishort\';

clear all
close all
load('s_functions.mat');

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


%%%% occupation probabilities no tat
prob1=0.7; %%% 
prob2=0.18; %%% 
prob3=0.004; %%% 
prob4=0.11 + 0.006; %%% 

ppvalno4=[prob1,prob2,prob3,prob4];

%%%% occupation probabilities high tat
prob1=0.21; %%% 0.2
prob2=0.31; %%% 0.31
prob3=0.06; %%% 0.057
prob4=0.29 + 0.13; %%% 0.13
ppvalhi4=[prob1,prob2,prob3,prob4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




dirwrite = 'results_mergekinimergep4p5ponkiniallshort45/';
mkdir(dirwrite);
xfile=[dirwrite,'all_storebis.xlsx'];

iline=1;
xlswrite(xfile,{'model','obj','logvarobj','k1','k2','k3','k4','k5','k6no','k7no','k8no','k9no','k10no','k11no','k12no',...
    'k6hi','k7hi','k8hi','k9hi','k10hi','k11hi','k12hi','p5no','p5hi',...
    'lambda1no','lambda2no','lambda3no','lambda4no','lambda5no',...
    'lambda1hi','lambda2hi','lambda3hi','lambda4hi','lambda5hi',...
    'A1no','A2no','A3no','A4no','A5no',...
    'A1hi','A2hi','A3hi','A4hi','A5hi',...
    'logvark1','logvark2','logvark3','logvark4','logvark5','logvark6no','logvark7no','logvark8no','logvark9no','logvark10no','logvark11no','logvark12no',...
    'logvark6hi','logvark7hi','logvark8hi','logvark9hi','logvark10hi','logvark11hi','logvark12hi','varp5no','varp5hi',....
    'logvarlambda1no','logvarlambda2no','logvarlambda3no','logvarlambda4no','logvarlambda5no',...
    'logvarlambda1hi','logvarlambda2hi','logvarlambda3hi','logvarlambda4hi','logvarlambda5hi',...
    'varA1no','varA2no','varA3no','varA4no','varA5no',...
    'varA1hi','varA2hi','varA3hi','varA4hi','varA5hi','ponkinino','ponkinihi',...
    'p1no','p2no','p3no','p4no','p5no','p1hi','p2hi','p3hi','p4hi','p5hi','Obj'},1,['A',num2str(iline)]);
iline= iline+1;
%for jjj=[1,2,3,4,5,6,7,8,13,14,17,18]    %%%% different models 
for jjj=1:18    
    switch jjj
        case 1,
directory='results_tat_MODEL1_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL1_30_all.txt'];
pfile = [dirwrite,'MODEL1_30.pdf'];
namem = 'MODEL1';
nstates=5;
        case 2,
directory='results_tat_MODEL2_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL2_30_all.txt'];
pfile = [dirwrite,'MODEL2_30.pdf'];
%directory='results_tat_MODEL2_95par_merge_onemovie\';
%results_file= [directory,'resultsMODEL2_95_all.txt'];
%pfile = [dirwrite,'MODEL2_95.pdf'];
namem = 'MODEL2';
nstates=5;
        case 3,
directory='results_tat_MODEL3_30par_merge_onemoviekini_mergep4p5ponkinishort\';
 results_file= [directory,'resultsMODEL3_30_all.txt'];
pfile = [dirwrite,'MODEL3_30.pdf'];
%directory='results_tat_MODEL3_95par_merge_onemovie\';
% results_file= [directory,'resultsMODEL3_95_all.txt'];
%pfile = [dirwrite,'MODEL3_95.pdf'];
namem = 'MODEL3';
nstates=5;
        case 4,
directory='results_tat_MODEL4_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL4_30_all.txt'];
pfile = [dirwrite,'MODEL4_30.pdf'];
namem = 'MODEL4';
nstates=5;
        case 5,
directory='results_tat_MODEL5_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL5_30_all.txt'];
pfile = [dirwrite,'MODEL5_30.pdf'];
namem = 'MODEL5';
nstates=5;
        case 6,
directory='results_tat_MODEL6_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL6_30_all.txt'];
pfile = [dirwrite,'MODEL6_30.pdf'];
namem = 'MODEL6';
nstates=5;
        case 7,
directory='results_tat_MODEL7_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL7_30_all.txt'];
pfile = [dirwrite,'MODEL7_30.pdf'];
namem = 'MODEL7';
nstates=5;
        case 8,
directory='results_tat_MODEL8_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL8_30_all.txt'];
pfile = [dirwrite,'MODEL8_30.pdf'];
namem = 'MODEL8';
nstates=5;
        case 9,
directory='results_tat_MODEL9_30par_merge_onemovieponkini\';
results_file= [directory,'resultsMODEL9_30_all.txt'];
pfile = [dirwrite,'MODEL9_30.pdf'];
namem = 'MODEL9';
nstates=4;
        case 10,
directory='results_tat_MODEL10_30par_merge_onemovieponkini\';
results_file= [directory,'resultsMODEL10_30_all.txt'];
pfile = [dirwrite,'MODEL10_30.pdf'];
namem = 'MODEL10';
nstates=4;
        case 11,
directory='results_tat_MODEL11_30par_merge_onemovieponkini\';
results_file= [directory,'resultsMODEL11_30_all.txt'];
pfile = [dirwrite,'MODEL11_30.pdf'];
namem = 'MODEL11';
nstates=4;
        case 12,
directory='results_tat_MODEL12_30par_merge_onemoviekini\';
results_file= [directory,'resultsMODEL12_30_all.txt'];
pfile = [dirwrite,'MODEL12_30.pdf'];
namem = 'MODEL12';
nstates=4;
        case 13,
directory='results_tat_MODEL13_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL13_30_all.txt'];
pfile = [dirwrite,'MODEL13_30.pdf'];
namem = 'MODEL13';
nstates=5;
        case 14,
directory='results_tat_MODEL14_30par_merge_onemoviekini_mergep4p5ponkinishortbis\';
results_file= [directory,'resultsMODEL14_30_all.txt'];
pfile = [dirwrite,'MODEL14_30.pdf'];
namem = 'MODEL14';
nstates=5;
       case 15,
directory='results_tat_MODEL15_30par_merge_onemoviekini\';
results_file= [directory,'resultsMODEL15_30_all.txt'];
pfile = [dirwrite,'MODEL15_30.pdf'];
namem = 'MODEL15';
nstates=4;
       case 16
directory='results_tat_MODEL16_30par_merge_onemoviekini\';
results_file= [directory,'resultsMODEL16_30_all.txt'];
pfile = [dirwrite,'MODEL16_30.pdf'];
namem = 'MODEL16';
nstates=4;
        case 17,
directory='results_tat_MODEL17_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL17_30_all.txt'];
pfile = [dirwrite,'MODEL17_30.pdf'];
namem = 'MODEL17';
nstates=5;
        case 18,
directory='results_tat_MODEL18_30par_merge_onemoviekini_mergep4p5ponkinishort\';
results_file= [directory,'resultsMODEL18_30_all.txt'];
pfile = [dirwrite,'MODEL18_30.pdf'];
namem = 'MODEL18';
nstates=5;
    end
   
opar   = zeros(1,35);   %%%% optimal par values 
logvar = zeros(1,35);  %%%% logvariance of par values
    
M=readtable(results_file,'Delimiter',',');    

n=size(M);
if nstates ==5
%%%% sort lambda values
%type, OBJ, lambda1no,lambda2no,lambda3no,lambda4no,lambda5no,A1no,A2no,A3no,A4no,A5no,S1no,lambda1hi,lambda2hi,lambda3hi,lambda4hi,lambda5hi,A1hi,A2hi,A3hi,A4hi,A5hi
% lambdano 3:7 Ano 8:12 lambdahi 14:18 Ahi 19:23
for ii=1:n(1)
    lno = table2array(M(ii,3:7));
    lnoc = table2cell(M(ii,3:7));anoc=table2cell(M(ii,8:12));
    [lnosorted,IX]=sort(lno);
    M(ii,3:7)=lnoc(IX);
    M(ii,8:12)=anoc(IX);
    lhi = table2array(M(ii,14:18));
    lhic = table2cell(M(ii,14:18));ahic=table2cell(M(ii,19:23));
    [lhisorted,IX]=sort(lhi);
    M(ii,14:18)=lhic(IX);
    M(ii,19:23)=ahic(IX);
end
else
for ii=1:n(1)
    lno = table2array(M(ii,3:6));
    lnoc = table2cell(M(ii,3:6));anoc=table2cell(M(ii,7:10));
    [lnosorted,IX]=sort(lno);
    M(ii,3:6)=lnoc(IX);
    M(ii,7:10)=anoc(IX);
    lhi = table2array(M(ii,12:15));
    lhic = table2cell(M(ii,12:15));ahic=table2cell(M(ii,16:19));
    [lhisorted,IX]=sort(lhi);
    M(ii,12:15)=lhic(IX);
    M(ii,16:19)=ahic(IX);
end    
    
end




hf=figure(jjj) %%% one figure per model

of = M.OBJ; %%% objective function
ofmin = of(1); %%%% best value 

%%%% overload for computing the logvariance
over = 5;
%over = 3;
%over=  2;
isel = find(of < ofmin *over);
logvarobj = var(log(of(isel)));


switch jjj
    case 1, dist1 = of(isel);
    case 2, dist2 = of(isel); 
    case 3, dist3 = of(isel); 
    case 4, dist4 = of(isel);
    case 5, dist5 = of(isel);
    case 6, dist6 = of(isel);
    case 7, dist7 = of(isel);
    case 8, dist8 = of(isel);  
    case 9, dist9 = of(isel);
    case 10, dist10 = of(isel);
    case 11, dist11 = of(isel);
    case 12, dist12 = of(isel);
    case 13, dist13 = of(isel);
    case 14, dist14 = of(isel);
    case 15, dist15 = of(isel);
    case 16, dist16 = of(isel);
    case 17, dist17 = of(isel);
    case 18, dist18 = of(isel);    
end



%%%% build data (set of par values), the first value is the optimal
for iii=1:41
    switch iii
        case 1,  data=M.k1;   parname= 'k1';
        case 2,  data=M.k2;   parname= 'k2';
        case 3,  data=M.k3;   parname= 'k3';
        case 4,  data = M.k4; parname= 'k4';
        case 5,  data = M.k5; parname= 'k5';
        case 6,  if any(strcmp(M.Properties.VariableNames,'k6no')), data = M.k6no; else, data =zeros(size(M.k1)), end; parname= 'k6no';
        case 7,  if any(strcmp(M.Properties.VariableNames,'k7no')), data = M.k7no; else, data =zeros(size(M.k1)), end;parname= 'k7no';
        case 8,  if any(strcmp(M.Properties.VariableNames,'k8no')), data = M.k8no; else, data =zeros(size(M.k1)), end;parname= 'k8no';
        case 9,  if any(strcmp(M.Properties.VariableNames,'k9no')), data = M.k9no; else data =zeros(size(M.k1)), end;parname= 'k9no';
        case 10,  if any(strcmp(M.Properties.VariableNames,'k10no')), data = M.k10no; else, data =zeros(size(M.k1)), end;parname= 'k10no';
        case 11,  if any(strcmp(M.Properties.VariableNames,'k11no')), data = M.k11no; else, data =zeros(size(M.k1)), end;parname= 'k11no';
        case 12,  if any(strcmp(M.Properties.VariableNames,'k12no')), data = M.k12no; else, data =zeros(size(M.k1)), end;parname= 'k12no';
        case 13,  if any(strcmp(M.Properties.VariableNames,'k6hi')), data = M.k6hi; else, data =zeros(size(M.k1)), end;parname= 'k6hi';
        case 14,  if any(strcmp(M.Properties.VariableNames,'k7hi')), data = M.k7hi; else, data =zeros(size(M.k1)), end;parname= 'k7hi';
        case 15,  if any(strcmp(M.Properties.VariableNames,'k8hi')), data = M.k8hi; else, data =zeros(size(M.k1)), end;parname= 'k8hi';
        case 16,  if any(strcmp(M.Properties.VariableNames,'k9hi')), data = M.k9hi; else, data =zeros(size(M.k1)), end;parname= 'k9hi';
        case 17,  if any(strcmp(M.Properties.VariableNames,'k10hi')), data = M.k10hi; else, data =zeros(size(M.k1)), end;parname= 'k10hi';
        case 18,  if any(strcmp(M.Properties.VariableNames,'k11hi')), data = M.k11hi; else, data =zeros(size(M.k1)), end;parname= 'k11hi';   
        case 19,  if any(strcmp(M.Properties.VariableNames,'k12hi')), data = M.k12hi; else, data =zeros(size(M.k1)), end;parname= 'k12hi';
        case 20,  if any(strcmp(M.Properties.VariableNames,'p5no')), data = M.p5no; else, data =zeros(size(M.k1)), end;parname= 'p5no';
        case 21,  if any(strcmp(M.Properties.VariableNames,'p5hi')), data = M.p5hi; else, data =zeros(size(M.k1)), end;parname= 'p5hi';
        case 22,  data=M.lambda1no;parname= 'lambda1no';    
        case 23,  data=M.lambda2no;parname= 'lambda2no';
        case 24,  data=M.lambda3no;parname= 'lamnda3no';
        case 25,  data=M.lambda4no;parname= 'lambda4no';
        case 26,  if any(strcmp(M.Properties.VariableNames,'lambda5no')), data = M.lambda5no; else, data =zeros(size(M.k1)), end;parname= 'lambda5no';
        case 27,  data=M.lambda1hi;parname= 'lambda1hi';    
        case 28,  data=M.lambda2hi;parname= 'lambda2hi';
        case 29,  data=M.lambda3hi;parname= 'lambda3hi';
        case 30,  data=M.lambda4hi;parname= 'lambda4hi';
        case 31,  if any(strcmp(M.Properties.VariableNames,'lambda5hi')), data = M.lambda5hi; else, data =zeros(size(M.k1)), end;parname= 'lambda5hi';
        case 32,  data=M.A1no; parname= 'A1no';   
        case 33,  data=M.A2no; parname= 'A2no';
        case 34,  data=M.A3no; parname= 'A3no';
        case 35,  data=M.A4no; parname= 'A4no';
        case 36,  if any(strcmp(M.Properties.VariableNames,'A5no')), data = M.A5no; else, data =zeros(size(M.k1)), end; parname= 'A5no';
        case 37,  data=M.A1hi; parname= 'A1hi';           
        case 38,  data=M.A2hi; parname= 'A2hi';
        case 39,  data=M.A3hi; parname= 'A3hi';
        case 40,  data=M.A4hi; parname= 'A4hi';
        case 41,  if any(strcmp(M.Properties.VariableNames,'A5hi')), data = M.A5hi; else, data =zeros(size(M.k1)), end; parname= 'A5hi';
    end
    


   
   if abs(data(1)) > 0
        if iii==20 | iii==21 | iii> 31
            vals = data(isel);
            oval = data(1);
        else
            vals = log(abs(data(isel)))/log(10); %%%% log values
            oval = log(abs(data(1)))/log(10);
        end
        %oval = log(abs(data(1)))/log(10); %%%% optimal log par value
        %vals = data(isel);
   else
       vals = 0;
   end
   
   opar(iii)= data(1); %%% optimal values
   logvar(iii)= var(vals) ; %%% log variance
   
   %%%% histogram of logvalues %%%%%%%%%%%%%%%%%
   [eff,centers]=hist(vals,20);efftot = sum(eff); 

if data(1) > 0 & iii < 16 
    subplot(3,5,iii)
    hold off
    col ='r';
for j=1:length(centers)
    if eff(j) >0
    w=0.5;
    x=centers(j)-w/2;
    y=0;
    h=eff(j)/efftot;
    rectangle('Position',[x y w h],'FaceColor',col,'EdgeColor',col);
    end
end
axis([-8,1,0,1])
hold on
plot( [oval,oval],[0,1],'--','linewidth',2,'color',col) 
xlabel(['log_{10}(',parname,')']) 
end


end %%% for iii

print(hf,pfile,'-dpdf');


%if nstates==5  %%% 5 states models
%    ponhi=0.13;
%    ponno=0.006;
%    ponkinihi = - 1/(M.A1hi(1)/M.lambda1hi(1) + M.A1hi(2)/M.lambda1hi(2) + M.A1hi(3)/M.lambda1hi(3) + M.A1hi(4)/M.lambda1hi(4) + M.A1hi(5)/M.lambda1hi(5));
%    ponkinino =- 1/(M.A1no(1)/M.lambda1no(1) + M.A1no(2)/M.lambda1no(2) + M.A1no(3)/M.lambda1no(3) + M.A1no(4)/M.lambda1no(4) + M.A1no(5)/M.lambda1no(5));
%else
%    ponhi=0.13+0.29;
%    ponno=0.006+0.11; 
%    ponkinihi = - 1/(M.A1hi(1)/M.lambda1hi(1) + M.A1hi(2)/M.lambda1hi(2) + M.A1hi(3)/M.lambda1hi(3) + M.A1hi(4)/M.lambda1hi(4) );
%    ponkinino =- 1/(M.A1no(1)/M.lambda1no(1) + M.A1no(2)/M.lambda1no(2) + M.A1no(3)/M.lambda1no(3) + M.A1no(4)/M.lambda1no(4) );
%end
%%%% compute PON kini



switch jjj
    case 1,
ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k8no(1)); 
ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1)); 

k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1)];
ee=exp_fitnessMODEL1mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);  
 [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL1(k);
    case 2,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7no(1) + M.k9no(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7hi(1) + M.k9hi(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 

  k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL2mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
    [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL2(k);
    case 3,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7no(1) + M.k9no(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k3(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7hi(1) + M.k9hi(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k3(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 

 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL3mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962); 
    [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL3(k);
    case 4,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 
 
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL4mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);   
 [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL4(k);  
    case 5,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k3(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k3(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 

k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL5mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);     
    [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL5(k);
    case 6,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 

k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k9hi(1),M.k10hi(1)];

ee=exp_fitnessMODEL6mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
        [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL6(k);
    case 7,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k3(1)*M.k5(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k3(1)*M.k5(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 
 
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL7mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);     
    [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL7(k);
    case 8,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1)); 
 
 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL8mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
    [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL8(k);
    case 9,
 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k9no(1),M.k10no(1),...
M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL9mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL9(k);
    ponno = p4no; ponhi = p4hi;
     case 10,
 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k9no(1),M.k10no(1),...
M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL10mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL10(k);
    ponno = p4no; ponhi = p4hi;   
    case 11,
 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k9no(1),M.k10no(1),...
M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL11mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL11(k);
    ponno = p4no; ponhi = p4hi;
    case 12,
 k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k9no(1),M.k10no(1),M.k11no(1),M.k12no(1),...
M.k9hi(1),M.k10hi(1),M.k11hi(1),M.k12hi(1)];
ee=exp_fitnessMODEL12mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL12(k);
    ponno = p4no; ponhi = p4hi;    
    case 13,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7no(1) + M.k9no(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1)); 
 phoni= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7hi(1) + M.k9hi(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1));    
  k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];

ee=exp_fitnessMODEL13mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
   [p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL13(k); 
    case 14,
 ponno= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7no(1) + M.k9no(1) + M.k11no(1) + M.k12no(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k11no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k12no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k11no(1) + M.k1(1)*M.k3(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k9no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k8no(1) + M.k1(1)*M.k3(1)*M.k6no(1)*M.k12no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k11no(1) + M.k1(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k11no(1) + M.k1(1)*M.k3(1)*M.k9no(1)*M.k10no(1) + M.k1(1)*M.k4(1)*M.k6no(1)*M.k12no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k11no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k11no(1) + M.k2(1)*M.k4(1)*M.k7no(1)*M.k10no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k9no(1) + M.k1(1)*M.k3(1)*M.k8no(1)*M.k12no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k11no(1) + M.k1(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k1(1)*M.k5(1)*M.k6no(1)*M.k12no(1) + M.k2(1)*M.k4(1)*M.k6no(1)*M.k12no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k11no(1) + M.k1(1)*M.k3(1)*M.k10no(1)*M.k11no(1) + M.k1(1)*M.k4(1)*M.k8no(1)*M.k12no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k11no(1) + M.k2(1)*M.k4(1)*M.k9no(1)*M.k10no(1) + M.k2(1)*M.k5(1)*M.k6no(1)*M.k12no(1) + M.k1(1)*M.k3(1)*M.k10no(1)*M.k12no(1) + M.k1(1)*M.k4(1)*M.k10no(1)*M.k11no(1) + M.k2(1)*M.k4(1)*M.k8no(1)*M.k12no(1) + M.k3(1)*M.k5(1)*M.k6no(1)*M.k12no(1) + M.k1(1)*M.k4(1)*M.k10no(1)*M.k12no(1) + M.k2(1)*M.k4(1)*M.k10no(1)*M.k11no(1) + M.k2(1)*M.k4(1)*M.k10no(1)*M.k12no(1)); 
 ponhi= (M.k1(1)*M.k3(1)*M.k5(1)*(M.k7hi(1) + M.k9hi(1) + M.k11hi(1) + M.k12hi(1)))/(M.k1(1)*M.k3(1)*M.k5(1)*M.k6hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k7hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k11hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k5(1)*M.k12hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k11hi(1) + M.k1(1)*M.k3(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k9hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k8hi(1) + M.k1(1)*M.k3(1)*M.k6hi(1)*M.k12hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k11hi(1) + M.k1(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k11hi(1) + M.k1(1)*M.k3(1)*M.k9hi(1)*M.k10hi(1) + M.k1(1)*M.k4(1)*M.k6hi(1)*M.k12hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k11hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k11hi(1) + M.k2(1)*M.k4(1)*M.k7hi(1)*M.k10hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k9hi(1) + M.k1(1)*M.k3(1)*M.k8hi(1)*M.k12hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k11hi(1) + M.k1(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k1(1)*M.k5(1)*M.k6hi(1)*M.k12hi(1) + M.k2(1)*M.k4(1)*M.k6hi(1)*M.k12hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k11hi(1) + M.k1(1)*M.k3(1)*M.k10hi(1)*M.k11hi(1) + M.k1(1)*M.k4(1)*M.k8hi(1)*M.k12hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k11hi(1) + M.k2(1)*M.k4(1)*M.k9hi(1)*M.k10hi(1) + M.k2(1)*M.k5(1)*M.k6hi(1)*M.k12hi(1) + M.k1(1)*M.k3(1)*M.k10hi(1)*M.k12hi(1) + M.k1(1)*M.k4(1)*M.k10hi(1)*M.k11hi(1) + M.k2(1)*M.k4(1)*M.k8hi(1)*M.k12hi(1) + M.k3(1)*M.k5(1)*M.k6hi(1)*M.k12hi(1) + M.k1(1)*M.k4(1)*M.k10hi(1)*M.k12hi(1) + M.k2(1)*M.k4(1)*M.k10hi(1)*M.k11hi(1) + M.k2(1)*M.k4(1)*M.k10hi(1)*M.k12hi(1)); 

k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),M.k9no(1),M.k10no(1),M.k11no(1),M.k12no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1),M.k11hi(1),M.k12hi(1)];

ee=exp_fitnessMODEL14mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL14(k);
    case 15,
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL15mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL15(k);
    ponno = p4no; ponhi = p4hi;
    case 16,
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL16mergeoneponkini(k, 0.3, xno, sno, xhi, shi,  ppvalno4,ppvalhi4, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p1hi,p2hi,p3hi,p4hi]=proba4MODEL16(k);
    ponno = p4no; ponhi = p4hi;    
    case 17,
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k7no(1),M.k8no(1),M.k10no(1),...
M.k6hi(1),M.k7hi(1),M.k8hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL17mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL17(k);
ponno = p5no;
ponhi = p5hi;

    case 18,
k = [M.lambda1no(1),M.lambda2no(1),M.lambda3no(1),M.lambda4no(1),M.lambda5no(1),...
M.A1no(1),M.A2no(1),M.A3no(1),M.A4no(1),...
M.lambda1hi(1),M.lambda2hi(1),M.lambda3hi(1),M.lambda4hi(1),M.lambda5hi(1),...
M.A1hi(1),M.A2hi(1),M.A3hi(1),M.A4hi(1),...
M.k1(1),M.k2(1),M.k3(1),M.k4(1),M.k5(1),...
M.k6no(1),M.k8no(1),M.k9no(1),M.k10no(1),...
M.k6hi(1),M.k8hi(1),M.k9hi(1),M.k10hi(1)];
ee=exp_fitnessMODEL18mergeonekini_mergep4p5ponkinishort(k, 0.3, xno, sno, xhi, shi,  ppvalno,ppvalhi, 0.0022138, 0.23962);    
[p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi]=proba5MODEL18(k);
ponno = p5no;
ponhi = p5hi;
end


switch jjj
    case 1 
        kinihi= M.k8hi(1);   
        kinino= M.k8no(1);  
    otherwise
        kinihi= M.k10hi(1); 
        kinino= M.k10no(1); 
end
ponkinihi=kinihi*ponhi;
ponkinino=kinino*ponno;
%line = [ofmin,logvarobj,opar,logvar,kiniponno,kiniponhi,ponkinino,ponkinihi];

%line = [ofmin,logvarobj,opar,logvar,ponkinino,ponkinihi,p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi,ee',sum(ee.*ee)];

if nstates == 5
 line = [ofmin,logvarobj,opar,logvar,ponkinino,ponkinihi,p1no,p2no,p3no,p4no,p5no,p1hi,p2hi,p3hi,p4hi,p5hi,sum(ee.*ee)];
else
  line = [ofmin,logvarobj,opar,logvar,ponkinino,ponkinihi,p1no,p2no,p3no,p4no,0,p1hi,p2hi,p3hi,p4hi,0,sum(ee.*ee)];  
end



xlswrite(xfile,{namem},1,['A',num2str(iline)]);
xlswrite(xfile,line,1,['B',num2str(iline)]);
iline = iline + 1;


source = [directory,'figure_optimal.png'];
destination = strrep(pfile,'.pdf','_merge.png');
copyfile(source,destination)




end %%% for jjj


fsz = 12;
h=figure(1000)
hold off
x=-2.25:0.1:-0.9;
[h1,c1]=hist(log(dist1)/log(10),x);h1=h1/length(dist1);
[h2,c2]=hist(log(dist2)/log(10),x);h2=h2/length(dist2);
[h3,c3]=hist(log(dist3)/log(10),x);h3=h3/length(dist3);
[h4,c4]=hist(log(dist13)/log(10),x);h4=h4/length(dist13);
%[h5,c5]=hist(log(dist9)/log(10),x);
%[h6,c6]=hist(log(dist10)/log(10),x);
%[h7,c7]=hist(log(dist11)/log(10),x);
bar(x,[h1',h2',h3',h4'])
%bar(x,[h1',h2',h3',h4'])
xlabel('log_{10}(Obj)','fontsize',fsz)
ylabel('Solutions','fontsize',fsz)
legend({'Model1','Model2','Model3','Model13'},'Location','northwest')
set(gca,'fontsize',fsz)

print(h,'-dpdf','objective_values.pdf')




