%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% reads results of genetic algorithm for short movies, data from long
%%%%% movies
%%%%% needs 1) short movies decomvolution results  result_tat_name_short.mat
%%%%% 2) long movie data name_long.mat or long movie data name_long_raw.mat
%%%%% performs unconstrained two exponentials fit, compute parameters of
%%%%% the two states model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global alpha xs xl ss sl;
alpha = 0.6; %%% relative importance linear
thresh=1200; %%% threshold long movie    
%%% parameters
fsz=16; msz=5;lw=1;
outliers_long=1;
outliers_short=0;
visualize = 1; 
%fact=0.75; %%% [fact/(1-fact)]^2 is the relative importance of linear scale in the objective function
%alpha=0.9;  %%%% if alpha = 1 linear 
%%%%%%%%%%%%%%%%%%%%%%% Parameters 
TaillePreMarq = 700; % 700 bases (ref publi)
TailleSeqMarq = 2900; % 2900 bases (ref publi)
EspaceInterPolyMin = 10; % space between polymerase in bases 
Polym_speed = 67; % average speed bases per second (Ref publi)
TaillePostMarq = 1600+100*Polym_speed; % 1600 bases (ref publi)
FreqEchImg = (1/3); % 1/3 image per second data time sampling
Intensity_for_1_Polym = 1; %%%% 
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed; % signal length in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, August
%%%%% 2024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% reads MS2 data computes survival functions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% needs data from mat files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%result_tat_clean_OMX_42_H128_PST_WTc2cPosPred10.mat
%result_tat_clean_OMX_44_H128_PST_WTc2cPosPred10.mat
%result_tat_clean_OMX_48_WT_calib_RA_IntIntcPosPred10.mat
%result_tat_clean_OMX_60_WT_calib_RA_IntIntcPosPred10.mat
%result_tat_clean_OMX_77_noTatcPosPred10.mat
%result_tat_clean_OMX_80_noTatcPosPred10.mat
%data_tat_WT_long_raw.mat'
%data_tat_no_tat_long_raw.mat'


for iname=1:2

    
    switch iname
        case 1
lframes=300;            
name_short={'OMX_42_H128_PST_WTc2','OMX_44_H128_PST_WTc2','OMX_48_WT_calib_RA_IntInt','OMX_60_WT_calib_RA_IntInt'};
name_long='data_tat_WT_long_raw'; %%%% will use raw long data
cens = 0;
isel_vis=[4,30,41];
isel_vis_long=[2,4, 6];
%%%% no tat  
        case 2
lframes=400;            
name_short={'OMX_77_noTat','OMX_80_noTat'};
name_long='data_tat_no_tat_long_raw'; %%%% will use raw long data
cens = 1;
isel_vis=[28,52,95];
isel_vis_long=[1,2,4];
    end

%%%%% other parameters short movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstep = 3; %%% frames period in seconds
DureeSimu = lframes*tstep; %%% movie length in seconds
frame_num = lframes; %%% number of frames
DureeAnalysee = DureeSignal + DureeSimu ; % 
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many intervals (possible poly start position)
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed));     
tmax=DureeAnalysee; %%% duration of timeseries (consider all the polymerases between 0 and tmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%% load deconvolution result; concatenate if there are several files %%%%
Gexp=[]; Gpred=[]; Gpos=[]; Gfit = [];
for jjj=1:length(name_short)  %%% different movies
    fname = ['result_tat_clean_',name_short{jjj},'cPosPred10.mat'];
    load(fname)
    Gexp=[Gexp,DataExp];
    Gpred=[Gpred,DataPred];
    Gpos=[Gpos,cPosPred];
    Gfit=[Gfit;Fit];
end



DataExp = Gexp;
DataPred = Gpred;
cPosPred = Gpos;
Fit = Gfit;
%%%%%%%%%%%%%%%% load long movie data

fname = [name_long,'.mat'];
load(fname); %% loads DataExpLong and Time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd=size(DataPred);
n_cells=sd(2); %%%% number of cells


%%%%%%%%%%%%% eliminate outliers handling short movie %%%%%%%%%%%%%
if outliers_short
    Ev=zeros(n_cells,1); %%% based on quartile of the numbers of polymerases 
for i=1:n_cells
    Ev(i) = length(cPosPred{i});
end
Q = quantile(Ev,3);
isel=find(Q(1)-1.5*(Q(3)-Q(1)) < Ev & Ev < Q(3) + 1.5*(Q(3)-Q(1)));
DataExp=DataExp(:,isel);
DataPred=DataPred(:,isel);
cPosPred=cPosPred(isel);
Fit=Fit(isel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if outliers_short
n_cells=length(isel);
end 

%%%%% compute distribution of spacings short movie %%%%%%
dt=[]; dtc=[]; 
for i=1:n_cells
    indices = (cPosPred{i})';
    times = indices / FreqEchSimu;
    lt = length(times);
    switch lt
        case 0
            dtc=[dtc;tmax]; 
        case 1
            dtc=[dtc;tmax-times(1)];
        otherwise
        dtimes = diff(times); %%%% uncensured intervals
        %dtimes=dtimes(dtimes>0);
        dt=[dt;dtimes(1:end)];
        if tmax > times(end)
           dtc=[dtc;tmax-times(end)];
        end
        if times(1) > 0
            dt=[dt;times(1)];
        end
    end
end


dtg=[dt;dtc];
censored_short=[zeros(size(dt));ones(size(dtc))];




%%%%%%%%%%%%% outliers handling long movie %%%%%%%%%%%%%
if outliers_long
n = size(DataExpLong);
FI=zeros(n(2),1);
for i=1:n(2)
    N0=length(find(DataExpLong(:,i)==0));
    N1=length(find(DataExpLong(:,i)==1));
    FI(i) = N0/(N0+N1);
end


Q = quantile(FI,3);
isel=find(Q(1)-1.5*(Q(3)-Q(1)) < FI & FI < min([Q(3) + 1.5*(Q(3)-Q(1)),1]) & FI>0);
%isel=find(Q(1)-5*(Q(3)-Q(1)) < FI & FI < 1 & FI < Q(3) + 5*(Q(3)-Q(1)));
DataExpLong=DataExpLong(:,isel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
   
    
store = DataExpLong;
n = size(store);
%%%%% estimate distributions for long movies %%%%%%%%%%%%%%%%%%%%%
wt=[];wtc=[];
Ninactive=0; Total=0; tmaxlomax=0;
for i=1:n(2) %%% all cells
    ilast = find(store(:,i) ==2);
    if isempty(ilast)
        tmaxlo=Time(end)*60;
    else
        tmaxlo=Time(min(ilast))*60;
    end

    
    ind=  find(store(:,i)==1);
    
    if ~isempty(ind)
    laps = diff(ind);
    Ninactive = Ninactive+length(find(laps>1))+1;
    wtimes=(laps(laps>1)-1)*tstep*60;
    wlast=tmaxlo - ind(end)*tstep*60; %%%% tmaxlo depends on cell
    wt=[wt;wtimes]; %%% waiting times in seconds
    if wlast > 0
    wtc=[wtc;wlast]; %%% add last interval
    end
    else
       %wtc=[wtc;tmaxlo]; %%% exclude 
    end 
    Total = Total + tmaxlo;
    if tmaxlo > tmaxlomax
        tmaxlomax = tmaxlo;
    end
end

time=0:0.1:tmaxlomax; %%%% time in seconds

censored=[zeros(size(wt));ones(size(wtc))]; %%%% censored long movie



    
    
    
omin = 100;
for cens =0:1

    
%%%% short movie    
if cens
    [fsg,xsg,logg,upg]=ecdf(dtg,'censoring',censored_short);
else
    [fsg,xsg,logg,upg]=ecdf(dtg);
end


for nnn=0:4
    
    wwt=[wt;wtc]+3*(nnn)*60/2; %%%% shifted distribution from 0 t0 6 min 
        shift=3*nnn/2; %%% shift in minutes
 
 
%%%% estimate of integral for computing p1 
a=3*60;
Tinactive =sum(wwt);
Tactive = Total-Tinactive;
Nactive = Ninactive - n(2); %%%% number of active periods
imax=max(find(xsg<a));
imin=min(find(xsg>a));
x1=xsg(imin);f1=fsg(imin);
x2=xsg(imax);f2=fsg(imax);
f = f1 + (f2-f1)*(a-x1)/(x2-x1);
esp=(-a*(1-f)+trapz(xsg(xsg<a),1-fsg(xsg<a)))/f;
estp1= Ninactive/(Ninactive + Tactive/esp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 

%%%%% long movie
if cens
    [fl,xl]=ecdf(wwt,'censoring',censored);
else
    [fl,xl]=ecdf(wwt);
end


%%%% fit p2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fact1=1;fact2=1;
[uxsg,iu]=unique(xsg);
ufsg=fsg(iu);
[uxl,iu]=unique(xl);
ufl=fl(iu);
M=max(uxsg)/fact2;
m=min(uxl(2:end))*fact1;

if M > m
%%%%% interpolate on [m,M];
nsteps=100;
step = (M-m)/nsteps;
x=m:step:M;
y1=interp1(uxsg,ufsg,x,'cubic'); %%% fast 
y2=interp1(uxl,ufl,x,'cubic'); %%% slow
nn=length(y1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=estp1;
opts = optimoptions(@lsqnonlin,'TolFun', 1e-15,'MaxIter',1e4, ...
        'MaxFunEvals',1e6,'TolX', 1e-10);
objective = @(k)(  (log((1-y2)*p1) - log( (1-y1)*(1-k(1)) + k(1) ))/nn );
     k0=0.001;    
     [k, obj] = lsqnonlin(objective,k0,0,1,opts);
p2=k;
if obj < omin
    omin = obj;
    p2min = p2;
    p1min = p1;
    xsgmin = xsg;
    fsgmin = fsg;
    xlmin  = xl; 
    flmin  = fl;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visualize
h=figure(200+nnn);
hold off
loglog(xsg,(1-fsg)*(1-p2)+p2,'og')
hold on
loglog(xl,(1-fl)*p1,'ok')
legend({'short movie','long movie'})
%loglog(time,s,'k')
axis([1e-1, 1e4, 1e-6, 1e0])
xlabel('Time [s]','fontsize',fsz)
ylabel('Frequency','fontsize',fsz)
set(gca,'fontsize',fsz)
title(['\Delta_0 =',num2str(3*nnn)],'fontsize',fsz)
end% if visualize

%%%% fit distribution of spacings using combination of 3 exponentials. 5 
%%%% params
xl=xl(1:end-1);
fl=fl(1:end-1);
NS = length(xsg); NL = length(xl);
sNS=sqrt(NS);sNL=sqrt(NL); 
fact1=sqrt(1-alpha);
fact2=sqrt(alpha);
 
eshort=(1-fsg)*(1-p2)+p2;
elong = (1-fl)*p1;
 
end %%%% if m < M    
end %%% for nnn

end %%% cens



xs = xsgmin;
ss = (1-fsgmin)*(1-p2min)+p2min;
xl = xlmin;
sl = (1-flmin)*p1min;

switch iname
    case 1,
        xshi=xs; 
        xlhi=xl; 
        sshi=ss;
        slhi=sl;
        xlhi = xlhi(1:end-1);
        slhi=  slhi(1:end-1);
    case 2,
        xsno=xs; 
        xlno=xl; 
        ssno=ss;
        slno=sl;
end
%xl = xl(1:end-1);
%sl=  sl(1:end-1);



end %%% iname

save('s_functions_sl.mat','xshi','xsno','xlhi','xlno','sshi','ssno','slhi','slno');

