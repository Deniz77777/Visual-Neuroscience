%% HW2 - Deniz Rezapour
% part A 
% fig 2
clc
close all
clear all

clearvars -except Dacc DRTA DRT1 DRT2
%%%% Parameters to be varied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.005;           % Time step in msec
thresh    = 15;        % Decision threshold     
% Coherence level
coh       = [0, 51.2];
mu0       = 30.0;      % External stimulus strength
w=1.7;

MU0=10:5:40;
COH = logspace(0, 2, 10);
W = [1.68, 1.70, 1.72];

[r1_traj r2_traj]= wang(coh(1),mu0,w); % for %coh 0.0
[r11_traj r22_traj]= wang(coh(2),mu0,w); % for %coh 51.2

figure()
plot(dt:dt:dt*size(r1_traj,2),r1_traj','b','linewidth',0.5)
hold on
plot(dt:dt:dt*size(r2_traj,2),r2_traj','r','linewidth',0.5)
hold on
plot(dt:dt:dt*size(r11_traj,2),r11_traj','m','linewidth',1)
hold on
plot(dt:dt:dt*size(r22_traj,2),r22_traj','k','linewidth',1)
yline(thresh,'--k','linewidth',2.5)
title('Firing Rate vs. Time')
ylabel('Firing Rate')
xlabel('Time[s]')

figure()
plot(dt:dt:dt*size(r1_traj,2),r1_traj','b')
hold on
plot(dt:dt:dt*size(r2_traj,2),r2_traj','r')
yline(thresh,'--k','linewidth',2)
title('Firing Rate vs. Time in %coh = 0%')
ylabel('Firing Rate')
xlabel('Time[s]')

figure()
plot(dt:dt:dt*size(r11_traj,2),r11_traj','m')
hold on
plot(dt:dt:dt*size(r22_traj,2),r22_traj','k')
yline(thresh,'--k','linewidth',2)
title('Firing Rate vs. Time in %coh = 51.2%')
ylabel('Firing Rate')
xlabel('Time[s]')
% fig 6-A
% %coh is 0
close all
finalRTs = zeros(3, 7);
maxRTs = zeros(3, 7);
minRTs = zeros(3, 7);
for i=1:7
    [r1 r2]= wang(coh(1),MU0(i),w); 
    [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
    finalRTs(1, i) = RT_all;
    finalRTs(2, i) = RT_pref;
    finalRTs(3, i) = RT_non;
    maxRTs(1, i) = max(RT_allArray);
    maxRTs(2, i) = max(RT_prefArray);
    maxRTs(3, i) = max(RT_nonArray);
    minRTs(1, i) = min(RT_allArray);
    minRTs(2, i) = min(RT_prefArray);
    minRTs(3, i) = min(RT_nonArray);
end

figure()
plot(MU0, finalRTs(2,:), '-o')
hold on
for j=1:7
    line([MU0(j) MU0(j)],  [maxRTs(2, j) minRTs(2, j)])
end
title('Reaction Time of Correct Decision Trials vs. mu0')
ylabel('RT[s]')
xlabel('mu0[Hz]')
xlim([5 45])
% fig 11
RTs11 = zeros(3, length(COH));
acc11 = zeros(3, length(COH));
for k=1:length(COH)
    for kk=1:3
        [r1 r2]= wang(COH(k),mu0,W(kk)); 
        [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
        RTs11(kk, k) = RT_pref;
        acc11(kk, k) = acc;
    end
end

figure()
semilogx(COH, RTs11(1,:), 'b');
hold on
semilogx(COH, RTs11(2,:), 'r');
hold on
semilogx(COH, RTs11(3,:), 'k');
legend('w+ = 1.68', 'w+ = 1.70', 'w+ = 1.72')
title('Reaction Time vs. %Coh')
xlabel('%coh (in log)')
ylabel('RT[s]')


fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(COH)],...
               'StartPoint',[0 0.5]);
ft = fittype('1-0.5*exp(-(x/b)^a)','options',fo);
[curve1,~] = fit(COH',acc11(1,:)',ft)
[curve2,~] = fit(COH',acc11(2,:)',ft)
[curve3,~] = fit(COH',acc11(3,:)',ft)
figure()
x = [0:0.01:max(COH)];
plot(x,1-0.5*exp(-(x./curve1.b).^curve1.a),'b','linewidth',1.5);
hold on
plot(x,1-0.5*exp(-(x./curve2.b).^curve2.a),'r','linewidth',1.5);
hold on
plot(x,1-0.5*exp(-(x./curve3.b).^curve3.a),'k','linewidth',1.5);
set(gca,'XScale','log')
legend('w+ = 1.68', 'w+ = 1.70', 'w+ = 1.72')
title('Accuracy vs. %Coh')
xlabel('%coh (in log)')
ylabel('Acc')

figure()
semilogx(COH, acc11(1,:), 'b');
hold on
semilogx(COH, acc11(2,:), 'r');
hold on
semilogx(COH, acc11(3,:), 'k');
legend('w+ = 1.68', 'w+ = 1.70', 'w+ = 1.72')
title('Accuracy vs. %Coh (not fitted with a function)')
xlabel('%coh (in log)')
ylabel('Acc')

%% Extra part (effect of threshold on performance)
close all
clc
ths = [5, 15, 25, 35];
RTs11 = zeros(4, length(COH));
acc11 = zeros(4, length(COH));
for k=1:length(COH)
    for kk=1:4
        [r1 r2]= wang_sup(COH(k),mu0,1.7, ths(kk)); 
        [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
        RTs11(kk, k) = RT_pref;
        acc11(kk, k) = acc;
    end
end

figure()
semilogx(COH, RTs11(1,:), 'b');
hold on
semilogx(COH, RTs11(2,:), 'r');
hold on
semilogx(COH, RTs11(3,:), 'k');
hold on
semilogx(COH, RTs11(4,:), 'm');
legend('threshold = 5', 'threshold = 15', 'threshold = 25', 'threshold = 35')
title('Reaction Time vs. %Coh')
xlabel('%coh (in log)')
ylabel('RT[s]')


fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(COH)],...
               'StartPoint',[0 0.5]);
ft = fittype('1-0.5*exp(-(x/b)^a)','options',fo);
[curve1,~] = fit(COH',acc11(1,:)',ft)
[curve2,~] = fit(COH',acc11(2,:)',ft)
[curve3,~] = fit(COH',acc11(3,:)',ft)
[curve4,~] = fit(COH',acc11(4,:)',ft)
figure()
x = [0:0.01:max(COH)];
plot(x,1-0.5*exp(-(x./curve1.b).^curve1.a),'b','linewidth',1.5);
hold on
plot(x,1-0.5*exp(-(x./curve2.b).^curve2.a),'r','linewidth',1.5);
hold on
plot(x,1-0.5*exp(-(x./curve3.b).^curve3.a),'k','linewidth',1.5);
hold on
plot(x,1-0.5*exp(-(x./curve4.b).^curve4.a),'m','linewidth',1.5);
set(gca,'XScale','log')
legend('threshold = 5', 'threshold = 15', 'threshold = 25', 'threshold = 35')
title('Accuracy vs. %Coh')
xlabel('%coh (in log)')
ylabel('Acc')

% RT vs threshold
thss = [5, 10, 15, 20, 25, 30, 35, 40];
finalRTs = zeros(3, 7);
for i=1:8
    [r1 r2]= wang_sup(coh(1),mu0,w, thss(i)); 
    [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
    finalRTs(1, i) = RT_all;
    finalRTs(2, i) = RT_pref;
    finalRTs(3, i) = RT_non;
end
figure()
plot(thss, finalRTs(2,:), '-o')
title('Reaction Time of Correct Decision Trials vs. Threshold')
ylabel('RT[s]')
xlabel('Threshold[Hz]')
xlim([0 45])
ylim([0.2 0.8])
%% part B
clc
clear all
dt = 0.005;           % Time step in msec


load('subData.mat');
data = table2array(subData);
prefered = find(data(:,3)==data(:,4));
data(prefered,5) = 1;
COHs = [1.6, 3.2, 6.4, 12.8, 25.6];

for c = 1:length(COHs)
    coh = COHs(c);
    cohtemp = find(data(:,1) == coh);
    corrtemp = find(data(:,5) == 1);
    incorrtemp = find(data(:,5)==0);
    tempCorr = intersect(cohtemp,corrtemp);
    tempInCorr = intersect(cohtemp,incorrtemp);
    RT.coh{c,1} = data(tempCorr,2);
    MeanRTCorrEx(c) = mean(RT.coh{c,1})/1000;
    StdRTCorrEx(c) = std(RT.coh{c,1})/1000;
    RT.coh{c,2} = data(tempInCorr,2);
    MeanRTInCorrEx(c) = mean(RT.coh{c,2})/1000;
    StdRTInCorrEx(c) = std(RT.coh{c,2})/1000;
    acc_ex(c) = length(tempCorr)/(length(tempCorr)+length(tempInCorr));
end

MU0=[31, 33, 35, 38, 40, 36.0011];
W = [1.5, 1.60, 1.65, 1.75, 1.80, 1.6779];

RTsPref = zeros(100, length(COHs));
RTsNonPref = zeros(100, length(COHs));
ACCs = zeros(100, length(COHs));
for i=1:100
    for j=1:length(COHs)
        coh = COHs(j);
        [r1 r2]= wang(coh,MU0(6),W(6)); 
        [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
        RTsPref(i, j) = RT_pref;
        ACCs(i, j) = acc;
    end
end
meanRTsPref = mean(RTsPref);
meanRTsNonPref = mean(RTsNonPref);
meanACCs = mean(ACCs,1);

acc_EX =  [0.6000,0.6617,0.8283,0.9383,0.9967];
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(COHs)],...
               'StartPoint',[0.1 0.1]);
ft = fittype('1-0.5*exp(-(x/b)^a)','options',fo);
[curve1,~] = fit(COHs',meanACCs',ft)
[curve2,~] = fit(COHs',acc_EX',ft)

figure()
plot(COHs, meanRTsPref, '-o')
hold on
scatter(COHs,MeanRTCorrEx,'r','filled');
hold on
plot(COHs,MeanRTCorrEx,'r')
set(gca,'XScale','log')
h=get(gca,'Children')
legend(h([3 2]),'Model','Experiment')
title('Reaction Time vs. %Coh in Model and Experiment')
ylabel('RT[s]')
xlabel('%coh')
ylim([0 1])

figure()
x = [0.1:0.01:max(COHs)];
scatter(COHs,meanACCs,'b','filled');
hold on
plot(x,1-0.5*exp(-(x./curve1.b).^curve1.a),'b','linewidth',1.5);
hold on
scatter(COHs,acc_EX,'r','filled')
hold on
plot(x,1-0.5*exp(-(x./curve2.b).^curve2.a),'r','linewidth',1.5);
set(gca,'XScale','log')
h=get(gca,'Children');
legend(h([4 2]),'Model','Experiment')
title('Accuracy vs. %Coh in Model and Experiment')
ylabel('ACC')
xlabel('%coh')

%% part C
clc
close all
clear all

W = 1.38;
dt = 0.005;           % Time step in msec

cohVec = [1.6, 3.2, 6.4, 12.8, 25.6];
mu0 = 37;     

fun= @(p) final_optimization(p(1),p(2),cohVec)
objFun = @(p) fun(p);
sol= ga(objFun,2,[],[],[],[],[1.5 25],[1.9 40],[])


function [y]=final_optimization(W,mu0,COHs)
RTsPref = zeros(1, length(COHs));
RTsNonPref = zeros(1, length(COHs));
ACCs = zeros(100, length(COHs));
dt = 0.005;           % Time step in msec
W
mu0
for j=1:length(COHs)
    coh = COHs(j);
    [r1 r2]= wang(coh,mu0,W); 
    [RT_all,RT_pref RT_non acc RT_prefArray RT_nonArray RT_allArray]=RT_calc(r1,r2,dt);
    RTsPref(1, j) = RT_pref;
    ACCs(1, j) = acc;
end
meanRTsPref = mean(RTsPref);
meanRTsNonPref = mean(RTsNonPref);
meanACCs = mean(ACCs,1);
acc_ex =  [0.6000,0.6617,0.8283,0.9383,0.9967];

RTs =  [0.658077777777778,0.649216624685139,0.630187122736419,0.577255772646537,0.510041806020067];
y=norm(RTs-meanRTsPref)+norm(acc_ex-meanACCs);
end

%% function
function [RTA,RT1 RT2 acc rt1 rt2 rtA ]=RT_calc(r1_traj,r2_traj,dt)

r1=r1_traj(:,101:400);
r2=r2_traj(:,101:400);

for ii=1:100
    temp1=find(r1(ii,:)>=15);
    temp2=find(r2(ii,:)>=15);
    
    if temp1~[];
        rt1(ii)=temp1(1);
        rtA(ii)=temp1(1);
    else
        rt1(ii)=0;
    end
    
    if temp2~[];
        rt2(ii)=temp2(1);
        rtA(ii)=temp2(1);
    else
        rt2(ii)=0;
    end
end
rt1=rt1*dt + 0.35;
rt2=rt2*dt + 0.35;
rtA=rtA*dt + 0.35;

RT1=mean(rt1);
RT2=mean(rt2);
RTA=mean(rtA);

acc=sum(rt2==0.35)/100;

end

function [r1_traj r2_traj]= wang(coh,mu0,jn)

%%%% FI curve parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 270; b = 108;
d = 0.1540;  % Parameters for excitatory cells

%%%% Vectorizing variables for evaluations at end of (block) loop %%%%%%%%%%

r1_traj = [];  r2_traj = [];
s1_traj = [];  s2_traj = [];
cross_r1 = []; cross_r2 = [];

%--------------------------------------------------------------------------
Tnmda = 100;    % NMDAr
Tampa = 2;      % AMPAr
gamma = 0.641;

%--------------------------------------------------------------------------
thresh    = 15;        % Decision threshold
thresh_s  = (gamma*Tnmda*thresh/1000)/(1+gamma*Tnmda*thresh/1000); % Threshold in s-space
noise_amp = 0.020;      % Noise amplitude into selective populations
N_trials  = 100 ;       % Total number of trials
Tstim     = 2500;      % Stimulus duration (ms)

%%%%% Stimulus input strengths %%%%%

mu1 = mu0*(1+coh/100); % input strength to pop 1 (for coherence)
mu2 = mu0*(1-coh/100); % input strength to pop 2 (for coherence)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Trials number and (block) loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ww = 1:N_trials % Trial loop
    
    trial_no = ww;
    
    %---- Vectorise variables with sliding window -------------
    nu1_wind = [] ; nu2_wind = [] ;
    s1_wind  = [] ; s2_wind  = [] ;
    
    %---- Initial conditions and clearing variables -----------
    s1_in=0.1; s2_in=0.1;
    clear nu1_in nu2_in I_eta1 I_eta2  ;
    nu1_in = 2; nu2_in = 2;
    I_eta1_in = noise_amp*randn ; I_eta2_in = noise_amp*randn ;
    
    %---- Time conditions -------------------------------------
    
    dt = 0.5;           % Time step in msec
    T_total = 3000/dt;  % Total number of steps
    time_wind = 50/dt;  % Temporal window size for averaging
    slide_wind = 5/dt;  % Sliding step for window
    
    %---- Intialise and vectorise variables to be used in loops below ------
    
    s1 = s1_in.*ones(1,T_total); s2 = s2_in.*ones(1,T_total);
    nu1 = nu1_in.*ones(1,T_total); nu2 = nu2_in.*ones(1,T_total);
    phi1 = nu1_in.*ones(1,T_total); phi2 = nu2_in.*ones(1,T_total);
    I_eta1 = I_eta1_in.*ones(1,T_total); I_eta2 = I_eta2_in.*ones(1,T_total);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for t = 1:T_total
        
        %---- Constant effective external current input (with inhibition taken into account)
        I0E1 = 0.3255; I0E2 = 0.3255;
        
        %---- External stimulus---------------------------------------------------------
        JAext = 0.00052; % Synaptic coupling constant to external inputs
        I_stim_1 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu1); % To population 1
        I_stim_2 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu2); % To population 2
        
        %---- Recurrent synaptic coupling constants-------------------------------------
        %JN11 = 0.2609; JN22 = 0.2609;
        JN11 = 0.1631*jn; JN22 = 0.1631*jn;
        JN12 = 0.057; JN21 = 0.057; 
        
        %---- Resonse function of competiting excitatory population 1 ------
        Isyn1(t) = JN11.*s1(t) - JN12.*s2(t) + I0E1 + I_stim_1 + I_eta1(t) ;
        phi1(t)  = (a.*Isyn1(t)-b)./(1-exp(-d.*(a.*Isyn1(t)-b)));
        
        %---- Response function of competiting excitatory population 2 -----
        Isyn2(t) = JN22.*s2(t) - JN21.*s1(t) + I0E2 + I_stim_2 + I_eta2(t) ;
        phi2(t)  = (a.*Isyn2(t)-b)./(1-exp(-d.*(a.*Isyn2(t)-b)));
        
        %---- Dynamical equations -------------------------------------------
        
        % Mean NMDA-receptor dynamics
        s1(t+1) = s1(t) + dt*(-(s1(t)/Tnmda) + (1-s1(t))*gamma*nu1(t)/1000);
        s2(t+1) = s2(t) + dt*(-(s2(t)/Tnmda) + (1-s2(t))*gamma*nu2(t)/1000);
        
        % Noise through synaptic currents of pop1 and 2
        I_eta1(t+1) = I_eta1(t) + (dt/Tampa)*(-I_eta1(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        I_eta2(t+1) = I_eta2(t) + (dt/Tampa)*(-I_eta2(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        
        % To ensure firing rates are always positive (noise may cause negative)
        if phi1(t) < 0
            nu1(t+1) = 0;
            phi1(t) = 0;
        else
            nu1(t+1) = phi1(t);
        end;
        if phi2(t) < 0
            nu2(t+1) = 0;
            phi2(t) = 0;
        else
            nu2(t+1) = phi2(t);
        end;
        
        %==============================================================================================
        
    end;  %---- End of time loop --------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %---- Calculating the mean rates and gating variables with sliding window -----
    
    nu1_wind = [nu1_wind (mean(nu1(1:time_wind)))] ;
    nu2_wind = [nu2_wind (mean(nu2(1:time_wind)))] ;
    s1_wind  = [s1_wind (mean(s1(1:time_wind)))] ;
    s2_wind  = [s2_wind (mean(s2(1:time_wind)))] ;
    
    for t = 1:(T_total-time_wind)/slide_wind
        
        nu1_wind = [nu1_wind (mean(nu1(slide_wind*t:slide_wind*t+time_wind)))] ;
        nu2_wind = [nu2_wind (mean(nu2(slide_wind*t:slide_wind*t+time_wind)))] ;
        s1_wind  = [s1_wind (mean(s1(slide_wind*t:slide_wind*t+time_wind)))] ;
        s2_wind  = [s2_wind (mean(s2(slide_wind*t:slide_wind*t+time_wind)))] ;
        
    end;
    
    r1_traj=[r1_traj; nu1_wind]; r2_traj=[r2_traj; nu2_wind];
    s1_traj=[s1_traj; s1_wind]; s2_traj=[s2_traj; s2_wind];
    clear nu1 nu2 s1 s2;
    clear nu1_wind nu2_wind s1_wind s2_wind ;
    
end; %---- End trial loop ---------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Plots of first few trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N_traj = 10; % Total number of trajectories
%
% subplot(2,2,1:2) % Plot timecourse of activity (population firing rates)
% for ww = 1:N_traj
%     plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r1_traj(ww,1:end-10),'b');hold on;
%     plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r2_traj(ww,1:end-10),'r--');hold on;
% end;
% %grid on;
% %axis tight;
% hold on; plot([1 dt*(T_total-time_wind)],[thresh thresh],'k--');
% xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
% title('Timecourse of firing rates r_1 (blue) and r_2 (red)');
%
% subplot(2,2,3) % Plot trajectories in phase-space
% for ww=1:N_traj
%     if s1_traj(ww,end)>thresh_s
%         plot(s2_traj(ww,:),s1_traj(ww,:),'b');hold on;
%     else
%         plot(s2_traj(ww,:),s1_traj(ww,:),'r');hold on;
%     end;
% end;
% hold on; plot([0 1],[thresh_s thresh_s],'k--');
% hold on; plot([thresh_s thresh_s],[0 1],'k--');
% xlabel('S_2'); ylabel('S_1');
%
% subplot(2,2,4) % Plot trajectories in phase-space
% for ww=1:N_traj
%     if r1_traj(ww,end)>thresh
%         plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'b');hold on;
%     else
%         plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'r');hold on;
%     end;
% end;
% hold on; plot([0 45],[thresh thresh],'k--');
% hold on; plot([thresh thresh],[0 45],'k--');
% xlabel('r_2 (Hz)'); ylabel('r_1 (Hz)');
% axis tight;




end

function [r1_traj r2_traj]= wang_sup(coh,mu0,jn, threshold)

%%%% FI curve parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 270; b = 108;
d = 0.1540;  % Parameters for excitatory cells

%%%% Vectorizing variables for evaluations at end of (block) loop %%%%%%%%%%

r1_traj = [];  r2_traj = [];
s1_traj = [];  s2_traj = [];
cross_r1 = []; cross_r2 = [];

%--------------------------------------------------------------------------
Tnmda = 100;    % NMDAr
Tampa = 2;      % AMPAr
gamma = 0.641;

%--------------------------------------------------------------------------
thresh    = threshold;        % Decision threshold
thresh_s  = (gamma*Tnmda*thresh/1000)/(1+gamma*Tnmda*thresh/1000); % Threshold in s-space
noise_amp = 0.020;      % Noise amplitude into selective populations
N_trials  = 100 ;       % Total number of trials
Tstim     = 2500;      % Stimulus duration (ms)

%%%%% Stimulus input strengths %%%%%

mu1 = mu0*(1+coh/100); % input strength to pop 1 (for coherence)
mu2 = mu0*(1-coh/100); % input strength to pop 2 (for coherence)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Trials number and (block) loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ww = 1:N_trials % Trial loop
    
    trial_no = ww;
    
    %---- Vectorise variables with sliding window -------------
    nu1_wind = [] ; nu2_wind = [] ;
    s1_wind  = [] ; s2_wind  = [] ;
    
    %---- Initial conditions and clearing variables -----------
    s1_in=0.1; s2_in=0.1;
    clear nu1_in nu2_in I_eta1 I_eta2  ;
    nu1_in = 2; nu2_in = 2;
    I_eta1_in = noise_amp*randn ; I_eta2_in = noise_amp*randn ;
    
    %---- Time conditions -------------------------------------
    
    dt = 0.5;           % Time step in msec
    T_total = 3000/dt;  % Total number of steps
    time_wind = 50/dt;  % Temporal window size for averaging
    slide_wind = 5/dt;  % Sliding step for window
    
    %---- Intialise and vectorise variables to be used in loops below ------
    
    s1 = s1_in.*ones(1,T_total); s2 = s2_in.*ones(1,T_total);
    nu1 = nu1_in.*ones(1,T_total); nu2 = nu2_in.*ones(1,T_total);
    phi1 = nu1_in.*ones(1,T_total); phi2 = nu2_in.*ones(1,T_total);
    I_eta1 = I_eta1_in.*ones(1,T_total); I_eta2 = I_eta2_in.*ones(1,T_total);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for t = 1:T_total
        
        %---- Constant effective external current input (with inhibition taken into account)
        I0E1 = 0.3255; I0E2 = 0.3255;
        
        %---- External stimulus---------------------------------------------------------
        JAext = 0.00052; % Synaptic coupling constant to external inputs
        I_stim_1 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu1); % To population 1
        I_stim_2 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu2); % To population 2
        
        %---- Recurrent synaptic coupling constants-------------------------------------
        %JN11 = 0.2609; JN22 = 0.2609;
        JN11 = 0.1631*jn; JN22 = 0.1631*jn;
        JN12 = 0.057; JN21 = 0.057; 
        
        %---- Resonse function of competiting excitatory population 1 ------
        Isyn1(t) = JN11.*s1(t) - JN12.*s2(t) + I0E1 + I_stim_1 + I_eta1(t) ;
        phi1(t)  = (a.*Isyn1(t)-b)./(1-exp(-d.*(a.*Isyn1(t)-b)));
        
        %---- Response function of competiting excitatory population 2 -----
        Isyn2(t) = JN22.*s2(t) - JN21.*s1(t) + I0E2 + I_stim_2 + I_eta2(t) ;
        phi2(t)  = (a.*Isyn2(t)-b)./(1-exp(-d.*(a.*Isyn2(t)-b)));
        
        %---- Dynamical equations -------------------------------------------
        
        % Mean NMDA-receptor dynamics
        s1(t+1) = s1(t) + dt*(-(s1(t)/Tnmda) + (1-s1(t))*gamma*nu1(t)/1000);
        s2(t+1) = s2(t) + dt*(-(s2(t)/Tnmda) + (1-s2(t))*gamma*nu2(t)/1000);
        
        % Noise through synaptic currents of pop1 and 2
        I_eta1(t+1) = I_eta1(t) + (dt/Tampa)*(-I_eta1(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        I_eta2(t+1) = I_eta2(t) + (dt/Tampa)*(-I_eta2(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        
        % To ensure firing rates are always positive (noise may cause negative)
        if phi1(t) < 0
            nu1(t+1) = 0;
            phi1(t) = 0;
        else
            nu1(t+1) = phi1(t);
        end;
        if phi2(t) < 0
            nu2(t+1) = 0;
            phi2(t) = 0;
        else
            nu2(t+1) = phi2(t);
        end;
        
        %==============================================================================================
        
    end;  %---- End of time loop --------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %---- Calculating the mean rates and gating variables with sliding window -----
    
    nu1_wind = [nu1_wind (mean(nu1(1:time_wind)))] ;
    nu2_wind = [nu2_wind (mean(nu2(1:time_wind)))] ;
    s1_wind  = [s1_wind (mean(s1(1:time_wind)))] ;
    s2_wind  = [s2_wind (mean(s2(1:time_wind)))] ;
    
    for t = 1:(T_total-time_wind)/slide_wind
        
        nu1_wind = [nu1_wind (mean(nu1(slide_wind*t:slide_wind*t+time_wind)))] ;
        nu2_wind = [nu2_wind (mean(nu2(slide_wind*t:slide_wind*t+time_wind)))] ;
        s1_wind  = [s1_wind (mean(s1(slide_wind*t:slide_wind*t+time_wind)))] ;
        s2_wind  = [s2_wind (mean(s2(slide_wind*t:slide_wind*t+time_wind)))] ;
        
    end;
    
    r1_traj=[r1_traj; nu1_wind]; r2_traj=[r2_traj; nu2_wind];
    s1_traj=[s1_traj; s1_wind]; s2_traj=[s2_traj; s2_wind];
    clear nu1 nu2 s1 s2;
    clear nu1_wind nu2_wind s1_wind s2_wind ;
    
end; %---- End trial loop ---------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Plots of first few trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N_traj = 10; % Total number of trajectories
%
% subplot(2,2,1:2) % Plot timecourse of activity (population firing rates)
% for ww = 1:N_traj
%     plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r1_traj(ww,1:end-10),'b');hold on;
%     plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r2_traj(ww,1:end-10),'r--');hold on;
% end;
% %grid on;
% %axis tight;
% hold on; plot([1 dt*(T_total-time_wind)],[thresh thresh],'k--');
% xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
% title('Timecourse of firing rates r_1 (blue) and r_2 (red)');
%
% subplot(2,2,3) % Plot trajectories in phase-space
% for ww=1:N_traj
%     if s1_traj(ww,end)>thresh_s
%         plot(s2_traj(ww,:),s1_traj(ww,:),'b');hold on;
%     else
%         plot(s2_traj(ww,:),s1_traj(ww,:),'r');hold on;
%     end;
% end;
% hold on; plot([0 1],[thresh_s thresh_s],'k--');
% hold on; plot([thresh_s thresh_s],[0 1],'k--');
% xlabel('S_2'); ylabel('S_1');
%
% subplot(2,2,4) % Plot trajectories in phase-space
% for ww=1:N_traj
%     if r1_traj(ww,end)>thresh
%         plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'b');hold on;
%     else
%         plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'r');hold on;
%     end;
% end;
% hold on; plot([0 45],[thresh thresh],'k--');
% hold on; plot([thresh thresh],[0 45],'k--');
% xlabel('r_2 (Hz)'); ylabel('r_1 (Hz)');
% axis tight;




end
