function Risk = Simu_Trajectory(nq,F_W1,F_W2,Week_Time,nGoats,mu_wf,sigma_wf,mu_np,sigma_np,Evap,Precipitation,Temperature,T_Zone,ZN,OPG_Goats,W,Spec)

%%%% WEEK 1 %%%%

%% Compute specificity 
% [W , Spec] = Simu_Specificity(nq,nGoats);

%% Simulate time in pasture for each goat and cumulative time in each quadrat for the first week
[gt , Q_Time] = Simu_CumTime(Week_Time/nGoats,0.1*(Week_Time/nGoats),nGoats,F_W1);

%% Simulate the time spend per quadrat q per goat g G_Time(g,q) for the first week
G_Time = Simu_TimeQuadratPerGoat(nq,W,gt,Q_Time,nGoats);

%% Simulate the OPG per goat
% OPG_Goats = Simu_OPG(OPG,nGoats);

%% Simulate spatial distribution of eggs
Polution_t = Simu_FecesDrop(nGoats,nq,mu_wf,sigma_wf,mu_np,sigma_np,G_Time,gt,OPG_Goats);

%% Simulate eggs development to L3
Polution_tp1 = Simu_Larva(Polution_t,Temperature,Evap,Precipitation);



%%%% WEEK 2 %%%%
%% Simulate time in pasture for each goat and cumulative time in each quadrat for the second week
[gt , Q_Time] = Simu_CumTime(Week_Time/nGoats,0.1*(Week_Time/nGoats),nGoats,F_W2);

%% Simulate the time spend per quadrat q per goat g G_Time(g,q) for the second week
G_Time = Simu_TimeQuadratPerGoat(nq,W,gt,Q_Time,nGoats);

%% Simulate risk
Risk = Simu_Risk(nGoats,nq,Polution_tp1,G_Time,T_Zone,ZN);