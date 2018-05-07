function [gt , Q_Time] = Simu_CumTime(m,sigma,nGoats,F)

%% Simulate time in pasture for each goat
gt = normrnd(m,sigma,[1 , nGoats]);

%% Compute the cumulative time in each quadrat
Q_Time = F.*sum(gt);