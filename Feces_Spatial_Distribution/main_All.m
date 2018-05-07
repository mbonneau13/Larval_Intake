%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G1
%% Load all the necesseray data
clear
load Data.mat
load Meteo.mat
load ZNumber.mat
load Distri_Weights_Feces.mat

%% Fix the number of goats in the flock and flock number (Groupe)
nGoats = 21;
Groupe = 1;

%% The field is divided into Zone, depending on the dominent weed species. 
if Groupe == 1
    ZN = ZN_G1;
else
    ZN = ZN_G2;
end

%% Compute the occurrence frequency for the first week
H = ComputeProportion(Groupe,Data,[1 4]);
Idx = find(H > 0);
F_W1 = H(Idx);
nq = length(F_W1);

%% Compute the occurrence frequency for the second week
H = ComputeProportion(Groupe,Data,[5 8]);
F_W2 = H(Idx);


%% Compute the average time spend in the field
T = ComputeTime(Data,1,1:8);

Week_Time = mean(sum(T)/2);

mu_wf = DistF.mu/5;
sigma_wf = DistF.sigma/5;
mu_np = 8.5/11;
sigma_np = 3.6/11;

%% Number of simulations for a fixed set of parameters
nIt = 1000;

%% Load the possible parameter for the initil FEC distribution
load ParameterOPG.mat
nc = 5;
p = mean(Par(2,:));
m = linspace(min(Par(3,:)),max(Par(3,:)),nc);
r = (m.*p)./(1 - p);

k_s = 1;

%% Number of simulations in total
nL = length(r)*nIt;

%% Storage of the ingestion risk for each goats and run
Res = zeros(nGoats,nL); 

%% Storage of the specificity for each goat
Specificity = zeros(nGoats,nL); 

%% Storage of the OPG
OPG = zeros(nGoats,nL); 

%% Storage of the feces distribution
SF = zeros(nq,nL);

%% Storage of the feces distribution, smaller scale
SF_plus = zeros(nq*25,nL);


%% It'll show a waitbar allowing to visualize the remaining number of 
%% simulation
h = waitbar(k_s/nL,'Simulation done');

for i = 1:length(r)
    
    for nS = 1:nIt
        
        %% W are the lambdas in the article. W(i,:) are the lambdas for
        %% goat i
        W = zeros(nGoats,nq);
        Spec = zeros(nGoats,1);
        
        %% We use simulation scheme 2 so first goat is not fixed
        for g = 1:nGoats
            
            a = 0.001 + rand*(0.1-0.001);
            %     a = 0.001;
            Spec(g) = a;
            Weights = pdf('geo',0:(nq-1),a); % the repartition of pasture time on the quadrats, when a=0.001 it's more or less uniform
            Weights = Weights./sum(Weights);
            W(g,:) = Weights;
        end
        
        %% We simulate the FEC when entering the first week
        OPG_Goats = random('nbin',r(i),p,[nGoats,1]);
        
        Specificity(:,k_s) = Spec;
        OPG(:,k_s) = OPG_Goats;
        
        %% Simulate risk for each goat
        [Risk , F_Dist , F_Dist_plus] = Simu_Trajectory(nq,F_W1,F_W2,Week_Time,nGoats,mu_wf,sigma_wf,mu_np,sigma_np,Evap,Precipitation,Temperature,T_Zone,ZN,OPG_Goats,W,Spec);
        Res(:,k_s) = Risk;
        SF(:,k_s) = F_Dist';
        SF_plus(:,k_s) = reshape(F_Dist_plus,[nq*25,1]);
        k_s = k_s + 1;
        waitbar(k_s/nL,h,['Running Simulations ' int2str(k_s)])
        if mod(500,k_s) == 0
            save('Res_G1.mat','Res','Specificity','OPG','SF','SF_plus')
        end
        %                 [i j k nS k_s]
    end
end


close(h) 
save('Res_G1.mat')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2
clear
load Data.mat
load Meteo.mat
load ZNumber.mat
load Distri_Weights_Feces.mat

nGoats = 21;
Groupe = 2;

if Groupe == 1
    ZN = ZN_G1;
else
    ZN = ZN_G2;
end


%% Compute the occurrence frequency for the first week
H = ComputeProportion(Groupe,Data,[1 4]);
Idx = find(H > 0);
F_W1 = H(Idx);
nq = length(F_W1);

%% Compute the occurrence frequency for the second week
H = ComputeProportion(Groupe,Data,[5 8]);
F_W2 = H(Idx);


%% Compute the average time spend in the field
T = ComputeTime(Data,Groupe,1:8);

Week_Time = mean(sum(T)/2);

mu_wf = DistF.mu/5;
sigma_wf = DistF.sigma/5;
mu_np = 8.5/11;
sigma_np = 3.6/11;


nIt = 1000;
load ParameterOPG.mat
nc = 5;
p = mean(Par(2,:));
m = linspace(min(Par(3,:)),max(Par(3,:)),nc);
r = (m.*p)./(1 - p);

k_s = 1;

nL = length(r)*nIt;
%% Storage of the ingestion risk for each goats and run
Res = zeros(nGoats,nL); 

%% Storage of the specificity for each goat
Specificity = zeros(nGoats,nL); 

%% Storage of the OPG
OPG = zeros(nGoats,nL); 

%% Storage of the feces distribution
SF = zeros(nq,nL);

%% Storage of the feces distribution, smaller scale
SF_plus = zeros(nq*25,nL);

h = waitbar(k_s/nL,'Simulation done');

for i = 1:length(r)
    
    for nS = 1:nIt
        W = zeros(nGoats,nq);
        Spec = zeros(nGoats,1);
        
        for g = 1:nGoats
            
            a = 0.001 + rand*(0.1-0.001);
            %     a = 0.001;
            Spec(g) = a;
            Weights = pdf('geo',0:(nq-1),a); % the repartition of pasture time on the quadrats, when a=0.001 it's more or less uniform
            Weights = Weights./sum(Weights);
            W(g,:) = Weights;
        end
        
        OPG_Goats = random('nbin',r(i),p,[nGoats,1]);
        
        Specificity(:,k_s) = Spec;
        OPG(:,k_s) = OPG_Goats;
        [Risk , F_Dist , F_Dist_plus] = Simu_Trajectory(nq,F_W1,F_W2,Week_Time,nGoats,mu_wf,sigma_wf,mu_np,sigma_np,Evap,Precipitation,Temperature,T_Zone,ZN,OPG_Goats,W,Spec);
        Res(:,k_s) = Risk;
        SF(:,k_s) = F_Dist';
        SF_plus(:,k_s) = reshape(F_Dist_plus,[nq*25,1]);
        k_s = k_s + 1;
        waitbar(k_s/nL,h,['Running Simulations ' int2str(k_s)])
        if mod(500,k_s) == 0
            save('Res_G2.mat','Res','Specificity','OPG','SF','SF_plus')
        end
        %                 [i j k nS k_s]
    end
end


close(h) 
save('Res_G2.mat')



