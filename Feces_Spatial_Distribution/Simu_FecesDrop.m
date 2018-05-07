function [Polution_t , Polution_plus] = Simu_FecesDrop(nGoats,nq,mu_wf,sigma_wf,mu_np,sigma_np,G_Time,gt,OPG_Goats)
%% Simulate the weights of feces
FEC_Weights = 4.*(random('Normal',mu_wf,sigma_wf,1,nGoats));
FEC_Weights = FEC_Weights.*1000;
FEC_Weights(FEC_Weights < 0) = zeros;
%     FEC_Weights = 4.*(600.*ones + (900 - 600).*rand(1,nGoats));
Polution_t = zeros(1,nq);
Polution_plus = zeros(25,nq);

for g = 1:nGoats
    
    %% Simulate the number of packs
    %         nP = sum(40.*ones + rand(1,4)*(60 - 40));
    nP = random('Normal',mu_np,sigma_np);
    nP(nP < 0) = zeros;
    nP = nP*24*4;
    propT = gt(g)/(24*60*4);
    w = FEC_Weights(g)/nP;
    nP = max(round(propT*nP),1); %% We only consider clump during time on the pasture
    %% Select the quadrats, respecting the propotion of time spend inside the quadrat
    Q = gendist(G_Time(g,:),1,nP);
    for i = 1:length(Q)
        
        Polution_t(Q(i)) = Polution_t(Q(i)) +  OPG_Goats(g)*w;
        Polution_plus(randi(25),Q(i)) = Polution_plus(randi(25),Q(i)) + OPG_Goats(g)*w;
    end
end