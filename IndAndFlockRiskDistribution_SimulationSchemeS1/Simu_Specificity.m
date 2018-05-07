function [W , Spec] = Simu_Specificity(nq,nGoats)

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