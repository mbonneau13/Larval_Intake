function Risk = Simu_Risk(nGoats,nq,Polution_tp1,G_Time,T_Zone,ZN)



Risk = zeros(nGoats,1);

for q = 1:nq
    
    Risk = Risk + NbEat(G_Time(:,q),Polution_tp1(q),T_Zone(T_Zone(:,1) == ZN(q),2));
end
