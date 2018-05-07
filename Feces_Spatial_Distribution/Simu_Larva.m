function Polution_tp1 = Simu_Larva(Polution_t,Temperature,Evap,Precipitation)

%% Simulate survival of larvae during 31 days

Nt = zeros(3,length(Polution_t));
Nt(1,:) = Polution_t;

for t = 5:length(Temperature)
    T = Temperature(t);
    E = Evap(t);
    P = Precipitation(t);
    
    delta = -0.09746 + 0.01063*T;
    mu1 = exp(-1.62026 -0.17771*T + 0.00629*T^2);
    mu2 = exp(-1.82300 -0.14180*T + 0.00405*T^2);
    mu3 = exp(-2.63080 -0.14407*T + 0.00463*T^2);
    
    m1 = 0.051;
    if P >= 2
        m1 = 0.25;
    end
    
    if (P<2) && (sum(Precipitation((t-4):t)./Evap((t-4):t)) < 1)
        
        m1 = 0;
    end
    
    A = [(1 - mu1 - 2*delta) , 0 , 0;...
        (2*delta) , (1 - mu2 - 2*delta) , 0;...
        0 , (2*delta) , (1 - mu3 - m1)];
    Nt = A*Nt;
end
Polution_tp1 = Nt(end,:);
