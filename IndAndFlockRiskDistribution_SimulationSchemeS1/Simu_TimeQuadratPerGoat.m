function G_Time = Simu_TimeQuadratPerGoat(nq,W,gt,Q_Time,nGoats)
%% Initialize time in each quadrat for each goat
G_Time = zeros(nGoats,nq);


%% T_All is the time spend per goat in some quadrat. Here the column index of T_All are fake.
%% The true distribution of time will be define below.
T_All = zeros(nGoats,nq);

for g = 1:nGoats
    
    T_All(g,:) = W(g,:).*gt(g);
end

%% Now we start to simulate the time distribution over the quadrat, using T_All
K = ones(1,nGoats);
T_Init = sum(Q_Time);

while ~isequal(floor(sum(sum(G_Time))*10),floor(10*T_Init))
    
    g = randi(nGoats);
    k = K(g);
    T = T_All(g,:);
    
    while (k <= nq) && (sum(Q_Time) > 0)
        
        IQ = find(Q_Time > T(k));
        if isempty(IQ)
            IQ = find(Q_Time > 0);
        end
        %% We choose a quadrat
        q = IQ(randi(numel(IQ)));
        
        if (Q_Time(q) >= T(k))
            
            G_Time(g,q) = G_Time(g,q) + T(k);
            Q_Time(q) = Q_Time(q) - T(k);
            T(k) = 0;
            k = k + 1;
        else
            G_Time(g,q) = G_Time(g,q) + Q_Time(q);
            T(k) = T(k) - Q_Time(q);
            Q_Time(q) = 0;
        end
    end
    K(g) = k;
    T_All(g,:) = T;
end