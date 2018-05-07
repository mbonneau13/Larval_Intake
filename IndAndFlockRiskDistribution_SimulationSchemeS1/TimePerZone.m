%% Compute Time To Eat All
load Data.mat

%% Groupe 1, Week 1
nGoats = 21;

%% Compute the occurrence frequency
H = ComputeProportion(1,Data,[1 4]);
Idx = find(H > 0);
F = H(Idx);
nq = length(F);

%% Compute the average time spend in the field
T = ComputeTime(Data,1,1:8);

Week_Time = mean(sum(T)/2);

%% Simulate time in pasture for each goat
gt = Week_Time/21.*ones(1,21);%normrnd(Week_Time/21,0.1*(Week_Time/21),[1 , 21]);
Q_Time = F.*sum(gt);

%% Determine zone number of each quadrat
U = unique(Data.Idx);
U = U(Idx);
ZN = zeros(size(U));
for zn = 1:length(U)
    
    I = find(Data.Idx == U(zn));
    ZN(zn) = mode(Data.ZNumber(I));
end

ZN_G1 = ZN;

U = unique(ZN);
H = hist(ZN,U);

Time_zn = zeros(size(U));
for zn = 1:length(H)
    
    if H(zn) > 25
        
        q = quantile(Q_Time(ZN == U(zn)),0.95);
        t = Q_Time(ZN == U(zn));
        t = t(t >= q);
        Time_zn(zn) = 1/mean(t);
    end
end

for zn = 1:length(H)
    
    if H(zn) <= 25
        
        Time_zn(zn) = mean(Time_zn(Time_zn > 0));
    end
end


T_Zone = [U , Time_zn];




%% Groupe 1, Week 1
nGoats = 21;

%% Compute the occurrence frequency
H = ComputeProportion(2,Data,[1 4]);
Idx = find(H > 0);
F = H(Idx);
nq = length(F);

%% Compute the average time spend in the field
T = ComputeTime(Data,2,1:8);

Week_Time = mean(sum(T)/2);

%% Simulate time in pasture for each goat
gt = Week_Time/21.*ones(1,21);%normrnd(Week_Time/21,0.1*(Week_Time/21),[1 , 21]);
Q_Time = F.*sum(gt);

%% Determine zone number of each quadrat
U = unique(Data.Idx);
U = U(Idx);
ZN = zeros(size(U));
for zn = 1:length(U)
    
    I = find(Data.Idx == U(zn));
    ZN(zn) = mode(Data.ZNumber(I));
end

ZN_G2 = ZN;

U = unique(ZN);
H = hist(ZN,U);

Time_zn = zeros(size(U));
for zn = 1:length(H)
    
    if H(zn) > 25
        
        q = quantile(Q_Time(ZN == U(zn)),0.95);
        t = Q_Time(ZN == U(zn));
        t = t(t >= q);
        Time_zn(zn) = 1/mean(t);
    end
end

for zn = 1:length(H)
    
    if H(zn) <= 25
        
        Time_zn(zn) = mean(Time_zn(Time_zn > 0));
    end
end

T_Zone = [T_Zone ; U , Time_zn];


U = unique(T_Zone(:,1));

T_Zone2 = [];
for i = 1:length(U)
    
    T_Zone2 = [T_Zone2 ; U(i) , mean(T_Zone(T_Zone(:,1) == U(i),2))];
end

T_Zone = T_Zone2;

save('ZNumber.mat','T_Zone','ZN_G1','ZN_G2')




