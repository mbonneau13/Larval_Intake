function Y = ComputeTime(Data,Groupe,Days)

D = Days;
Y = zeros(size(Days));
k = 1;
m = 0;
while ~isempty(D)
    
    I = intersect(find(Data.Day == D(1)),find(Data.Groupe == Groupe));
    Time = unique(Data.Time(I));
    y = zeros(1,numel(Time));
    x = zeros(1,numel(Time));
    
    for i = 1:numel(Time)
        
        y(i) = sum(Data.Time(I) == Time(i));
        x(i) = Time(i);
    end
    if max(y) > m
        m = max(y);
    end
    Y(k) = trapz(x,y);
    k = k + 1;
    D = setdiff(D,D(1));
end
% disp(['Maximal Number of Carbis -> ' num2str(m)])

