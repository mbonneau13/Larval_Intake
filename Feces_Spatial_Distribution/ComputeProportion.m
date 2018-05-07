function H = ComputeProportion(Groupe,Data,Days)


U = unique(Data.Idx);

I = find(Data.Groupe == Groupe);
Time = unique(Data.Time(I));
H = zeros(1,length(U));
H2 = 0;

for Day = Days(1):Days(2)
    
    ind = intersect(I,find(Data.Day == Day));
    Time = unique(Data.Time(ind));
    for i = 1:numel(Time)
        
        F = find(Data.Time(ind) == Time(i));
        h = hist(Data.Idx(I(F)),U);
        if sum(h > 0) > 0
            H2 = H2 + 1;
        end
        h(h>0) = ones;
        H = H + h;
    end
end

if Groupe == 2
    XY = [2987 , 2988 , 2989 , 3061 , 3062 , 3063 , 3064 , 3136 , 3137 , 3138 , 3139 , 3212 , 3213];
    for i = 1:length(XY)
        H(U == XY(i)) = zeros;
    end
end

H = H./sum(H);

% J = U(H>0)';
% c = H(H>0);

% J = U(H>0)';
% c = H(H>0);
% x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
% y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];
% patch(x,y,c,'EdgeColor','w')
% C = rot90(rot90(colormap(gray(50))));
% colormap(C)
% caxis([0 0.15])
% colorbar