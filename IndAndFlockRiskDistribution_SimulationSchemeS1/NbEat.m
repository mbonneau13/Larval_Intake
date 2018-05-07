function Res = NbEat(Time,NbL,t_zn)

%% V3
nL = 50;
Time = round(Time);
Res = zeros(nL,numel(Time));

if NbL > 0
    
    Position = [];
    for g = 1:numel(Time)
        
        if Time(g) > 0
            Position = [Position , g.*ones(1,Time(g))];
        end
    end
    
    for i = 1:nL
        
        
        P = Position(randperm(numel(Position)));
        
        Q = (t_zn*NbL).*ones(size(Position));
        Qcs = cumsum(Q);
        I = find(Qcs <= NbL,1,'Last');
        P = P(1:I);
        
        for g = 1:numel(Time)
            Res(i,g) = sum(P == g)*t_zn*NbL;
        end
    end
end

Res = mean(Res,1)';


%% V2
% Time = round(Time)*60;
%
% nL = 100;
%
% Res = zeros(nL,numel(Time));
%
%
% if NbL > 0
%     for i = 1:nL
%         Position = [];
%
%         for g = 1:numel(Time)
%
%             if Time(g) > 0
%                 Position = [Position , g.*ones(1,Time(g))];
%             end
%         end
%
%         Position = Position(randperm(numel(Position)));
%         N = binornd(1,p,1,numel(Position));
%
%         CN = cumsum(N);
%         idx = find(CN <= NbL,1,'last');
%         N = N(1:idx);
%         Position = Position(1:idx);
%
%         for g = 1:numel(Time)
%
%             Res(g) = sum(N(Position == g));
%         end
%     end
%     Res= mean(Res,1);
%
% else
%     Res = zeros(1,numel(Time));
% end

%% V1

%
% [~ , I] = sort(-Time);
%
% for i = 1:length(I)
%
%     [~ , idx] = min(abs(x-Time(i).*ones));
%     Res(I(i)) = k_Prop(idx)*Polution;
%     Polution = max(0,Polution - Res(I(i)));
% end





