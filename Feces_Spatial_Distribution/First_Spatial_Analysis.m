% %% Compute the (x,y) axis of the quadrats, store it in x and y
% H = ComputeProportion(Groupe,Data,[1 4]);
% Idx = find(H > 0);
% F_W1 = H(Idx);
% nq = length(F_W1);
% 
% U = unique(Data.Idx);
% U = U(Idx);
% 
% 
% J = U(F_W1>0)';
% c = H(F_W1>0);
% x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
% y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];
% 
% x = mean(x,1)'; %% Take the middle of the quadrat as a coordinate
% y = mean(y,1)';
% 
% %% Fit the intensity of the inhomogeneous poisson point process
% % Simulate feces distribution
% % load Res_G1.mat
% % mdl = fitlm([x , y],log(SF(:,1) + ones));
% % lambda = SF(1,:);
% 
% 
% % 
% % %% Compute the Cartesian coordinates of the pasture's border
% Border = [P1.Position(1) , P1.Position(2) ;  P2.Position(1) , P2.Position(2) ; P3.Position(1) , P3.Position(2) ;...
%     P4.Position(1) , P4.Position(2) ; P5.Position(1) , P5.Position(2) ; P6.Position(1) , P6.Position(2)];
% Border = [P1.Position(1) , P1.Position(2) ;  P2.Position(1) , P2.Position(2) ; P3.Position(1) , P3.Position(2) ;...
%     P4.Position(1) , P4.Position(2) ];
% % load Border_G1.mat
% % pgon_P = polyshape(Border(:,1),Border(:,2));
% % Surf_P = polyarea(Border(:,1),Border(:,2));
% % 
% Dist = zeros(length(x),length(x));
% 
% for i = 1:length(x)
%     for j = (i+1):length(x)
%         x1 = [x(i),y(i)];
%         x2 = [x(j) , y(j)];
%         Dist(i,j) = norm(x1 - x2);
%         Dist(j,i) = Dist(i,j);
%     end
% end
% % 
% % theta = 0:(pi/100):(2*pi);
% % W = zeros(length(x),length(x));
% % 
% % for i = 1:length(x)
% %     for j = 1:length(x)
% %         
% %         if ~isequal(i,j)
% %             [x_c , y_c] = pol2cart(theta,Dist(i,j).*ones(size(theta)));
% %             pgon_C = polyshape(x_c + x(i),y_c + y(i));
% %             in = inpolygon(pgon_C.Vertices(:,1),pgon_C.Vertices(:,2),pgon_P.Vertices(:,1),pgon_P.Vertices(:,2));
% %             W(i,j) = sum(in)/numel(in);
% %             if W(i,j)==0
% %                 disp('u')
% %             end
% %         end
% %     end
% % end
% % W(W==0) = min(min(W(W > 0)));
% load W_G1.mat
% 
% 
% %% Compute the F-Function
% K = zeros(1,70);
% NN = find(lambda > 0);
% for t = 1:70
%     s = 0;
%     for i = 1:length(x)
%         
%         if lambda(i) > 0
%             I = intersect(setdiff(find(Dist(i,:) <= t),i),NN);
%             if ~isempty(I)
%                 %         L_i = exp(mdl.feval([x(i) , y(i)]));
%                 %         L_j = exp(mdl.feval([x(I) , y(I)]));
%                 L_i = lambda(i);
%                 L_j = lambda(I);
%                 w_ij = W(i,I);
%                 
%                 r = sum(w_ij./(L_j.*L_i));
%                 s = s + r/L_i;
%             end
%         end
% %         [i s]
%     end
%     K(t) = s/Surf_P;
% end
% plot(1:70,K)
% 
% 
% %% Compute the experimental variogram
% figure
% hold on
% d=1;
% for k = 1:100
% n = randi(5000);
% dmax = 30;
% F = SF(:,n);
% % F = (F - mean(F).*ones)./std(F);
% F = F./std(F);
% C = -ones(length(x),length(x));
% Value = zeros(size(C));
% d1 = (d - d/2):d:70.5;
% d1 = [0 d1];
% 
% 
% for i = 1:length(x)
%     if F(i) > 0
%         for j = i:length(x)
%             if F(j) > 0
%                 x1 = [x(i),y(i)];
%                 x2 = [x(j) , y(j)];
%                 z = norm(x1 - x2);
%                 C(i,j) = find(z >= d1,1,'last');
%                 Value(i,j) = (F(i) - F(j))^2;
%             end
%         end
%     end
% end
% 
% h = 0:dmax;
% V = zeros(dmax+1,3);
% 
% for h = 0:dmax
% 
%     V(h + 1,1) = h;
%     V(h + 1,2) = sum(Value(C == (h + 1)));
%     V(h + 1,3) = sum(sum(C == (h + 1)));
% end
% 
% % figure('name',int2str(n))
% plot(0:dmax,V(:,2)./(2.*V(:,3)))
% end
% 
% 
% % Conclusion, it's a difficult spatial distribution to characterize due to
% % the extreme discontinuity. 
% 
% 
%% G1
load Res_G1.mat
load Dist_G1.mat
[A,R] = geotiffread('CarteBase.tif');

%% We first divide the pasture into subplot for regular sampling
nsz = 25; % number of sampling zone
load Border_G1.mat
load Tree.mat
qx = [Data.Cx1(:)' ; Data.Cx2(:)' ; Data.Cx3(:)' ; Data.Cx4(:)'];
qy = [Data.Cy1(:)' ; Data.Cy2(:)' ; Data.Cy3(:)' ; Data.Cy4(:)'];
qx = mean(qx,1)';
qy = mean(qy,1)';

in = inpolygon(qx,qy,Border(:,1),Border(:,2));
qx = qx(in);
qy = qy(in);
in = inpolygon(qx,qy,Tree(:,1),Tree(:,2));
qx = qx(~in);
qy = qy(~in);
SZN = kmeans([qx , qy],nsz); % Sampling Zone Number

% It works !
figure
mapshow(A,R);hold on
V = zeros(1,nsz);
Pog = cell(nsz,1);
for i = 1:nsz
    
    J = find(SZN == i);
    [K , v] = convhull(qx(J),qy(J));
    V(i) = v;
    po = polyshape(qx(J(K)),qy(J(K)));
    plot(po)
    Pog{i,1} = po;
%     pause
end
V


%% We store the szn for each record

% First load the coordinates of each quadrat
H = ComputeProportion(Groupe,Data,[1 4]);
Idx = find(H > 0);
F_W1 = H(Idx);
nq = length(F_W1);

U = unique(Data.Idx);
U = U(Idx);


J = U(F_W1>0)';
c = H(F_W1>0);
x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];

x = mean(x,1)'; %% Take the middle of the quadrat as a coordinate
y = mean(y,1)';

ZN = -ones(size(x));

for i = 1:length(x)
    
    in = 0;
    k = 1;
    while (k <= nsz) && (in == 0)

        po = Pog{k,1};
        in = inpolygon(x(i),y(i),po.Vertices(:,1),po.Vertices(:,2));
        k = k + 1;
    end
    if in == 1
        ZN(i) = k - 1;
    end 
end

I = find(ZN == -1);
for i = 1:sum(ZN == -1)

    [~ , J] = sort(Dist(I(i),:));
    ZN(I(i)) = ZN(J(find(ZN(J) ~= -1,1,'first')));
end

% It works !
mapshow(A,R);hold on
for i = 1:length(x)
    text(x(i),y(i),num2str(ZN(i)))
end

nbSmax = 5;
nbR = 500;
S = zeros(nL,nbSmax,nbR);
h = waitbar(0,'running');
SF(isnan(SF)) = zeros;
Pos = reshape(1:(nq*25),[25 , nq]);

for i = 1:nL
    Fec = reshape(SF_plus(:,i)',[25 , nq]);
    Fec(isnan(Fec)) = zeros;
    for nbS = 1:nbSmax
        
        s = zeros(1,nbR);
        for n = 1:nbR
            for z = 1:nsz
                
                I = find(ZN == z);
                Pos_p = Pos(:,I);
                q = randperm(numel(Pos_p),nbS);
                A = Fec(Pos_p(q));
                s(1,n) = s(1,n) + sum(A);
            end
        end
        S(i,nbS,:) = s./(nbS*nsz);
    end
    waitbar(i/nL,h)
end
close(h)

save('Sampling_G1.mat','S')

R_1 = [];
for i = 1:1000
    R_1 = [R_1 ; transpose(squeeze(S(i,:,:).*25))];
end

R_2 = [];
for i = 1001:2000
    R_2 = [R_2 ; transpose(squeeze(S(i,:,:).*25))];
end

R_3 = [];
for i = 2001:3000
    R_3 = [R_3 ; transpose(squeeze(S(i,:,:).*25))];
end

R_4 = [];
for i = 3001:4000
    R_4 = [R_4 ; transpose(squeeze(S(i,:,:).*25))];
end

R_5 = [];
for i = 4001:5000
    R_5 = [R_5 ; transpose(squeeze(S(i,:,:).*25))];
end

figure(1)
boxplot(R_1)

figure(2)
boxplot(R_2)

figure(3)
boxplot(R_3)

figure(4)
boxplot(R_4)

figure(5)
boxplot(R_5)

100 - (100*sum(sum(S(:,1,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,2,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,3,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,4,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,5,:) == 0)))/(size(S,1)*size(S,3))



figure
hold on
boxplot(R')
plot(1:5,mean(SF(:,1000)).*ones(1,5),'r')

M = repmat(mean(SF,1)',[1 , nbSmax , nbR]);
Diff = abs(M - S);
Diff = (100.*Diff)./M;
R = squeeze(Diff(1,:,:));
boxplot(R')












%% G2
clear all
load Res_G2.mat
load Dist_G2.mat
[A,R] = geotiffread('CarteBase.tif');

%% We first divide the pasture into subplot for regular sampling
nsz = 25; % number of sampling zone
load Border_G2.mat
qx = [Data.Cx1(:)' ; Data.Cx2(:)' ; Data.Cx3(:)' ; Data.Cx4(:)'];
qy = [Data.Cy1(:)' ; Data.Cy2(:)' ; Data.Cy3(:)' ; Data.Cy4(:)'];
qx = mean(qx,1)';
qy = mean(qy,1)';

in = inpolygon(qx,qy,Border(:,1),Border(:,2));
qx = qx(in);
qy = qy(in);
load Tache_G2.mat
in = inpolygon(qx,qy,Tree(:,1),Tree(:,2));
qx = qx(~in);
qy = qy(~in);

SZN = kmeans([qx , qy],nsz); % Sampling Zone Number

% It works !
figure
mapshow(A,R);hold on
V = zeros(1,nsz);
Pog = cell(nsz,1);
for i = 1:nsz
    
    J = find(SZN == i);
    [K , v] = convhull(qx(J),qy(J));
    V(i) = v;
    po = polyshape(qx(J(K)),qy(J(K)));
    plot(po)
    Pog{i,1} = po;
%     pause
end
V


%% We store the szn for each record

% First load the coordinates of each quadrat
H = ComputeProportion(2,Data,[1 4]);
Idx = find(H > 0);
F_W1 = H(Idx);
nq = length(F_W1);

U = unique(Data.Idx);
U = U(Idx);


J = U(F_W1>0)';
c = H(F_W1>0);
x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];

x = mean(x,1)'; %% Take the middle of the quadrat as a coordinate
y = mean(y,1)';

ZN = -ones(size(x));

for i = 1:length(x)
    
    in = 0;
    k = 1;
    while (k <= nsz) && (in == 0)

        po = Pog{k,1};
        in = inpolygon(x(i),y(i),po.Vertices(:,1),po.Vertices(:,2));
        k = k + 1;
    end
    if in == 1
        ZN(i) = k - 1;
    end 
end

I = find(ZN == -1);
for i = 1:sum(ZN == -1)

    [~ , J] = sort(Dist(I(i),:));
    ZN(I(i)) = ZN(J(find(ZN(J) ~= -1,1,'first')));
end

% It works !
mapshow(A,R);hold on
for i = 1:length(x)
    text(x(i),y(i),num2str(ZN(i)))
end

nbSmax = 5;
nbR = 500;
S = zeros(nL,nbSmax,nbR);
h = waitbar(0,'running');
SF(isnan(SF)) = zeros;
Pos = reshape(1:(nq*25),[25 , nq]);

for i = 1:nL
    Fec = reshape(SF_plus(:,i)',[25 , nq]);
    Fec(isnan(Fec)) = zeros;
    for nbS = 1:nbSmax
        
        s = zeros(1,nbR);
        for n = 1:nbR
            for z = 1:nsz
                
                I = find(ZN == z);
                Pos_p = Pos(:,I);
                q = randperm(numel(Pos_p),nbS);
                A = Fec(Pos_p(q));
                s(1,n) = s(1,n) + sum(A);
            end
        end
        S(i,nbS,:) = s./(nbS*nsz);
    end
    waitbar(i/nL,h)
end
close(h)
save('Sampling_G2.mat','S')


R_1 = [];
for i = 1:1000
    R_1 = [R_1 ; transpose(squeeze(S(i,:,:).*25))];
end

R_2 = [];
for i = 1001:2000
    R_2 = [R_2 ; transpose(squeeze(S(i,:,:).*25))];
end

R_3 = [];
for i = 2001:3000
    R_3 = [R_3 ; transpose(squeeze(S(i,:,:).*25))];
end

R_4 = [];
for i = 3001:4000
    R_4 = [R_4 ; transpose(squeeze(S(i,:,:).*25))];
end

R_5 = [];
for i = 4001:5000
    R_5 = [R_5 ; transpose(squeeze(S(i,:,:).*25))];
end

figure(1)
boxplot(R_1)

figure(2)
boxplot(R_2)

figure(3)
boxplot(R_3)

figure(4)
boxplot(R_4)

figure(5)
boxplot(R_5)

100 - (100*sum(sum(S(:,1,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,2,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,3,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,4,:) == 0)))/(size(S,1)*size(S,3))
100 - (100*sum(sum(S(:,5,:) == 0)))/(size(S,1)*size(S,3))



figure
hold on
boxplot(R')
plot(1:5,mean(SF(:,1000)).*ones(1,5),'r')

M = repmat(mean(SF,1)',[1 , nbSmax , nbR]);
Diff = abs(M - S);
Diff = (100.*Diff)./M;
R = squeeze(Diff(1,:,:));
boxplot(R')