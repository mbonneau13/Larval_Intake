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

