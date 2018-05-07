%% Show example
ShowPolution(Data,SF(:,1)',Idx,'Ex 1')


%% Aggregation coefficient 
%% G1
load Res_G1.mat
m = SF;
m(isnan(m)) = zeros;
M = mean(m,1);
X2 = zeros(1,nL);
P = zeros(1,nL);
for i = 1:nL
%     m(:,i) = (m(:,i) - mean(m(:,i)))./std(m(:,i));
    X2(i) = sum((m(:,i) - M(i).*ones).^2)./(M(i));
    P(i) = 1 - chi2cdf(X2(i),nq-1);
end
plot(M,X2,'o')
I = find(P > 0.5);
ShowPolution(Data,SF(:,I(1))',Idx,'Ex 1')

%% G2
load Res_G2.mat
m = SF;
m(isnan(m)) = zeros;
M = mean(m,1);
X2 = zeros(1,nL);
P = zeros(1,nL);
for i = 1:nL
%     m(:,i) = (m(:,i) - mean(m(:,i)))./std(m(:,i));
    X2(i) = sum((m(:,i) - M(i).*ones).^2)./(M(i));
    P(i) = 1 - chi2cdf(X2(i),nq-1);
end
plot(M,X2,'o')





Par = reshape(repmat(r',[1,1000])',[1,nL]);
boxplot(X2,Par)


U = unique(Par);
V = zeros(2,length(U));

for i = 1:length(U)
    
    V(1,i) = U(i);
    V(2,i) = mean(X2(Par == U(i)));
end
plot(V(1,:),V(2,:))
    
plot(Par,X2,'o')    


%% Show example of feces distribution
load Workspace_G2.mat
load Res_G2.mat
load h.mat
i = 2;

ShowPolution(Data,SF(:,3003)',Idx,'Ex 1')
set(gca,'FontSize',24,'FontName','TimeNewsRoman')
axis off
set(gca,h)


load Sampling_G1.mat
%% Pct of empty samples
(100*sum(sum(S(:,1,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,2,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,3,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,4,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,5,:) == 0)))/(size(S,1)*size(S,3))

%% Average pct of difference
nbSmax = 5;
nbR = 500;
m = SF;
m(isnan(SF))=0;
M = repmat(mean(m,1)',[1 , nbSmax , nbR]);
Diff = abs(M - S);
Diff = (100.*Diff)./M;
mean(mean(Diff,3),1)


load Sampling_G2.mat
%% Pct of empty samples
(100*sum(sum(S(:,1,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,2,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,3,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,4,:) == 0)))/(size(S,1)*size(S,3))
(100*sum(sum(S(:,5,:) == 0)))/(size(S,1)*size(S,3))

%% Average pct of difference
nbSmax = 5;
nbR = 500;
m = SF;
m(isnan(SF))=0;
M = repmat(mean(m,1)',[1 , nbSmax , nbR]);
Diff = abs(M - S);
Diff = (100.*Diff)./M;
mean(mean(Diff,3),1)


figure
load Sampling_G1.mat
load Res_G1.mat
% load Sampling_G2.mat
% load Res_G2.mat
m = SF;
m(isnan(SF)) = zeros;
m =mean(m,1);
V = [];
G = [];

for i = 1:5

    M = squeeze(S(:,i,:));
    [l , c] = find(M == 0);
    G = [G ; (25*i).*ones(length(l),1)];
    V = [V ; m(l)'];
end


subplot(1,2,1)
% subplot(1,2,2)
hold on
boxplot(V,G,'Labels',{'25','50','75','100','125'});
% for i = 1:numel(bh)
% set(bh(i),'Linewidth',2.5')
% end
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 

for i = 1:length(t)
    box1 = a(i);   % The 7th object is the first box
    name = t{i};
    switch name(1)
        case 'O'
            set(box1,'Color',[204 37 41]./255,'Linewidth',2.5');
            
        case 'M'
%             set(box1,'Color',[62 150 81]./255,'Linewidth',2.5');
            set(box1,'Color',[57 106 177]./255,'Linewidth',2.5');
            
%         case 'U'
%             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
            
%         case 'L'
%             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
            
        otherwise
            set(box1,'Color',[83 81 84]./255,'Linewidth',2.5');
    end
end      
grid on
xlabel('Number of samples')
ylabel('Average nb of L3/m^2')
title('First Flock')
% title('Second Flock')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255)


GridColor


figure
hold on
boxplot(R')
plot(1:5,mean(SF(:,1000)).*ones(1,5),'r')




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

m = SF;
m(isnan(SF))=0;
M = repmat(mean(m,1)',[1 , nbSmax , nbR]);
Diff = abs(M - S);
Diff = (100.*Diff)./M;
mean(mean(Diff,3),1)
R = squeeze(Diff(1,:,:));
boxplot(R')

