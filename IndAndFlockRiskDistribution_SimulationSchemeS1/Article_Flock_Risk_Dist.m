%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% G1 
clc
clear
load Res_G1.mat
load Parameters_G1.mat

% Test = zeros(1,nL);
% h = waitbar(0,'Running KS test G1');
% for i = 1:nL 
% 
%     O = log(Res(:,i));
%     O = (O - mean(O))./std(O);
%     Test(i) = double(kstest(O));
%     waitbar(i/nL,h)
% end
% sum(Test)



%% Same but with jet color
figure
subplot(1,2,1)
hold on
Q = zeros(5,length(r));
M = zeros(1,length(r));
k = 1;
% x = -100:10:100;
Class = [0 100 1000 3000 5000];

k = 1;
for j = 1:length(r)
    for z = 1:length(p)
      
        M(k) = (r(j)*(1-p(z)))/p(z);
        k = k + 1;
    end
end
[M , In] = sort(M);

Color = colormap(jet(101));

k = 1;
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = mean(Res(:,J),1);
%         M(k) = (r(j)*(1-p(z)))/p(z);
%         pct = (ones + x./100).*M(k);
        k_S = find(In == k);
        for i = 2:length(Class)
            
            Q(i-1,k_S) = numel(intersect(find(R >= Class(i-1)),find(R < Class(i))))./numel(R);
            c = Q(i-1,k_S);
            c = round(100*c) + 1;
            rectangle('Position',[(k_S-0.5) (i - 2 - 0.5) 1 1],'FaceColor',Color(c,:),'EdgeColor',[0 0 0])
%             rectangle('Position',[log(M(k)) (i - 2 - 0.05) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%            plot(log(M(k)),i-1,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%             plot(log(M(k)),i-2,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        end
        Q(end,k_S) = numel(find(R >= Class(end)))./numel(R);
        c = Q(end,k_S);
        c = round(100*c) + 1;
        rectangle('Position',[(k_S-0.5) (4 - 0.5) 1 1],'FaceColor',Color(c,:),'EdgeColor',[0 0 0])
%         rectangle('Position',[log(M(k)) (4 - 0.5) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%         plot(log(M(k)),4,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%         plot(log(M(k)),4,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        k = k + 1;
    end
end
xlabel('Mean Initial FEC')
ylabel('Flock Risk Classes')
grid on
axis image
title('First Flock')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','YTick',0:4,'YTickLabel',{'0 - 100','100 - 1000','1000 - 3000','3000 - 5000','+5000'},'XTick',1:(k-1),'XTickLabel',...
    num2str(round(M)','%g'))
hcb = colorbar;
title(hcb,'Probability');
hcb.Location = 'Northoutside';




%% Test different distribution probability
% nS = 100;
% T = zeros(nL,4);
% h = waitbar(0,'Performing Permutation Test G1');
% k = 1;
% All = nL*nS;
% Par = cell(nL,4);
% tic
% for i = 1:nL
%     
%     S = Res(:,i);
%     
%     pd1 = fitdist(S,'lognormal');
%     n1 = random(pd1,[length(S) nS]);
%     
%     pd2 = fitdist(S,'Weibull');
%     n2 = random(pd2,[length(S) nS]);
%     
%     pd3 = fitdist(S,'Gamma');
%     n3 = random(pd3,[length(S) nS]);
%     
%     n4 = zeros(length(S),nS);
%     try
%         pd4 = fitdist(round(S),'nbin');
%         n4 = random(pd4,[length(S) nS]);
%     catch
%     end
%     
%     for j = 1:nS
%         
%         T(i,1) = T(i,1) + (1 - double(kstest2(n1(:,j),S)));
%         
%         T(i,2) = T(i,2) + (1 - double(kstest2(n2(:,j),S)));
%         
%         T(i,3) = T(i,3) + (1 - double(kstest2(n3(:,j),S)));
%         
%         T(i,4) = T(i,4) + (1 - double(kstest2(n4(:,j),S)));
%         
%         k = k + 1;
%         waitbar(k/All,h)
%     end
%     
%     Par{i,1} = pd1;
%     Par{i,2} = pd2;
%     Par{i,3} = pd3;
%     Par{i,4} = pd4;
% end
% toc
% 
% 
% save('FRD_G1.mat','Par','T')


% X = (PValue(1,:).*(ones - PValue(2,:)))./PValue(2,:);
X = mean(OPG,1);
Y = mean(Res,1);
I = find(Y < (mean(Y) + 10*std(Y)));
mdl = fitlm(X(I),Y(I))
figure('Name','G2')
hold on
plot(X(I),Y(I),'ko') 
plot(X(I),mdl.feval(X(I)),'r','Linewidth',2.5)
grid on
xlabel('Average FEC')
ylabel('Average Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman') 
legend('Data','Linear Regression')

x = 0:500:25000;
figure
hold on
for z = 1:length(r)
    for j = 1:length(p)
    
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        
        M = zeros(2,length(J));
        for k = 1:length(J)
            
            pd = fitdist(Res(:,J(k)),'lognormal');
            M(1,k) = pd.mu;
            M(2,k) = pd.sigma;
        end


        plot(x,pdf('lognormal',x,mean(M(1,:)),mean(M(2,:))))
    end
end

%% Examples
% j=1;z=1;
j=5;z=1;
J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));

M = zeros(2,length(J));
for k = 1:length(J)
    
    pd = fitdist(Res(:,J(k)),'lognormal');
    M(1,k) = pd.mu;
    M(2,k) = pd.sigma;
end
nS = 100;
X = random('lognormal',mean(M(1,:)),mean(M(2,:)),[20 , nS]);
x = round(linspace(min(X(:)),max(X(:)),7));
H = zeros(nS,6);
for i = 1:nS
    figure(1)
    h = histogram(X(:,i),x,'FaceColor',[132 186 91]./255,'LineWidth',1.5);
    H(i,:) = h.Values;
    close
%     i
end
% x = x(1:(end-1)) + (x(2:end) - x(1:(end-1)))/2;

figure(10)
subplot(1,2,2);hold on
bar(x(1:(end-1)),(100.*mean(H,1))./20,'FaceColor',[132 186 91]./255,'LineWidth',1.5);
grid on
xlabel('Risk')
ylabel('Pct of the flock')
% title('Average FEC = 306')
title('Average FEC = 7237')
T = [];
for i = 1:length(x)
    T = strvcat(T,num2str(round(x(i))));
end
set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255,'XTickLabel',T,'Xtick',round(x))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2
%% G2 
clc
clear
load Res_G2.mat
load Parameters_G2.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Test = zeros(1,nL);
h = waitbar(0,'Running KS test G2');
for i = 1:nL 

    O = log(Res(:,i));
    O = (O - mean(O))./std(O);
    Test(i) = double(kstest(O));
    waitbar(i/nL,h)
end
100-(100*sum(Test))/numel(Test)



%% Same but with jet color
subplot(1,2,2)
hold on
Q = zeros(5,length(r));
M = zeros(1,length(r));
k = 1;
% x = -100:10:100;
Class = [0 100 1000 3000 5000];

k = 1;
for j = 1:length(r)
    for z = 1:length(p)
      
        M(k) = (r(j)*(1-p(z)))/p(z);
        k = k + 1;
    end
end
[M , In] = sort(M);

Color = colormap(jet(101));

k = 1;
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = mean(Res(:,J),1);
%         M(k) = (r(j)*(1-p(z)))/p(z);
%         pct = (ones + x./100).*M(k);
        k_S = find(In == k);
        for i = 2:length(Class)
            
            Q(i-1,k_S) = numel(intersect(find(R >= Class(i-1)),find(R < Class(i))))./numel(R);
            c = Q(i-1,k_S);
            c = round(100*c) + 1;
            rectangle('Position',[(k_S-0.5) (i - 2 - 0.5) 1 1],'FaceColor',Color(c,:),'EdgeColor',[0 0 0])
%             rectangle('Position',[log(M(k)) (i - 2 - 0.05) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%            plot(log(M(k)),i-1,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%             plot(log(M(k)),i-2,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        end
        Q(end,k_S) = numel(find(R >= Class(end)))./numel(R);
        c = Q(end,k_S);
        c = round(100*c) + 1;
        rectangle('Position',[(k_S-0.5) (4 - 0.5) 1 1],'FaceColor',Color(c,:),'EdgeColor',[0 0 0])
%         rectangle('Position',[log(M(k)) (4 - 0.5) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%         plot(log(M(k)),4,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%         plot(log(M(k)),4,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        k = k + 1;
    end
end
xlabel('Mean Initial FEC')
ylabel('Flock Risk Classes')
grid on
axis image
title('Second Flock')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','YTick',0:4,'YTickLabel',{'0 - 100','100 - 1000','1000 - 3000','3000 - 5000','+5000'},'XTick',1:(k-1),'XTickLabel',...
    num2str(round(M)','%g'))
hcb = colorbar;
title(hcb,'Probability');
hcb.Location = 'Northoutside';



%% Load parameters with Article.m
nS = 100;
T = zeros(nL,4);
h = waitbar(0,'Performing Permutation Test G2');
k = 1;
All = nL*nS;
Par = cell(nL,4);
tic
for i = 1:nL
    
    S = Res(:,i);
    
        pd1 = fitdist(S,'lognormal');
        n1 = random(pd1,[length(S) nS]);
        
        pd2 = fitdist(S,'Weibull');
        n2 = random(pd2,[length(S) nS]);
        
        pd3 = fitdist(S,'Gamma');
        n3 = random(pd3,[length(S) nS]);
        
        
        n4 = zeros(length(S),nS);
        try
            pd4 = fitdist(round(S),'nbin');
            n4 = random(pd4,[length(S) nS]);
        catch
        end
        
        for j = 1:nS
            
            T(i,1) = T(i,1) + (1 - double(kstest2(n1(:,j),S)));
            
            T(i,2) = T(i,2) + (1 - double(kstest2(n2(:,j),S)));
            
            T(i,3) = T(i,3) + (1 - double(kstest2(n3(:,j),S)));
            
            T(i,4) = T(i,4) + (1 - double(kstest2(n4(:,j),S)));
            
            k = k + 1;
            waitbar(k/All,h)
        end
        
        Par{i,1} = pd1;
        Par{i,2} = pd2;
        Par{i,3} = pd3;
        Par{i,4} = pd4;
end
toc

save('FRD_G2.mat','Par','T')

% X = (PValue(1,:).*(ones - PValue(2,:)))./PValue(2,:);
X = mean(OPG,1);
Y = mean(Res,1);
% I = find(Y < 200000);
I = find(Y < (mean(Y) + 10*std(Y)));
mdl = fitlm(X(I),Y(I))
figure('Name','G2')
hold on
plot(X(I),Y(I),'ko') 
plot(X(I),mdl.feval(X(I)),'r','Linewidth',2.5)
grid on
xlabel('Average FEC')
ylabel('Average Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')
legend('Data','Linear Regression')