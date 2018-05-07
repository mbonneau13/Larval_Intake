%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G1
clc
clear
load Res_G1.mat
load Parameters_G1.mat

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %% Examples
% j = 3;z=3;
% % j = 4;z=1;
% J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
% R = Res(:,J);
% V = log(R(:));
% % V = R(:); 
% histogram(V,'FaceColor',[132 186 91]./255,'LineWidth',1.5)
% grid on
% xlabel('log(Risk)')
% ylabel('Frequency')
% % legend('Average FEC = 34464')
% legend('Average FEC = 1942')
% set(gca,'FontSize',16,'FontName','TimeNewsRoman')
% 
% 
% 
%% K-S test
k = 1;
Test = zeros(length(r)* ,1);
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = Res(1,J);
        V = log(R(:));
        V = (V - mean(V))./std(V);
        Test(k) = double(kstest(V));
        k = k + 1;

    end
end

Test
% 
k = 1;
M = zeros(length(T),3);

for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = Res(:,J);
        V = V(:);

        M(k,1) = mean(V);
%         M(k,2) = mean(mean(OPG(:,J),1));
        M(k,2) = ((r(j)*(1-p(z)))/p(z));
        M(k,3) = ((r(j)*(1-p(z)))/(p(z)^2));

        k = k + 1;
    end
end
X = M(:,2);
Y = M(:,1);
mdl = fitlm(X,Y)
figure('Name','G2')
hold on
plot(X,Y,'ko') 
plot(X,mdl.feval(X),'r','Linewidth',1.5)
grid on
xlabel('Expected FEC')
ylabel('Average Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')
legend('Data','Linear Regression')



% Q = zeros(5,length(r)*length(p));
% M = zeros(1,length(r)*length(p));
% k = 1;
% for j = 1:length(r)
%     for z = 1:length(p)
%         
%         J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
%         R = Res(1,J);
%         Q(1,k) = quantile(R,0.1);
%         Q(2,k) = quantile(R,0.25);
%         Q(3,k) = quantile(R,0.5);
%         Q(4,k) = quantile(R,0.75);
%         Q(5,k) = quantile(R,0.9);
%         M(k) = (r(j)*(1-p(z)))/p(z);
%         k = k + 1;
%     end
% end
% [M , I] = sort(M);
% Q(1,:) = Q(1,I);
% Q(2,:) = Q(2,I);
% Q(3,:) = Q(3,I);
% Q(4,:) = Q(4,I);
% Q(5,:) = Q(5,I);
% 
% plot(M,Q(1,:),M,Q(2,:),M,Q(3,:),M,Q(4,:),M,Q(5,:),'LineWidth',2.5)
% legend('10% quantile','25% quantile','50% quantile','75% quantile','90% quantile')
% xlabel('Mean Initial FEC')
% ylabel('Risk value')
% grid on 
% set(gca,'FontSize',16,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,...
%     'YColor',[10 10 10]./255,'Color',[189 189 189]./255,'XTick',linspace(min(M),max(M),5),'XTickLabel',num2str(linspace(min(M),max(M),5)',1))


%% Individual risk, in and out
figure 
hold on
Q = zeros(21,length(r)*2);
M = zeros(1,length(r)*2);
k = 1;
% x = -100:10:100;
Class = [0 100 1000 3000 5000];

k = 1;
for j = 1:length(r)
    for z = 1:4:length(p)
      
        M(k) = (r(j)*(1-p(z)))/p(z);
        k = k + 1;
    end
end
[M , In] = sort(M);

k = 1;
for j = 1:length(r)
    for z = 1:4:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = Res(1,J);
%         M(k) = (r(j)*(1-p(z)))/p(z);
%         pct = (ones + x./100).*M(k);
        for i = 2:length(Class)
            
            Q(i-1,In(k)) = numel(intersect(find(R >= Class(i-1)),find(R < Class(i))))./numel(R);
            c = 1 - Q(i-1,In(k));
            rectangle('Position',[(In(k)-0.5) (i - 2 - 0.5) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%             rectangle('Position',[log(M(k)) (i - 2 - 0.05) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%            plot(log(M(k)),i-1,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%             plot(log(M(k)),i-2,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        end
        Q(end,In(k)) = numel(find(R >= Class(end)))./numel(R);
        c = 1 - Q(end,In(k));
        rectangle('Position',[(In(k)-0.5) (4 - 0.5) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%         rectangle('Position',[log(M(k)) (4 - 0.5) 1 1],'FaceColor',[c c c],'EdgeColor',[0 0 0])
%         plot(log(M(k)),4,'s','MarkerFaceColor',[c c c],'MarkerEdgeColor',[0 0 0],'Markersize',10)
%         plot(log(M(k)),4,'s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1],'Markersize',max(1,floor(10*c)))
        k = k + 1;
    end
end
xlabel('Mean Initial FEC')
ylabel('Individual Risk Classes')
grid on
axis image
set(gca,'FontSize',16,'FontName','TimeNewsRoman','YTick',0:4,'YTickLabel',{'0 - 100','100 - 1000','1000 - 3000','3000 - 5000','+5000'},'XTick',1:(k-1),'XTickLabel',...
    num2str(M','%g'))

%% Same but with jet color
figure 
subplot(1,2,1);
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
        R = Res(1,J);
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
ylabel('Individual Risk Classes')
grid on
axis image
title('First Flock')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','YTick',0:4,'YTickLabel',{'0 - 100','100 - 1000','1000 - 3000','3000 - 5000','+5000'},'XTick',1:(k-1),'XTickLabel',...
    num2str(round(M)','%g'))
hcb = colorbar;
title(hcb,'Probability');
hcb.Location = 'Northoutside';


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2
clc
clear
load Res_G2.mat
load Parameters_G2.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

k = 1;
Test = zeros(length(r)*length(p),1);
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = Res(1,J);
        V = log(R(:));
        V = (V - mean(V))./std(V);
        Test(k) = double(kstest(V));
        k = k + 1;

    end
end

Test


k = 1;
M = zeros(length(r),2);

for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = Res(:,J);
        V = V(:);

        M(k,1) = mean(V);
        M(k,2) = ((r(j)*(1-p(z)))/p(z));

        k = k + 1;
    end
end
X = M(:,2);
Y = M(:,1);
mdl = fitlm(X,Y)
figure('Name','G2')
hold on
plot(X,Y,'ko') 
plot(X,mdl.feval(X),'r','Linewidth',1.5)
grid on
xlabel('Expected FEC')
ylabel('Average Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')
legend('Data','Linear Regression')




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
        R = Res(1,J);
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
ylabel('Individual Risk Classes')
grid on
axis image
title('Second Flock')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','YTick',0:4,'YTickLabel',{'0 - 100','100 - 1000','1000 - 3000','3000 - 5000','+5000'},'XTick',1:(k-1),'XTickLabel',...
    num2str(round(M)','%g'))
hcb = colorbar;
title(hcb,'Probability');
hcb.Location = 'Northoutside';
