%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G1
clear
load Res_G1.mat
load Parameters_G1.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = ttest2(x,y) returns a test decision for the null hypothesis that the 
% data in vectors x and y comes from independent random samples from normal 
% distributions with equal means and equal but unknown variances, 
% using the two-sample t-test. The alternative hypothesis is that the data 
% in x and y comes from populations with unequal means. The result h is 1 
% if the test rejects the null hypothesis at the 5% significance level, 
% and 0 otherwise.


%% Code is used to determine if the specificity has an impact on the individual average mean risk
k = 1;
Test = -ones(length(r),1);
M = zeros(5,5);
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = [];
        G = [];
        for i = 1:length(a_Init)
            
            I = intersect(find(Specificity(1,:) == a_Init(i)),J);
            R = Res(1,I)';
            V = [V , R];
%             R = log(R(:));
%             V = [V ; R];
%             G = [G ; round(a_Init(i),2).*ones(size(R))];
            M(j,i) = quantile(R(:),0.75);
        end
        Test(j) = kruskalwallis(V,[],'off');

%% Comment if no boxplot
%         subplot(3,2,j);hold on
% %                 bh = boxplot(V,G,'Labels',num2str(a_Init')); % With outliers
%         bh = boxplot(V,G,'Labels',num2str(a_Init'),'Symbol',''); % Without outliers
% %         Without outliers
%         w = 1.5;
%         q3 = quantile(V,0.75);
%         q1 = quantile(V,0.25);
%         L = q3 + w*(q3-q1) + 1;
%         ylim([0 L]) % Comment when showing outliers
% 
%         a = get(get(gca,'children'),'children');   % Get the handles of all the objects
%         t = get(a,'tag');   % List the names of all the objects
%         
%         for i = 1:length(t)
%             box1 = a(i);   % The 7th object is the first box
%             name = t{i};
%             switch name(1)
%                 case 'O'
%                     set(box1,'Color',[204 37 41]./255,'Linewidth',2.5');
%                     
%                 case 'M'
%                     %             set(box1,'Color',[62 150 81]./255,'Linewidth',2.5');
%                     set(box1,'Color',[57 106 177]./255,'Linewidth',2.5');
%                     
%                     %         case 'U'
%                     %             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
%                     
%                     %         case 'L'
%                     %             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
%                     
%                 otherwise
%                     set(box1,'Color',[83 81 84]./255,'Linewidth',2.5');
%             end
%         end
%         grid on
%         xlabel('Specificity')
%         ylabel('log(Risk)')
%         set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255)
%         title(['Mean Initial FEC = ' num2str(((r(j)*(1-p(z)))/p(z)),1)])

    end
end

Test

figure
subplot(1,2,1);hold on
title('First Flock')
plot(a_Init,M(1,:),'-s',a_Init,M(2,:),'-s',a_Init,M(3,:),'-s',a_Init,M(4,:),'-s',a_Init,M(5,:),'-s','LineWidth',2)
grid on
xlabel('Specificity')
% ylabel('Mean Individual Risk')
ylabel('75% Quantile Individual Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255,...
    'XTick',a_Init)

%%% Mean risk as a function of FEC parameters
k = 1;
M = zeros(length(r)*length(p),2);

for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = Res(:,J);
        R = R(:);

        M(k,1) = mean(R);
%         M(k,1) = quantile(R,0.95);
        M(k,2) = ((r(j)*(1-p(z)))/p(z));

        k = k + 1;
    end
end
plot(M(:,2),M(:,1),'o')
mdl = fitlm(M(:,2),M(:,1))






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2
clc
clear all
load Res_G2.mat
load Parameters_G2.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% h = ttest2(x,y) returns a test decision for the null hypothesis that the 
% data in vectors x and y comes from independent random samples from normal 
% distributions with equal means and equal but unknown variances, 
% using the two-sample t-test. The alternative hypothesis is that the data 
% in x and y comes from populations with unequal means. The result h is 1 
% if the test rejects the null hypothesis at the 5% significance level, 
% and 0 otherwise.

%% Code is used to determine if the specificity has an impact on the individual average mean risk
k = 1;
Test = -ones(length(r),1);
M = zeros(5,5);
P1 = zeros(5,5);
P2 = zeros(5,5);
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = [];
        G = [];
        for i = 1:length(a_Init)
            
            I = intersect(find(Specificity(1,:) == a_Init(i)),J);
            R = Res(1,I)';
            V = [V , R];
%             R = log(R(:));
%             V = [V ; R];
%             G = [G ; round(a_Init(i),2).*ones(size(R))];
            M(j,i) = quantile(R(:),0.75);
            P1(j,i) = a_Init(i);
            P2(j,i) = (r(j)*(1-p))/p;
        end
        Test(j) = kruskalwallis(V,[],'off');

%% Comment if no boxplot
%         subplot(3,2,j);hold on
% %                 bh = boxplot(V,G,'Labels',num2str(a_Init')); % With outliers
%         bh = boxplot(V,G,'Labels',num2str(a_Init'),'Symbol',''); % Without outliers
% %         Without outliers
%         w = 1.5;
%         q3 = quantile(V,0.75);
%         q1 = quantile(V,0.25);
%         L = q3 + w*(q3-q1) + 1;
%         ylim([0 L]) % Comment when showing outliers
% 
%         a = get(get(gca,'children'),'children');   % Get the handles of all the objects
%         t = get(a,'tag');   % List the names of all the objects
%         
%         for i = 1:length(t)
%             box1 = a(i);   % The 7th object is the first box
%             name = t{i};
%             switch name(1)
%                 case 'O'
%                     set(box1,'Color',[204 37 41]./255,'Linewidth',2.5');
%                     
%                 case 'M'
%                     %             set(box1,'Color',[62 150 81]./255,'Linewidth',2.5');
%                     set(box1,'Color',[57 106 177]./255,'Linewidth',2.5');
%                     
%                     %         case 'U'
%                     %             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
%                     
%                     %         case 'L'
%                     %             set(box1,'Color',[218 124 48]./255,'Linewidth',2.5');
%                     
%                 otherwise
%                     set(box1,'Color',[83 81 84]./255,'Linewidth',2.5');
%             end
%         end
%         grid on
%         xlabel('Specificity')
%         ylabel('log(Risk)')
%         set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255)
%         title(['Mean Initial FEC = ' num2str(((r(j)*(1-p(z)))/p(z)),1)])

    end
end

Test

subplot(1,2,2);hold on
title('Second Flock')
plot(a_Init,M(1,:),'-s',a_Init,M(2,:),'-s',a_Init,M(3,:),'-s',a_Init,M(4,:),'-s',a_Init,M(5,:),'-s','LineWidth',2)
grid on
xlabel('Specificity')
% ylabel('Mean Individual Risk')
ylabel('75% quantile Individual Risk')
h = legend(num2str((r(1)*(1-p))/p,1),num2str((r(2)*(1-p))/p,1),num2str((r(3)*(1-p))/p,1),num2str((r(4)*(1-p))/p,1),num2str((r(5)*(1-p))/p,1));
title(h,'Initial Mean FEC')
set(gca,'FontSize',16,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,'YColor',[10 10 10]./255,'Color',[195 195 195]./255,...
    'XTick',a_Init)


%%% Mean risk as a function of FEC parameters
k = 1;
M = zeros(length(r)*length(p),2);

for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        R = Res(1,J);
        R = R(:);

%         M(k,1) = mean(R);
        M(k,1) = quantile(R,0.95);
        M(k,2) = ((r(j)*(1-p(z)))/p(z));

        k = k + 1;
    end
end
plot(M(:,2),M(:,1),'o')
mdl = fitlm(M(:,2),M(:,1))