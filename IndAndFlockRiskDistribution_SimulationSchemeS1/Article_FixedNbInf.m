load Res_G1_FixedInfected.mat

Par = zeros(nL,2);
cpt = 1;
for nI = 1:5
    for k = 1:length(a_Init)
        for nS = 1:nIt

            Par(cpt,1) = nI;
            Par(cpt,2) = a_Init(k);
            cpt = cpt + 1;
        end
    end
end

[~ , Bin]=hist(Res(1,:),50);

figure
cpt = 1;
for nI = 1:5
    for k = 1:length(a_Init)
%         figure('name',['nI = ' int2str(nI) ', k = ' num2str(a_Init(k))])
        subplot(5,5,cpt);hold on
        title(['nI = ' int2str(nI) ', k = ' num2str(a_Init(k))])
        I = intersect(find(Par(:,1) == nI),find(Par(:,2) == a_Init(k)));
        histogram(Res(1,I))
        cpt = cpt + 1;
    end
end

q25 = zeros(5,5);
qM = zeros(5,5);
q75 = zeros(5,5);

cpt = 1;
for nI = 1:5
    for k = 1:length(a_Init)

        I = intersect(find(Par(:,1) == nI),find(Par(:,2) == a_Init(k)));
        q25(nI,k) = quantile(Res(1,I),0.25);
        q75(nI,k) = quantile(Res(1,I),0.75);
        qM(nI,k) = mean(Res(1,I));
    end
end

figure
subplot(1,3,1);hold on
title('25% quantile')
plot(a_Init,q25(1,:),a_Init,q25(2,:),a_Init,q25(3,:),a_Init,q25(4,:),a_Init,q25(5,:),'LineWidth',3)
h = legend('1','2','3','4','5');
title(h,'Number of infected goat')
grid on

subplot(1,3,2);hold on
title('Mean')
plot(a_Init,qM(1,:),a_Init,qM(2,:),a_Init,qM(3,:),a_Init,qM(4,:),a_Init,qM(5,:),'LineWidth',3)
h = legend('1','2','3','4','5');
title(h,'Number of infected goat')
grid on

subplot(1,3,3);hold on
title('75% quantile')
plot(a_Init,q75(1,:),a_Init,q75(2,:),a_Init,q75(3,:),a_Init,q75(4,:),a_Init,q75(5,:),'LineWidth',3)
h = legend('1','2','3','4','5');
title(h,'Number of infected goat')
grid on



q25 = zeros(5,5);
qM = zeros(5,5);
q75 = zeros(5,5);

cpt = 1;
for nI = 1:5
    for k = 1:length(a_Init)

        I = intersect(find(Par(:,1) == nI),find(Par(:,2) == a_Init(k)));
        q25(k,nI) = quantile(Res(1,I),0.5);
        q75(k,nI) = quantile(Res(1,I),0.75);
        qM(k,nI) = mean(Res(1,I));
    end
end

figure
subplot(1,3,1);hold on
title('25% quantile')
plot(1:5,q25(1,:),1:5,q25(2,:),1:5,q25(3,:),1:5,q25(4,:),1:5,q25(5,:),'LineWidth',3)
h = legend('0.001','0.02575','0.0505','0.07525','0.1');
title(h,'Specificity')
grid on

subplot(1,3,2);hold on
title('Mean')
plot(1:5,qM(1,:),1:5,qM(2,:),1:5,qM(3,:),1:5,qM(4,:),1:5,qM(5,:),'LineWidth',3)
grid on

subplot(1,3,3);hold on
title('75% quantile')
plot(1:5,q75(1,:),1:5,q75(2,:),1:5,q75(3,:),1:5,q75(4,:),1:5,q75(5,:),'LineWidth',3)
grid on







%% Test for log normality, ok
% 
% load Res_G1.mat
% PValue = repmat(r',[1 , nS]);
% PValue = reshape(PValue',[1,numel(PValue)]);
% 
% 
% 
%  %% Test for lognormal distribution
%  
% h = waitbar(0,'Running KS test G1');
% 
% Test = zeros(length(R_Risk),length(r));
% 
% for nG = 1:length(R_Risk)
%     
%     Res = R_Risk{nG};
%     for i = 1:nL
%         
%         if sum(Res(:,i)) ~= 0
%             
%             O = log(Res(:,i) + ones);
%             O = (O - mean(O))./std(O);
%             Test(nG,i) = double(kstest(O));
%             waitbar((nG*i)/(nL*length(R_Risk)),h)
%         
%         else
%             
%             Test(nG,i) = -1;
%         end
%     end
% end
% numel(Test)
% sum(sum(Test==-1))
% sum(sum(Test == 0))
% sum(sum(Test == 1))
% 
% [l,c] = find(Test == 1);
% V = [];
% G = [];
% for i = 1:length(l)
%     
%     Res = R_Risk{l(i)};
%     V = [V ; Res(:,c(i))];
%     G = [G ; i.*ones(length(Res(:,c(i))),1)];
% %     histogram(Res(:,c(i)))
% end
% 
% 
% %% G2
% load Res_G2.mat
% PValue = repmat(r',[1 , nS]);
% PValue = reshape(PValue',[1,numel(PValue)]);
% 
% 
% 
%  %% Test for lognormal distribution
%  
% h = waitbar(0,'Running KS test G2');
% 
% Test = zeros(length(R_Risk),length(r));
% 
% for nG = 1:length(R_Risk)
%     
%     Res = R_Risk{nG};
%     for i = 1:nL
%         
%         if sum(Res(:,i)) ~= 0
%             
%             O = log(Res(:,i) + ones);
%             O = (O - mean(O))./std(O);
%             Test(nG,i) = double(kstest(O));
%             waitbar((nG*i)/(nL*length(R_Risk)),h)
%         
%         else
%             
%             Test(nG,i) = -1;
%         end
%     end
% end
% numel(Test)
% sum(Test==-1)
% sum(Test == 0)
% sum(Test == 1)
% 
% [l,c] = find(Test == 1);
% V = [];
% G = [];
% for i = 1:length(l)
%     
%     Res = R_Risk{l(i)};
%     V = [V ; Res(:,c(i))];
%     G = [G ; i.*ones(length(Res(:,c(i))),1)];
% %     histogram(Res(:,c(i)))
% end

clc
clear all
load Res_G2_FixedInfected.mat
PValue = repmat(r',[1 , nS]);
PValue = reshape(PValue',[1,numel(PValue)]);

%% Test for lognormal distribution

M = zeros(length(R_Risk),length(r));
X1 = zeros(length(R_Risk),length(r));
X2 = zeros(length(R_Risk),length(r));
G = 10:30;

for nG = 1:length(R_Risk)
    
    Res = R_Risk{nG};
    for i = 1:length(r)
        X1(nG,i) = (r(i)*(1-p))/p;
        X2(nG,i) = G(nG);
        M(nG,i) = mean(mean(Res(:,PValue == r(i))));
    end
    
end

mdl = fitlm([X1(:) , X2(:)],M(:),'interactions')
input.data = mdl.Coefficients;
% we want a complete LaTex document
input.makeCompleteLatexDocument = 1;

% generate LaTex code
latex = latexTable(input);

% save LaTex code as file
fid=fopen('MyLatex.tex','w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');




latex=latexTable(input)


%% Display mean risk as a function of initial FEC and Nb of goats
clc
clear all
load Res_G1.mat
PValue = repmat(r',[1 , nS]);
PValue = reshape(PValue',[1,numel(PValue)]);

figure
subplot(1,2,1);
hold on
M = zeros(length(R_Risk),length(r));
m = zeros(1,length(r));
LegendInfo = [];
k = 1;
for nG = 1:4:length(R_Risk)
    Res = R_Risk{nG};
    for i = 1:length(r)
        m(i) = (r(i)*(1-p))/p;
        M(nG,i) = mean(mean(Res(:,PValue == r(i))));
    end
    plot(m,M(nG,:),'LineWidth',2)
    LegendInfo{k} = int2str(nG);
    k = k + 1;
end

grid on
title('First Flock')
xlabel('Mean Initial FEC')
ylabel('Mean Flock Risk')
h = legend(LegendInfo);
title(h,'Number of goats')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,...
    'YColor',[10 10 10]./255,'Color',[195 195 195]./255)

load Res_G2.mat
PValue = repmat(r',[1 , nS]);
PValue = reshape(PValue',[1,numel(PValue)]);

subplot(1,2,2);
hold on
M = zeros(length(R_Risk),length(r));
m = zeros(1,length(r));
LegendInfo = [];
k = 1;
for nG = 1:4:length(R_Risk)
    Res = R_Risk{nG};
    for i = 1:length(r)
        m(i) = (r(i)*(1-p))/p;
        M(nG,i) = mean(mean(Res(:,PValue == r(i))));
    end
    plot(m,M(nG,:),'LineWidth',2)
    LegendInfo{k} = int2str(nG);
    k = k + 1;
end

grid on
title('Second Flock')
xlabel('Mean Initial FEC')
ylabel('Mean Flock Risk')
set(gca,'FontSize',20,'FontName','TimeNewsRoman','GridColor',[83 81 84]./255,'XColor',[10 10 10]./255,...
    'YColor',[10 10 10]./255,'Color',[195 195 195]./255)