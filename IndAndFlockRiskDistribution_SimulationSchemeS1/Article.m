% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Run it first to have all parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G1
% load Data.mat
% load Meteo.mat
% load ZNumber.mat
% load Distri_Weights_Feces.mat
% 
% nGoats = 21;
% Groupe = 1;
% 
% if Groupe == 1
%     ZN = ZN_G1;
% else
%     ZN = ZN_G2;
% end
% 
% %% Keep only farms with more than nGoats
% load OPG_Maurice.mat
% OPG = OPG(~isnan(OPG(:,1)),:);
% U = unique(OPG(:,2));
% nb = zeros(1,numel(U));
% O = [];
% for i = 1:numel(U)
%     
%     nb(i) = sum(OPG(:,2)==U(i));
%     if nb(i) > nGoats
%         O = [O ; OPG(OPG(:,2)==U(i),:)];
%     end 
%     
% end
% 
% OPG_I = O;
% 
% 
% %% Compute the occurrence frequency for the first week
% H = ComputeProportion(Groupe,Data,[1 4]);
% Idx = find(H > 0);
% F_W1 = H(Idx);
% nq = length(F_W1);
% 
% %% Compute the occurrence frequency for the second week
% H = ComputeProportion(Groupe,Data,[5 8]);
% F_W2 = H(Idx);
% 
% 
% %% Compute the average time spend in the field
% T = ComputeTime(Data,1,1:8);
% 
% Week_Time = mean(sum(T)/2);
% 
% mu_wf = DistF.mu/5;
% sigma_wf = DistF.sigma/5;
% mu_np = 8.5/11;
% sigma_np = 3.6/11;
% 
% 
% nIt = 1000;
% load ParameterOPG.mat
% nc = 5;
% r = linspace(min(Par(1,:)),max(Par(1,:)),nc);
% p = linspace(min(Par(2,:)),max(Par(2,:)),nc);
% 
% a_Init = linspace(0.001,0.1,5);
% k_s = 1;
% 
% nL = length(r)*length(p)*length(a_Init)*nIt;
% 
% 
% load Res_G1.mat
% load Stat_G1.mat
% load ParameterOPG.mat
% nc = 5;
% nIt = 1000;
% r = linspace(min(Par(1,:)),max(Par(1,:)),nc);
% p = linspace(min(Par(2,:)),max(Par(2,:)),nc);
% a_Init = linspace(0.001,0.1,5);
% k_s = 1;
% nL = length(r)*length(p)*length(a_Init)*nIt;
% 
% 
% PValue = zeros(3,nL);
% k_s = 1;
% for i = 1:length(r)
%     for j = 1:length(p)
%         for k = 1:length(a_Init)
%             for nS = 1:nIt
%                 PValue(:,k_s) = [r(i) , p(j) , a_Init(k)]';
%                 k_s = k_s + 1;
%             end
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% Examples of neg. binomial distribution for the schematic
% x = 0:80;
% plot(x,pdf('nbin',x,r(1),p(1)),'k','LineWidth',3)
% grid on
% xlabel('FEC')
% ylabel('Probability')
% set(gca,'FontSize',18,'FontName','TimeNewsRoman')
% 
% x = 0:10000;
% plot(x,pdf('nbin',x,r(end),p(end)),'k','LineWidth',3)
% grid on
% xlabel('FEC')
% ylabel('Probability')
% set(gca,'FontSize',18,'FontName','TimeNewsRoman')

%% Example of risk distribution
figure
j = 2;
z = 5;
J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
V = [];
G = [];
U = [];
for i = 1:length(a_Init)
    
    I = intersect(find(Specificity(1,:) == a_Init(i)),J);
    R = Res(:,I);
    R = R(:);
%     R = mean(((R - mean(R).*ones)/std(R)).^3);
    %             R = log(R);
    %             R = (R - mean(R))/std(R);
    V = [V ; R];
    G = [G ; round(a_Init(i),2).*ones(size(R))];
    U = [U , round(a_Init(i),2)];
end
subplot(1,2,1);hold on
bh = boxplot(V,G,'Symbol','w');
for i = 1:numel(bh)
    set(bh(i),'Linewidth',2.5')
end
grid on
q1 = quantile(V,0.25);
q3 = quantile(V,0.75);
w1 = q3 + 1.5*(q3 - q1);
ylim([(w2 - 0.1*w2) (w1 + 0.1*w1)])
title(['E(FEC) = ' num2str(round((r(j)*(1-p(z)))/p(z)))],'FontSize',12,'FontWeight','normal')
xlabel('Specificity')
ylabel('Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')

j = 5;
z = 2;
J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
V = [];
G = [];
U = [];
for i = 1:length(a_Init)
    
    I = intersect(find(Specificity(1,:) == a_Init(i)),J);
    R = Res(:,I);
    R = R(:);
    %             R = log(R);
    %             R = (R - mean(R))/std(R);
    V = [V ; R];
    G = [G ; round(a_Init(i),2).*ones(size(R))];
    U = [U , round(a_Init(i),2)];
end
subplot(1,2,2);hold on
bh = boxplot(V,G,'Symbol','w');
for i = 1:numel(bh)
    set(bh(i),'Linewidth',2.5')
end
grid on
q1 = quantile(V,0.25);
q3 = quantile(V,0.75);
w1 = q3 + 1.5*(q3 - q1);
ylim([(w2 - 0.1*w2) (w1 + 0.1*w1)])
title(['E(FEC) = ' num2str(round((r(j)*(1-p(z)))/p(z)))],'FontSize',12,'FontWeight','normal')
xlabel('Specificity')
ylabel('Risk')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')


%% Idem mais avec cdf
figure
j = 2;
z = 5;
J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
for i = 1:length(a_Init)
    
    I = intersect(find(Specificity(1,:) == a_Init(i)),J);
    R = Res(:,I);
    R = R(:);
%     R = mean(((R - mean(R).*ones)/std(R)).^3);
    %             R = log(R);
    %             R = (R - mean(R))/std(R);
    
    subplot(1,2,1);hold on
    [f , x] = ecdf(R);
    plot(log(x),f,'LineWidth',2)
    LI{i,1} = num2str(round(a_Init(i),2));
end
legend(LI)
grid on
xlim([0 14])
title(['E(FEC) = ' num2str(round((r(j)*(1-p(z)))/p(z)))],'FontSize',12,'FontWeight','normal')
xlabel('Risk (log scale)')
ylabel('Frequency')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')




j = 5;
z = 2;
J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
for i = 1:length(a_Init)
    
    I = intersect(find(Specificity(1,:) == a_Init(i)),J);
    R = Res(:,I);
    R = R(:);
%     R = mean(((R - mean(R).*ones)/std(R)).^3);
    %             R = log(R);
    %             R = (R - mean(R))/std(R);
    
    subplot(1,2,2);hold on
    [f , x] = ecdf(R);
    plot(log(x),f,'LineWidth',2)
    LI{i,1} = num2str(round(a_Init(i),2));
end
legend(LI)
grid on

title(['E(FEC) = ' num2str(round((r(j)*(1-p(z)))/p(z)))],'FontSize',12,'FontWeight','normal')
xlabel('Risk (log scale)')
ylabel('Frequency')
set(gca,'FontSize',16,'FontName','TimeNewsRoman')



%% risk distribution as a function of mean FEC and Specificity
figure
Here = pwd;
There = 'C:\Users\mtbonneau\Documents\Article\2017\StageMehdy';

T = zeros(1,25);
k = 1;
for j = 1:length(r)
    for z = 1:length(p)
        T(k) = round((r(j)*(1-p(z)))/p(z));
        k = k + 1;
    end
end

[T , Ind] = sort(T);

k = 1;
% Test = zeros(length(T),numel(a_Init)*(numel(a_Init)-1));
Test = zeros(length(T),length(r));
Test2 = zeros(length(T),length(r));
for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = [];
        G = [];
        U = [];
        for i = 1:length(a_Init)
            
            I = intersect(find(Specificity(1,:) == a_Init(i)),J);
            R = Res(:,I);
            R = R(:);
%             R = mean(((R - mean(R).*ones)./std(R)).^3);
            R = log(R);
%             R = (R - mean(R))/std(R);
            V = [V ; R];
            G = [G ; round(a_Init(i),2).*ones(size(R))];
            U = [U , round(a_Init(i),2)];
        end
        
%         V = log(V);
%         V = (V - mean(V))/std(V);
        %% kstest
%         for i = 1:length(a_Init)
%             
%             O = V(G == U(i));
%             Test(k,i) = mean(O);
%         end
        o = 1;
        for i = 1:length(a_Init)
            for ii = 1:length(a_Init)
                if ~isequal(i,ii)
                    O1 = V(G == U(i));
                    O2 = V(G == U(ii));
                    O1 = (O1 - mean(O1))/std(O1);
                    O2 = (O2 - mean(O2))/std(O2);
                    Test(k,o) = ttest2(O1,O2,'Alpha',0.05);
                    o = o + 1;
                end
            end
        end
        
        
        %% Try to test neg binomial, no! But log normal, yep
%         for i = 1:length(a_Init)     
%             n = 100;
%             pd = fitdist((V(G == U(i))),'lognormal');
%             E = random(pd,[1 length(V)]);
%             T0 = zeros(1,n);
%             for ii = 1:n
%                 s1 = randperm(numel(V));
%                 s2 = randperm(numel(V));
%                 T0(ii) = kstest2(E(s1(1:1000)),V(s1(1:1000)));
%             end
%             Test(k,i) = 100*(sum(T0 == true)/n)
%         end
        
        %% Try to test lognormality of the data -> no !
%         V = log(V);
%         V = (V - mean(V))/std(V);
        %% kstest
%         for i = 1:length(a_Init)
                 
%             O1 = V(G == U(i));
%                 Test(k,i) = chi2gof(V(G == U(i)));
%             [a b c d] = kstest(V(G == U(i)));
%             Test(k,i) = c;
%             Test2(k,i) = d;
%         end
%         o = 1;
%         for i = 1:length(a_Init)
%             for ii = 1:length(a_Init)
%                 if ~isequal(i,ii)
%                     
%                     Test(k,o) = kstest2(V(G == U(i)),V(G == U(ii)));
%                     o = o + 1;
%                 end
%             end
%         end
               
%         V = log(V);
%         subplot(5,5,k);hold on
% %         figure;hold on
%         bh = boxplot(V,G,'Symbol','w');
%         for i = 1:numel(bh)
%             set(bh(i),'Linewidth',2.5')
%         end
%         grid on
%         q1 = quantile(V,0.25);
%         q3 = quantile(V,0.75);
%         w1 = q3 + 1.5*(q3 - q1);
%         w2 = q1 - 1.5*(q3 - q1);
%         ylim([0 11])
%         ylim([(w2 - 0.1*w2) (w1 + 0.1*w1)])
% %         if k >= 21
% %             xlabel('Specificity')
% %         end
% %         if ~isempty(intersect(k,[1 6 11 16 21]))
% %             ylabel('Risk')
% %         end
% %         text(a_Init(2),w,['mean(OPG) = ' num2str(round((r(j)*(1-p(z)))/p(z)))])
%         title(['E(FEC) = ' num2str(round((r(j)*(1-p(z)))/p(z)))],'FontSize',12,'FontWeight','normal')
% title(int2str(k))
%         set(gca,'FontSize',16,'FontName','TimeNewsRoman')
%         cd(There)
%             print('-dpng',['FEC_Spec_' int2str(find(Ind ==k)) '.png'])
%         cd(Here)
%         close all
        k = k + 1;
    end
end
% 100*(sum(Test==true,2)./(k-1))
Test
Test2



%%% Mean risk as a function of FEC parameters
k = 1;
M = zeros(length(T),2);

for j = 1:length(r)
    for z = 1:length(p)
        
        J = intersect(find(PValue(1,:) == r(j)),find(PValue(2,:) == p(z)));
        V = [];
        G = [];
        U = [];
        for i = 1:length(a_Init)
            
            I = intersect(find(Specificity(1,:) == a_Init(i)),J);
            R = Res(:,I);
            R = R(:);
            V = [V ; R];
            G = [G ; round(a_Init(i),2).*ones(size(R))];
            U = [U , round(a_Init(i),2)];
        end

%         M(k,1) = quantile(V,0.95);
        M(k,1) = var(V);
        M(k,2) = ((r(j)*(1-p(z)))/p(z));

        k = k + 1;
    end
end
plot(M(:,2),M(:,1),'o')
mdl = fitlm(M(:,2),M(:,1))


%%% risk distribution at the flock's scale*


%% It's not exactly that. We should do one test per simulated flock. Each time you look if the distribution is neg bino or any other thing. 
%% Because in this configuration, within one boxplot we have all the individuals coming from different flocks. 
%% If we want to aggregate all, it'll be the estimated average risk value for each flock

V = [];
G = [];
k = [];
for i = 1:length(r)
    for j = 1:length(p)
        
        I = intersect(find(PValue(1,:) == r(i)),find(PValue(2,:) == p(j)));
        R = var(Res(:,I),1)';
%         R = R(:);
        V = [V ; R];
        G = [G ; (r(i)*(1-p(j)))/p(j).*ones(size(R))];
        k = [k , (r(i)*(1-p(j)))/p(j)];
    end
end

figure
bh = boxplot(log(V),round(G),'Symbol','w','Positions',sort(log(k)),'Labels',strvcat(num2str(sort(log(k')))),'LabelOrientation','inline','PlotStyle','traditional');
for i = 1:numel(bh)
set(bh(i),'Linewidth',1.5')
end
grid on
% ylim([0 17000])
ylim([0 10])
xlabel('LOG(mean OPG)')
ylabel('Log(Risk)')
set(gca,'FontSize',14,'FontName','TimeNewsRoman')

U = unique(G);
M = zeros(2,length(U));
for i = 1:length(U)
%     M(i) = mean(V(G == U(i)));
    pd = fitdist(V(G == U(i)),'lognormal');
    M(:,i) = [pd.mu ; pd.sigma];
end
plot(U,M)

figure
bh = boxplot(log(V),round(G),'Symbol','w');
for i = 1:numel(bh)
set(bh(i),'Linewidth',2.5')
end
grid on
ylim([0 10])
xlabel('mean OPG')
ylabel('Risk')
set(gca,'FontSize',14,'FontName','TimeNewsRoman')


bh = boxplot(V,G);
for i = 1:numel(bh)
set(bh(i),'Linewidth',2.5')
end
grid on
xlabel('mean OPG')
ylabel('Risk')
set(gca,'FontSize',18,'FontName','TimeNewsRoman')


