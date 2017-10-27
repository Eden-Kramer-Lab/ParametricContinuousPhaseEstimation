%% What to do
% Generate EEG data with larger noise and more sample data
% Subset the first 500 global coherence
% Spectrum image
% Time-domain signals

% Histogram of the whole data
% Subset the coherent part and non-coherent part

% Generate EEG 
% additive noise is 40
% fully coherent with different phase
% EEG=generate_channel_signal();
% 
% %% Calculate coherence
% % sample size of specified frequency is 200
% % Subset half og it
% COH = ay_global_coherence(EEG,1000,8000);
% COH = COH(:,1:500);

%% Load data
load('subset_COH.mat')
%% Visualization
p1 = figure()
imagesc(COH)
colorbar
title('Global Coherence');
xlabel('Frequency');
ylabel('Time');
%saveas(p1,sprintf('COH.png'));


p2 = figure()
subplot(2,1,1)
plot(COH(:,36));
title('Coherent Signals');
xlabel('Time');
ylabel('Global Coherence')
grid minor

subplot(2,1,2)
plot(COH(:,100));
title('Non-Coherent Signals');
xlabel('Time');
ylabel('Global Coherence')
grid minor
%saveas(p2,sprintf('Time-domain_Non_Coherent-Coherent.png'));
% 
% p3 = figure()
% n = 8;
% for i = 1:n
% subplot(2,n/2,i)
% specgram(EEG(i,:));
% title(['Channel ' num2str(i)]);
% end
% %saveas(p4,sprintf('Spectrum.png'));

%% Store coherent signals in Yn and store all signals in X
X2 = zeros(1,length(COH(:,36)));
for i = 1:length(COH(:,36))
  X2(i) = COH(i,36);
end

X = zeros(1,numel(COH));
for k = 1:numel(COH)
  X(k) = COH(k);
end

%% Subset the non-coherent data X1 and the coherent data X2 
% Check histogram of the whole data 
% Find the threshold T = 0.618
% Split the data from threshold
% X1 is the non-coherent data
p4 = figure()
histogram(X);
title(['General Distribution of EEG']);
xlabel('Global Coherence');
ylabel('EEG data')
grid minor
%saveas(p4,sprintf['Hist_whole.png'])

%%
X11 =  COH;
X11(:,36) = [];
X1 = zeros(1,numel(X11));
for k = 1:numel(X11)
  X1(k) = COH(k);
end

p5 = figure()
subplot(1,2,1)
hist(X1,0:0.01:1);
h = hist(X1,0:0.01:1);
hold on
plot(0:0.01:1,h,'r','LineWidth',1.2);
title('Non Coherent Data with PDF');
xlabel('Gobal Coherence');
xlim([0.25 0.6]);
grid minor

subplot(1,2,2)
hist(X2,0:0.02:1);
h = hist(X2,0:0.02:1);
hold on
plot(0:0.02:1,h,'r','LineWidth',1.2);
title('Coherent Data with PDF');
xlabel('Gobal Coherence');
xlim([0.45 0.9])
ylim([0 41]);
grid minor
%saveas(p5,sprintf('Histfit_Coherent_Non_Coherent.png'));

% Store coherent in COH to an array
Yn = zeros(1,length(COH(:,36)));
for i = 1:length(COH(:,36))
  Yn(i) = COH(i,36);
end
% figure()
% plot(Yn)

%% Fitting distribution
% Try Gamma, Beta and Log Normal on non-coherent data and coheret data
% Assume d is the shifter and find the upperbound of shifter
d1_max = min(X1);
d2_max = min(X2);
l1= 0:0.01:0.25;
l2 = (-0.50):0.01:0.61;
q1 = length(l1);
q2 = length(l2);
OrgX1 = X1;
OrgX2 = X2;

%% Non-coherent part 

%% Gamma distribution
h1_gamma = zeros(length(l1),1);
p1_gamma = zeros(length(l1),1);
max_y1_diff_gamma = zeros(length(l1),1);
y1_emp_gamma  = zeros(length(l1),1);
y1_hyp_gamma = zeros(length(l1),1);

for i = 1:q1
d = l1(i);
X1 = OrgX1-d;

%% Fitting data to the Gamma model 
p =gamfit(X1)
%% Compare the cdf of {global coherence} and their estimated distribution
m1 = figure(i)
y = cdfplot(X1)
hold on
X0 = y.XData;
y_emp = y.YData;
y1 = gamcdf(X0,p(1),p(2));
m11 = plot(X0,y1)


X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y1_diff_gamma(i) = max(y_diff);
index = find(y_diff == max_y1_diff_gamma(i));
y1_ks = y1(index);
y_emp_ks = y_emp(index);
m12=plot([X0(index) X0(index)],[y1_ks y_emp_ks],'LineWidth',1.5);
plot([X0(index) X0(index)],[y1_ks 0],'r--')
plot([X0(index) 0],[y1_ks y1_ks],'r--')
plot([X0(index) 0],[y_emp_ks y_emp_ks],'r--')
%title(['The shifter is ' num2str(l1(i)) ' with K-S distance  ' num2str(max_y1_diff_gamma(i))]);
title(['Fitting Gamma Distribution with shifter '  num2str(l1(i)) ' with K-S distance  is ' num2str(max_y1_diff_gamma(i))]);
hold off
legend([y m11],{'Empirical cdf','Hypothesis Gamma cdf'});
xlabel('Global Coherence of Non-coherent Data');
ylabel('Culmulative Distribution Function');
saveas(m1,sprintf('Gamma_CDF%d.png',i));

%% pdf of {global coherence} and their estimated distribution
m2=figure(i+27)
h = hist(X1,0:0.01:1);
hold on
h = h/sum(h);
plot(0:0.01:1,h);
hold on
pdf = gampdf(0:0.01:1,p(1),p(2));
pdf = pdf/100;
plot(0:0.01:1,pdf);
legend('Empirical pdf','Hypothesis Gamma pdf');
hold off;
xlabel('Global Coherence of Non-coherent Data');
ylabel('Probability Distribution Function');
title(['Fitting Gamma Distribution with shifter '  num2str(l1(i))]);
grid minor
saveas(m2,sprintf('Gamma_PDF%d.png',i));

% Plot additional figure with KS-distance
m3 = figure(i+54)
plot(y1,y_emp);
hold on
plot(y1,y1,'k--');
hold off
legend('Empirical','Reference Line');
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
title(['Fitting Gamma Distribution with shifter '  num2str(l1(i))]);
grid minor
saveas(m3,sprintf('Gamma_Comparison%d.png',i));

%% goodness of fit
pd_gamma = fitdist(X1','Gamma')
[h_gamma_chi,p_gamma_chi,st_gamma_chi] = chi2gof(X1,'CDF',pd_gamma)
% Store whether reject in the matrix
h1_gamma(i) = h_gamma_chi;
p1_gamma(i) = p_gamma_chi;
y1_emp_gamma  = y.YData;
y1_hyp_gamma = y1;
end
close all

%% Beta
h1_beta = zeros(length(l1),1);
p1_beta = zeros(length(l1),1);
max_y1_diff_beta = zeros(length(l1),1);
y1_emp_beta  = zeros(length(l1),1);
y1_hyp_beta = zeros(length(l1),1);

for i = 1:q1
d = l1(i);
X1 = OrgX1-d;
%% Fitting data to the Gamma model 
p =betafit(X1);
%% Compare the cdf of {global coherence} and their estimated distribution
 figure(i)
y = cdfplot(X1)
hold on
X0 = y.XData;
y1 = betacdf(X0,p(1),p(2));
y_emp = y.YData;
m2 = plot(X0,y1);

X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y1_diff_beta(i) = max(y_diff);

index = find(y_diff == max_y1_diff_beta(i));
y1_ks = y1(index);
y_emp_ks = y_emp(index);
m12=plot([X0(index) X0(index)],[y1_ks y_emp_ks],'LineWidth',1.5);
plot([X0(index) X0(index)],[y1_ks 0],'r--')
plot([X0(index) 0],[y1_ks y1_ks],'r--')
plot([X0(index) 0],[y_emp_ks y_emp_ks],'r--')
%title(['The max distance is ' num2str(max_y1_diff_beta(i))]);
title(['Fitting Beta Distribution with shifter '  num2str(l1(i)) ' with K-S distance  is ' num2str(max_y1_diff_beta(i))]);
legend([y m2],{'Empirical cdf','Hypothesis Beta cdf'});
xlabel('Global Coherence of Non-coherent Data');
ylabel('Culmulative Distribution Function');
hold off
%saveas(m1,sprintf('Beta_CDF%d.png',i));

%% pdf of {global coherence} and their estimated distribution
m2 = figure(i+27)
h = hist(X1,0:0.01:1);
hold on
h = h/sum(h);
plot(0:0.01:1,h);
pdf = betapdf(0:0.01:1,p(1),p(2));
pdf = pdf/100;
plot(0:0.01:1,pdf);
legend('Empirical pdf','Hypothesis Beta pdf');
xlabel('Global Coherence of Non-coherent Data');
ylabel('Probability Distribution Function');
title(['Fitting Beta Distribution with shifter '  num2str(l1(i))]);
grid minor
hold off
%saveas(m2,sprintf('Beta_PDF%d.png',i));

% Plot additional figure
m3 = figure(i+54)
plot(y1,y_emp);
hold on
plot(y1,y1,'k--');
hold off
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
legend('Empirical','Reference Line');
title(['Fitting Beta Distribution with shifter '  num2str(l1(i))]);
grid minor
%saveas(m3,sprintf('Beta_Comparison%d.png',i));

%% goodness of fit
pd_beta = fitdist(X1','Beta')
[h_beta_chi,p_beta_chi,st_beta_chi] = chi2gof(X1,'CDF',pd_beta)
% Store whether reject in the matrix
h1_beta(i) = h_beta_chi;
p1_beta(i) = p_beta_chi;
y1_emp_beta  = y.YData;
y1_hyp_beta = y1;
end
close all

%% Lognormal distribution
d = 0;
X1 = OrgX1-d;

%% Fitting data to the lognormal model 
p =lognfit(X1)
%% Compare the cdf of {global coherence} and their estimated distribution
figure(i)
y = cdfplot(X1)
hold on
X0 = y.XData;
y_emp = y.YData;
y1 = logncdf(X0,p(1),p(2));
m2= plot(X0,y1)

X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y1_diff_logn = max(y_diff);
index = find(y_diff == max_y1_diff_logn);
y1_ks = y1(index);
y_emp_ks = y_emp(index);
plot([X0(index) X0(index)],[y1_ks y_emp_ks],'LineWidth',1.5);
plot([X0(index) X0(index)],[y1_ks 0],'r--')
plot([X0(index) 0],[y1_ks y1_ks],'r--')
plot([X0(index) 0],[y_emp_ks y_emp_ks],'r--')
%title(['The max distance is ' num2str(max_y1_diff_beta(i))]);
title(['Fitting Lognomal Distribution with K-S distance  is ' num2str(max_y1_diff_logn)]);
hold off
legend([y m2],{'Empirical cdf','Hypothesis Beta cdf'});
xlabel('Global Coherence of Non-coherent Data');
ylabel('Culmulative Distribution Function');
hold off
saveas(m1,sprintf('Lognormal_CDF.png'));

%% pdf of {global coherence} and their estimated distribution
m2 = figure(i+27)
h = hist(X1,0:0.01:1);
hold on
h = h/sum(h);
plot(0:0.01:1,h);
hold on
pdf = lognpdf(0:0.01:1,p(1),p(2));
pdf = pdf/100;
plot(0:0.01:1,pdf);
legend('Empirical pdf','Hypothesis Lognormal pdf');
hold off;
xlabel('Global Coherence of Non-coherent Data');
ylabel('Probability Distribution Function');
title('Fitting Lognormal Distribution');
grid minor
hold off
saveas(m2,sprintf('Lognormal_PDF.png'));

% Plot additional figure
m3 = figure(i+54)
plot(y1,y_empirical);
hold on
plot(y1,y1,'k--');
legend('Empirical','Reference Line');
title('Fitting Lognormal Distribution');
grid minor
hold off
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
saveas(m3,sprintf('Lognormal_Comparison.png'));

%% goodness of fit
pd_logn = fitdist(X1','Lognormal')
[h_logn_chi,p_logn_chi,st_logn_chi] = chi2gof(X1,'CDF',pd_logn)
% Store whether reject in the matrix
h1_logn = h_logn_chi;
p1_logn =p_logn_chi;
y1_emp_logn  = y.YData;
y1_hyp_logn = y1;

close all
%% Check the max distance
gamma_d1 = min(max_y1_diff_gamma);
beta_d1 = min(max_y1_diff_beta);
logn_d1 = min(max_y1_diff_logn);

%% Coherent Part
X0 = [];
y_empirical = [];
y1 = [];

%% Gamma distribution
h2_gamma = zeros(length(l2),1);
p2_gamma = zeros(length(l2),1);
max_y2_diff_gamma = zeros(length(l2),1);
y2_emp_gamma  = zeros(length(l2),1);
y2_hyp_gamma = zeros(length(l2),1);

for i = 1:q2
d = l2(i);
X2 = OrgX2-d;

%% Fitting data to the Gamma model 
p =gamfit(X2)
%% Compare the cdf of {global coherence} and their estimated distribution
m1 = figure(i)
y = cdfplot(X2)
hold on
X0 = y.XData;
y_emp = y.YData;
y1 = gamcdf(X0,p(1),p(2));
plot(X0,y1)
legend('Empirical cdf','Hypothesis Gamma cdf');

X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y2_diff_gamma(i) = max(y_diff);
title(['The max distance is ' num2str(max_y2_diff_gamma(i))]);
hold off
%saveas(m1,sprintf('Coherent_Gamma_CDF%d.png',i));

%% pdf of {global coherence} and their estimated distribution
m2=figure(i+27)
h = hist(X2,0:0.01:1);
hold on
h = h/sum(h);
plot(0:0.01:1,h);
hold on
pdf = gampdf(0:0.01:1,p(1),p(2));
pdf = pdf/100;
plot(0:0.01:1,pdf);
legend('Empirical pdf','Hypothesis Gamma pdf');
hold off;
%saveas(m2,sprintf('Coherent_Gamma_PDF%d.png',i));

% Plot additional figure
m3 = figure(i+54)
plot(y1,y_emp);
hold on
plot(y1,y1,'k--');
hold off
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
%saveas(m3,sprintf('Coherent_Gamma_Comparison%d.png',i));

%% goodness of fit
pd_gamma = fitdist(X2','Gamma')
[h_gamma_chi,p_gamma_chi,st_gamma_chi] = chi2gof(X2,'CDF',pd_gamma)
% Store whether reject in the matrix
h2_gamma(i) = h_gamma_chi;
p2_gamma(i) = h_gamma_chi;
y2_emp_gamma  = y.YData;
y2_hyp_gamma = y1;
end

%% Beta
h2_beta = zeros(length(l2),1);
p2_beta = zeros(length(l2),1);
max_y2_diff_beta = zeros(length(l2),1);
y2_emp_beta  = zeros(length(l2),1);
y2_hyp_beta = zeros(length(l2),1);

for i = 1:q2
d = l2(i);
X2 = OrgX2-d;

%% Fitting data to the Gamma model 
p =betafit(X2);
%% Compare the cdf of {global coherence} and their estimated distribution
m1 = figure(i)
y = cdfplot(X2)
hold on
X0 = y.XData;
y1 = betacdf(X0,p(1),p(2));
y_emp = y.YData;
plot(X0,y1)
legend('Empirical cdf','Hypothesis Beta cdf');

X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y2_diff_beta(i) = max(y_diff);
title(['The max distance is ' num2str(max_y2_diff_beta(i))]);
hold off
% saveas(m1,sprintf('Coherent_Beta_CDF%d.png',i));

% %% pdf of {global coherence} and their estimated distribution
% m2 = figure(i+27)
% h = hist(X2,0:0.01:1);
% hold on
% h = h/sum(h);
% plot(0:0.01:1,h);
% hold on
% pdf = betapdf(0:0.01:1,p(1),p(2));
% pdf = pdf/100;
% plot(0:0.01:1,pdf);
% hold off
% saveas(m2,sprintf('Coherent_Beta_PDF%d.png',i));

% Plot additional figure
m3 = figure(i+54)
plot(y1,y_emp);
hold on
plot(y1,y1,'k--');
hold off
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
% saveas(m3,sprintf('Coherent_Beta_Comparison%d.png',i));

%% goodness of fit
pd_beta = fitdist(X2','Beta')
[h_beta_chi,p_beta_chi,st_beta_chi] = chi2gof(X2,'CDF',pd_beta)
% Store whether reject in the matrix
h2_beta(i) = h_beta_chi;
p2_beta(i) = p_beta_chi;
y2_emp_beta  = y.YData;
y2_hyp_beta = y1;
end

%% Lognormal distribution
d = 0;
X2 = OrgX2-d;

% Fitting data to the lognormal model 
p =lognfit(X2)

% Compare the cdf of {global coherence} and their estimated distribution
m1 = figure(i)
y = cdfplot(X2)
hold on
X0 = y.XData;
y_empirical = y.YData;
y1 = logncdf(X0,p(1),p(2));
plot(X0,y1)
legend('Empirical cdf','Hypothesis Lognormal cdf');

X0_min = min(X0);
X0_max = max(X0);
y_diff = abs(y1-y.YData);
max_y2_diff_logn = max(y_diff);
title(['The max distance is ' num2str(max_y2_diff_logn)]);
hold off
saveas(m1,sprintf('Coherent_Lognormal_CDF.png'));

%% pdf of {global coherence} and their estimated distribution
m2 = figure(i+27)
h = hist(X2,0:0.01:1);
hold on
h = h/sum(h);
plot(0:0.01:1,h);
hold on
pdf = lognpdf(0:0.01:1,p(1),p(2));
pdf = pdf/100;
plot(0:0.01:1,pdf);
hold off
saveas(m2,sprintf('Coherent_Lognormal_PDF.png'));

% Plot additional figure
m3 = figure(i+54)
plot(y1,y_empirical);
hold on
plot(y1,y1,'k--');
hold off
xlabel('Hypothesis CDF');
ylabel('Empirical CDF');
saveas(m3,sprintf('Coherent_Lognormal_Comparison.png'));

%% goodness of fit
pd_logn = fitdist(X2','Lognormal')
[h_logn_chi,p_logn_chi,st_logn_chi] = chi2gof(X2,'CDF',pd_logn)
% Store whether reject in the matrix
h2_logn = h_logn_chi;
p2_logn =p_logn_chi;
y2_emp_logn  = y.YData;
y2_hyp_logn = y1;


%% Check the max distance
gamma_d2 = min(max_y2_diff_gamma);
beta_d2 = min(max_y2_diff_beta);
logn_d2 = min(max_y2_diff_logn);