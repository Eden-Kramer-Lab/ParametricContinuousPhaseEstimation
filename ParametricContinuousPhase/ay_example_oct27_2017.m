clc
clear 
close all

%% Generate data BY Dr code
% ch_no       = 32;
% sample_rate = 1000;
% data_length = 160;
% base_noise  = 0.01;
% EEG = ay_generate_sample_signal(ch_no,sample_rate,data_length,base_noise);

%%  Generate Data By Reza Code
W_cfg = [];
W_cfg.W = [7 19 55];
W_cfg.detuning = [0.05 0.15 0.05];
W_cfg.mean = 0;
W_cfg.sigma = 0.5;

k_config = [];
k_config.K = [1.5 0.05 0.5];
k_config.mean = 0;
k_config.sigma = 0.04;

chanal_num = 2;
sample_rate = 1000;
data_length = 10;
scale_noise = 0.09;
SNR = 2;
EEG = generate_data(chanal_num , sample_rate , data_length , scale_noise , SNR ,W_cfg , k_config);

%% Build Analytical Filter
disp off
[gabor,base_filt] = generate_gabor(19,5,sample_rate,2000);
% N   = 1000;        % FIR filter order
% Fp  = 3;          % 5 Hz passband-edge frequency
% Fs  = 1e3;        % 1 kHz sampling frequency
% Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
% Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation
% 
% base_filt = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
% %fvtool(base_filt,'Fs',Fs,'Color','White') % Visualize filter
% disp on
% 
%% Take One Channel, Build Param
ch_no = 2;
%[Re_temp,Im_temp]=filter_fun(EEG(1,:),gabor/1000);
Iter = 20;
mode = 1;
if mode==1
    [Re_temp,Im_temp]=filter_fun(1,EEG(ch_no,:),gabor);
else
    [Re_temp,Im_temp]=filter_fun(2,EEG(ch_no,:),base_filt,19,sample_rate);
end
[rXSmt,rSSmt,Param,rXPos,rSPos] = ay_fft_filter([Re_temp' Im_temp'],Iter);

%% Find Probability of Specific Phase
p_in_range = ay_phase(rXSmt,rSSmt,[-1 cos(pi*170/180)]);

%% Filtering Result
figure(1)
subplot(2,2,1)
plot(Re_temp); hold on;
plot(Im_temp); hold off;
axis tight
xlabel('Time Index')
ylabel('Real/Imag')
title('Baseband signal')

subplot(2,2,3)
plot(p_in_range);
xlabel('Time Index')
ylabel('Proabbility of 180 deg')

subplot(2,2,2)
K  = length(rXSmt);
Xm = zeros(K,1);
Xb = zeros(K,1);
for k=1:K
    Xm(k) = rXSmt{k}(1);
    Xb(k) = rSSmt{k}(1,1);
end
ay_plot_bound(1,(1:K),Xm',(Xm-2*sqrt(Xb))',(Xm+2*sqrt(Xb))');
xlabel('Time Index')
ylabel('Filtered Real')

subplot(2,2,4)
Xm = zeros(K,1);
Xb = zeros(K,1);
for k=1:K
    Xm(k) = rXSmt{k}(2);
    Xb(k) = rSSmt{k}(2,2);
end
ay_plot_bound(1,(1:K),Xm',(Xm-2*sqrt(Xb))',(Xm+2*sqrt(Xb))');
xlabel('Time Index')
ylabel('Filtered Imag')

%% Filter
figure(2)
subplot(1,2,1);
plot(base_filt,'LineWidth',2);
title('Low Pass Filter')
xlabel('Index')
subplot(1,2,2);
plot(real(gabor),'LineWidth',2);hold on;plot(imag(gabor),'LineWidth',2);hold off;
title('Equivalent Real & Imaginary Filters')
xlabel('Index')

%% Singal 
figure(3)
plot(EEG(1,:),'LineWidth',2)
hold on;
plot(11000:12000,EEG(1,11000:12000),'r','LineWidth',2)
hold off
xlabel('Time Index')
ylabel('Amplitude')
title('Signal')
axis tight



