function  EEG = generate_channel_signal()
%% Setting parameters
n = 16;             

freq1 = zeros(1,n);  % reference frequency
freq2 = zeros(1,n);  % partially coherent by changing one freq
freq3 = zeros(1,n);  % completely non-coherent by big jump in freq

freq1_ini = 35;
freq2_ini = 80;
freq3_ini = 40;

for i = 1:n
    freq1(i) = freq1_ini + 0.01*(i-1);
    freq2(i) = freq2_ini + 0.01*(i-1);
    freq3(i) = freq3_ini + 20*(i-1);
end

noise_var = 200;
phase_noise = 3;
phase_noise5 = 8; % large phase noise

sampling_rate = 1000;   % sample per second
data_length   = 8000;    % sec
ratio = 0.75; % rate for partially coherent channels

%% Generate samples for this 300 sec
Phase = 2*pi*rand(n,1);

%% generate data
EEG = zeros(n,1*data_length*sampling_rate);

%% Fully coherent

for i = 1:n
    ind = 1:sampling_rate*data_length;
    dt  = 1/sampling_rate;
    EEG(i,ind)= cos(2*pi*freq1(i)*dt*ind+Phase(i))+sqrt(noise_var)*randn(size(ind));
end

% %% Add noise to the phase
% for i = 1:n
%     ind = (1:sampling_rate*data_length)+sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq1(i)*dt*ind+Phase(i)+sqrt(phase_noise)*randn(size(ind)))+sqrt(noise_var)*randn(size(ind));
% end
% 
% %% Partially coherent - different frequency
% for i = 1:ceil(ratio*n)
%     ind = (1:sampling_rate*data_length)+2*sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq1(i)*dt*ind+Phase(i))+sqrt(noise_var)*randn(size(ind));   
% end
% 
% for i = (ceil(ratio*n)+1):n
%     ind = (1:sampling_rate*data_length)+2*sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq3(i)*dt*ind+Phase(i))+sqrt(noise_var)*randn(size(ind));   
% end
% 
% %% Partially coherent - adding noise to phase
% for i = 1:ceil(ratio*n)
%     ind = (1:sampling_rate*data_length)+3*sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq1(i)*dt*ind+Phase(i))+sqrt(noise_var)*randn(size(ind));
% end
% 
% for i = (ceil(ratio*n)+1):n
%     ind = (1:sampling_rate*data_length)+3*sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq1(i)*dt*ind+Phase(i)+sqrt(phase_noise)*randn(size(ind)))+sqrt(noise_var)*randn(size(ind));
% end
% 
% %% Completely non-coherent - large noise
% for i = 1:n
%     ind = (1:sampling_rate*data_length)+4*sampling_rate*data_length;
%     dt = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq1(i)*dt*ind+Phase(i)+sqrt(phase_noise5)*randn(size(ind)))+sqrt(noise_var)*randn(size(ind));
% end    
% 
% %% Completely non-coherent - larger jump in frequency
% for i = 1:n
%     ind = (1:sampling_rate*data_length)+5*sampling_rate*data_length;
%     dt  = 1/sampling_rate;
%     EEG(i,ind)=cos(2*pi*freq3(i)*dt*ind+Phase(i))+sqrt(noise_var)*randn(size(ind));
% end
end 
