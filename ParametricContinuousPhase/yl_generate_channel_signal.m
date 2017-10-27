function EEG = yl_generate_channel_signal(n_channel,add_noise,sampling_rate,data_length)
%% Setting parameters
n = n_channel;             

freq1 = zeros(1,n);  % reference frequency
freq1_ini = 35;


for i = 1:n
    freq1(i) = freq1_ini + 0.01*(i-1);
end

noise_var = add_noise;


% sampling_rate : sample per second
% data_length   : sec


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

end 
