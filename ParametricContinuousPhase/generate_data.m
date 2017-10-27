function data  =  generate_data(chanal_num , sample_rate , data_length , scale_noise , SNR ,W_cfg , k_config)

%% chanal_num = number of chanal that we want to generate *** ex=32


%% sample rate = frequence of generating data *** ex=1000


%%%% data length = seconds of data that we want have synchrony in a
%%%% specific band , current data is 3X more than data length and our
%%%% synchrony varies by W_cfg and K_cfg  **** ex = 16


%%%% scale_noise = use this option to change power of noise *** ex = .5


%%%% SNR = signal to noise ratio *** ex = 250



%%%% W_cfg is an array  with 4 option that specifies centerfrequency and
%%%% detunitng  and a distribution for generating other chanal delta W

%%%% W_cfg.W is a vector with size 1*3 specifies ceterfrequncy for each part
%%%%of 3 part **** ex = [25 50 80]

%%%% W_cfg.detuning is a vector with size 1*3 specifies detuning of these
%%%% frequnecies **** ex = [3 5 1]

%%%% W_cfg.mean we need a distribution to generate new dutuing for other
%%%% chanal (exept 2 first chaal) , for this we use a guaussian 
%%%% distribution that this porperty specifies mean of it *** ex =0

%%%% W_cfg.sigma this property specifies sigma of gaussian distribtion for 
%%%% other chanals (except 2 first of them)  ex = 1;



%%%% K_cfg is an array  with 3 option that specifies couupling rength 
%%%% length and a distributio for other chanlas (except first chanal) 

%%%% K_cfg.K is a vector with size 1*3 specifies coupling strngth for 
%%%% each part of 3 part **** ex = [3 .2 1]

%%%% K_cfg.mean we need a distribution to generate new coupling strength 
%%%% for other chanal (exept first chaal), for this we use a guaussian 
%%%% distribution that this porperty specifies mean of it *** ex =0

%%%% K_cfg.sigma this property specifies sigma of gaussian distribtion for 
%%%% other chanals (excet first of them)  ex = 1;

number_of_chanals=chanal_num; 
% Randomize initial phases
initial_phase = randn(number_of_chanals,1)*2;
tim_sec=data_length;
dt = 1/sample_rate; % step size (here 1ms)
phases = zeros(number_of_chanals,(3*tim_sec)./dt);
phases(1:number_of_chanals,1)= initial_phase(1:number_of_chanals,1);

for ind=1:number_of_chanals
    noiseterm(ind,:)=powernoise(1,((3*tim_sec)./dt)+1,'normalize').*scale_noise;
end

%%%% making data for 3X more
KK = [];
WW = [];
for iii=1:3
    %%%%% gaussian distribuation
    delta_K = randn(number_of_chanals-1,1) * sqrt(k_config.sigma) + k_config.mean;
    delta_W = randn(number_of_chanals-1,1) * sqrt(W_cfg.sigma) + W_cfg.mean;
    KK(:,iii) = [k_config.K(iii); k_config.K(iii)+delta_K];
    
    WW(:,iii)  = [ W_cfg.W(iii)  ;W_cfg.W(iii)-W_cfg.detuning(iii)-delta_W ];
end

%%%%% smoothig changes over time
gus_win = gausswin(5000); % <-- this value determines the width of the smoothing window
gus_win = gus_win/sum(gus_win);


[row_num, col_num] = size(KK);

for ii=1:row_num
    temp_k = [];
    temp_w = [];
    for jj=1:col_num
        temp_k = [temp_k ones(1,data_length*sample_rate).*KK(ii,jj)];
        temp_w = [temp_w ones(1,data_length*sample_rate).*WW(ii,jj)];
    end
    
    K(ii,:) = conv(temp_k,gus_win,'same');
    W(ii,:) = conv(temp_w,gus_win,'same');
end


W=W.*(2*pi); %radians per sec
K=K.*(2*pi); %scaling of coupling


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for time=2:(3*tim_sec)/dt %%% SIMULATION of the phase-oscillators   %%%%%%%%%%
%     for ind=1:number_of_chanals
%         phases(ind,time)= phases(ind,time-1) + (dt*(W(ind,time))+ sum(-dt.*K(:,time).* (sin((phases(ind,time-1) -phases(:,time-1)))) ) )+  noiseterm(ind,time-1) ;
%     end
% end
phases(:,1) = 0.1*randn(number_of_chanals,1);
for time=2:(3*tim_sec)/dt %%% SIMULATION of the phase-oscillators   %%%%%%%%%%
     ind = 1;
     phases(ind,time) = phases(ind,time-1) + dt *W(ind,time-1) +  noiseterm(ind,time-1) ;
     for ind=2:number_of_chanals
         phases(ind,time)= phases(ind,time-1) + dt * W(ind,time-1)  + dt * K(ind,time-1)* sin(phases(1,time-1)-phases(ind,time-1)) +  noiseterm(ind,time-1) ;
     end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph=  (mod(phases',2*pi))';
xx=(exp(1i*(ph(:,1:1:end)))');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  NOISE  %%%%%%%
noiseterm=randn(number_of_chanals,3*tim_sec*sample_rate);%.*((3*tim_sec*sample_rate)/(2*sqrt(3*tim_sec*sample_rate)));  % white noise properly scaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = ((real(xx)'))+((noiseterm).* 1./sqrt(SNR)); %creating composite signal


%%%% additional functions
function x = powernoise(alpha, N, varargin)
% Generate samples of power law noise. The power spectrum
% of the signal scales as f^(-alpha).
% Useage:
%  x = powernoise(alpha, N)
%  x = powernoise(alpha, N, 'option1', 'option2', ...)
% Inputs:
%  alpha - power law scaling exponent
%  N     - number of samples to generate
% Output:
%  x     - N x 1 vector of power law samples
% With no option strings specified, the power spectrum is
% deterministic, and the phases are uniformly distributed in the range
% -pi to +pi. The power law extends all the way down to 0Hz (DC)
% component. By specifying the 'randpower' option string however, the
% power spectrum will be stochastic with Chi-square distribution. The
% 'normalize' option string forces scaling of the output to the range
% [-1, 1], consequently the power law will not necessarily extend
% right down to 0Hz.
% (cc) Max Little, 2008. This software is licensed under the
% Attribution-Share Alike 2.5 Generic Creative Commons license:
% http://creativecommons.org/licenses/by-sa/2.5/
% If you use this work, please cite:
% Little MA et al. (2007), "Exploiting nonlinear recurrence and fractal
% scaling properties for voice disorder detection", Biomed Eng Online, 6:23
% As of 20080323 markup
% If you use this work, consider saying hi on comp.dsp
% Dale B. Dalrymple

opt_randpow = false;
opt_normal = false;
for j = 1:(nargin-2)
    switch varargin{j}
        case 'normalize', opt_normal = true;
        case 'randpower', opt_randpow = true;
    end
end
N2 = floor(N/2)-1;
f = (2:(N2+1))';
A2 = 1./(f.^(alpha/2));
if (~opt_randpow)
    p2 = (rand(N2,1)-0.5)*2*pi;
    d2 = A2.*exp(i*p2);
else
    % 20080323
    p2 = randn(N2,1) + i * randn(N2,1);
    d2 = A2.*p2;
end
d = [1; d2; 1/((N2+2)^alpha); flipud(conj(d2))];
x = real(ifft(d));
if (opt_normal)
    x = ((x - min(x))/(max(x) - min(x)) - 0.5) * 2;
end

function r =  circ_dist(x,y)
%
% r = circ_dist(alpha, beta)
%   Pairwise difference x_i-y_i around the circle computed efficiently.
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable or one single angle
%   Output:
%     r       matrix with differences
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
% PHB 3/19/2009
% Circular Statistics Toolbox for Matlab
% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
    error('Input dimensions do not match.')
end
r = angle(exp(1i*x)./exp(1i*y));