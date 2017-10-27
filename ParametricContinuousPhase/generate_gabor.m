function [out_filt,Y] = generate_gabor(center_fr,band_width,sample_rate,peak_ratio)

T  = -100:1/sample_rate:100;
sigma = 6/(2*pi*band_width);
Y  = normpdf(T,0,sigma);
out_filt = Y .* exp(sqrt(-1)*2*pi*center_fr.*T);
[~,left_ind]  = min(abs(Y(1:end/2)*peak_ratio-max(Y)));
[~,right_ind] = min(abs(Y(end/2:end)*peak_ratio-max(Y)));
out_filt  = out_filt(left_ind:right_ind+length(Y)/2);
Y  = Y(left_ind:right_ind+length(Y)/2);
end