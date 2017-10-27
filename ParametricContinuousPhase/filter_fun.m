function [real_out,imag_out] = filter_fun(mode,data,gab_filt,center_frq,sample_rate)
if mode==1
    [m,n] = size(data);
    for i=1:m
        data_temp = data(i,:);
        output = filter(gab_filt,1,data_temp);
        real_out(i,:) = real(output);
        imag_out(i,:) = imag(output);
    end
else
    [m,n]     = size(data);
    dt        = (0:n-1)/sample_rate;
    data_temp = data(1,:).*exp(sqrt(-1)*2*pi*center_frq*dt);
    output = filter(gab_filt,1,data_temp);
    real_out(1,:) = real(output);
    imag_out(1,:) = imag(output);
    
end

end