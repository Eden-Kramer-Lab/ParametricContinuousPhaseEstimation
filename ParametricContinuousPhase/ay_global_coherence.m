function COH = ay_global_coherence(EEG,Ls,Wnd,fs_a,fs_b)
% global coherence
[m,n]=size(EEG);
for i=1:m
    XF{i}=ay_fft(1,EEG(i,:),Ls,Wnd);
    i
end
% we have multiple measure of FFT per Wnd - number: Wnd/Ls
[m,z] = size(XF{1});
fs    = length(XF{1}{1});
COH   = zeros(m,fs_a-fs_b+1);
for i = 1:m
    i
    %% generate covariance matrix
    for j = fs_a:fs_b
        cov = zeros(size(EEG,1),size(EEG,1));
    
        for p=1:size(EEG,1)
            for q=1:size(EEG,1)
                Fa = XF{p}(i,:);
                Fb = XF{q}(i,:);
                for t=1:z
                     % covariance matrix
                     Ta(t) = Fa{t}(j);
                     Tb(t) = Fb{t}(j);
                end
                rTa = real(Ta)-mean(real(Ta));
                iTa = imag(Ta)-mean(imag(Ta));
                Ta = rTa + iTa*sqrt(-1);
                
                rTb = real(Tb)-mean(real(Tb));
                iTb = imag(Tb)-mean(imag(Tb));
                Tb  = rTb - iTb*sqrt(-1);
                               
                cov(p,q) = sum(Ta.*Tb);
            end
        end
        eigs    = svd(cov);
        COH(i,j-fs_a+1)= eigs(1)/sum(eigs);
    end
end