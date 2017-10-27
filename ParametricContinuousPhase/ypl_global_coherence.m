function COH= ay_global_coherence_methods(method,EEG,Ls,Wnd,n_channel,Iter,fs_a,fs_b,Ns)
% global coherence
[m,n]=size(EEG);
for i=1:m
    XF{i} = ay_fft(2,EEG(i,:),Ls,Wnd);
end

[m,z] = size(XF{1});
fs    = length(XF{1}{1});
COH   = zeros(m,fs_b-fs_a+1);

%% filter each element in XF
parfor p = 1:n_channel
    for i=fs_a:fs_b
        Re_temp = [];
        Im_temp = [];
        for j = 1:m
            for k=1:z
                Re_temp = [Re_temp;real(XF{p}{j,k}(i))];
                Im_temp = [Im_temp;imag(XF{p}{j,k}(i))];
            end
        end
        disp(['Iteration ' num2str(p) ' ' num2str(i) ' ' num2str(j)])
        [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb] = yl_filter(Re_temp,Iter);
        % Real_x = rXPos;
        Param.Bk = 0;
        Real_x   = ay_state_sample(Ns,rXSmt,rSSmt,rXPos,rSPos,Param,zeros(length(Re_temp),1));
        
        [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb] = yl_filter(Im_temp,Iter);
        % Img_x = rXPos;
        Param.Bk = 0;
        Img_x = ay_state_sample(Ns,rXSmt,rSSmt,rXPos,rSPos,Param,zeros(length(Im_temp),1));
        
        for samp=1:Ns
            ind = 1;
            for j = 1:m
                for k=1:z
                    XF_filter{samp,p}{j,k}(i) = complex(Real_x(samp,ind),Img_x(samp,ind));
                    ind = ind + 1;
                end
            end
        end
    end
    
end

% we have multiple measure of FFT per Wnd - number: Wnd/Ls
for i = 1:m
    i
    %% generate covariance matrix
    for j = fs_a:fs_b
        cov = zeros(size(EEG,1),size(EEG,1));
        for p=1:size(EEG,1)
            for q=p:size(EEG,1)
                Fa = [];
                Fb = [];
                for samp=1:Ns
                    temp = cell2mat(XF_filter{samp,p}(i,:)');
                    Fa   = [Fa;temp(:,j)];
                    temp = cell2mat(XF_filter{samp,q}(i,:)');
                    Fb = [Fb;temp(:,j)];
                end
                % method 1:take one sample and center by the mean of phase
                if method==1
                    Ta = [];
                    Tb = [];
                    for t=1:length(Fa)
                         % covariance matrix
                         Ta(t) = Fa(t);  
                         Tb(t) = Fb(t); 
                    end
                    rTa = real(Ta)-mean(real(Ta));
                    iTa = imag(Ta)-mean(imag(Ta));
                    Ta = rTa + iTa*sqrt(-1);

                    rTb = real(Tb)-mean(real(Tb));
                    iTb = imag(Tb)-mean(imag(Tb));
                    Tb  = rTb - iTb*sqrt(-1);

                    cov(p,q) = sum(Ta.*Tb);
                    cov(q,p) = conj(cov(p,q));
                else

                    for t=1:length(Fa)
                         % method 2:draw Ns samples and normalized
                         % take the mean, go back to amp
                         Ta(t) = Fa(t)/abs(Fa(t));  
                         Tb(t) = Fb(t)/abs(Fb(t));
                    end
                    mTa = mean(Ta); pTa=angle(mTa);
                    mTb = mean(Tb); pTb=angle(mTb);
                    for t=1:length(Fa)
                         % method 2:draw Ns samples and normalized
                         % take the mean, go back to amp
                         Ta(t) = complex(abs(Fa(t))*cos(angle(Fa(t))-pTa),abs(Fa(t))*sin(angle(Fa(t))-pTa));
                         Tb(t) = complex(abs(Fb(t))*cos(angle(Fb(t))-pTb),-abs(Fb(t))*sin(angle(Fb(t))-pTb));
                    end
                    
                    cov(p,q) = sum(Ta.*Tb);
                    cov(q,p) = conj(cov(p,q));
                    
                end
                    
            end
        end
        eigs    = svd(cov);
        COH(i,j-fs_a+1)= eigs(1)/sum(eigs);
    end
end
