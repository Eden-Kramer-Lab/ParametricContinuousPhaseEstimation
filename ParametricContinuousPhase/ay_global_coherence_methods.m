function [COH,COH_F,Extra] = ay_global_coherence_methods(method,EEG,Ls,Wnd,Iter,fs_a,fs_b,Ns)
%% Setting
% "method 1" is PNAS paper with Gamma Filter (no energy normalization, hamming window)
% "method 2" is PNAS preceded by FFT Filter (no energy normalization, no hamming window) with Gamma filter
% "method 3" is PLV based (energy normalization, no hamming window) followed by Gamma filter
% "method 4" is PLV based preceded by FFT filter(energy normalization, no hamming window) followed by Gamma filter


% global coherence calculation
n_channel=size(EEG,1);
XF = cell(n_channel,1);
for i=1:n_channel
    XF{i} = ay_fft(2,EEG(i,:),Ls,Wnd);
end

% filter each element in XF
[m,z] = size(XF{1});
COH   = zeros(m,fs_b-fs_a+1);

if method == 1 || method == 3  % PNAS Method & my PLV method
    for i = 1:m
      for j = fs_a:fs_b
            cov = zeros(size(EEG,1),size(EEG,1));
            for p=1:size(EEG,1)
                for q=p:size(EEG,1)
                    Fa = XF{p}(i,:);
                    Fb = XF{q}(i,:);
                    if method==1 || method==2
                        for t=1:length(Fa)
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
                    else
                    
                        for t=1:length(Fa)
                            % method 2:draw Ns samples and normalized
                            % take the mean, go back to amp
                            Ta(t) = Fa{t}(j)/abs(Fa{t}(j));  
                            Tb(t) = Fb{t}(j)/abs(Fb{t}(j));
                        end
                        mTa = mean(Ta); pTa=angle(mTa);
                        mTb = mean(Tb); pTb=angle(mTb);
                        for t=1:length(Fa)
                            % method 2:draw Ns samples and normalized
                            % take the mean, go back to amp
                            % Ta(t) = complex(abs(Fa{t}(j))*cos(angle(Fa{t}(j))-pTa),abs(Fa{t}(j))*sin(angle(Fa{t}(j))-pTa));
                            % Tb(t) = complex(abs(Fb{t}(j))*cos(angle(Fb{t}(j))-pTb),-abs(Fb{t}(j))*sin(angle(Fb{t}(j))-pTb));
                            Ta(t) = complex(cos(angle(Fa{t}(j))-pTa),sin(angle(Fa{t}(j))-pTa));
                            Tb(t) = complex(cos(angle(Fb{t}(j))-pTb),-sin(angle(Fb{t}(j))-pTb));
                        end
                    end
                    cov(p,q) = sum(Ta.*Tb);
                    cov(q,p) = conj(cov(p,q));
                end
            end
        eigs    = svd(cov);
        COH(i,j-fs_a+1)= eigs(1)/sum(eigs);
      end
    end
else
    % FFT Filter First
    for p = 1:n_channel
   % for p = 1:n_channel
        for i=fs_a:fs_b
            Re_temp = [];
            Im_temp = [];
            for j = 1:m
                for k=1:z
                    Re_temp = [Re_temp;real(XF{p}{j,k}(i))];
                    Im_temp = [Im_temp;imag(XF{p}{j,k}(i))];
                end
            end
            [rXSmt,rSSmt,Param,rXPos,rSPos] = ay_fft_filter([Re_temp Im_temp],Iter);
            % Real_x = rXPos;
            xy     = ay_multi_state_sample(Ns,rXSmt,rSSmt,rXPos,rSPos,Param);
            Real_x = squeeze(xy(:,:,1));   
            Img_x  = squeeze(xy(:,:,2));

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
    % Global coherence next
    % we have multiple measure of FFT per Wnd - number: Wnd/Ls
    for i = 1:m
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
                    if method==1 || method==2
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
                            % Ta(t) = complex(abs(Fa(t))*cos(angle(Fa(t))-pTa),abs(Fa(t))*sin(angle(Fa(t))-pTa));
                            % Tb(t) = complex(abs(Fb(t))*cos(angle(Fb(t))-pTb),-abs(Fb(t))*sin(angle(Fb(t))-pTb));
                            Ta(t) = complex(cos(angle(Fa(t))-pTa),sin(angle(Fa(t))-pTa));
                            Tb(t) = complex(cos(angle(Fb(t))-pTb),-sin(angle(Fb(t))-pTb));
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
end


COH_F   = COH;
Extra.S = zeros(1,size(COH_F,2));
Extra.Xs = zeros(size(COH)); % smoother
Extra.Ws = zeros(size(COH));
Extra.Xf = zeros(size(COH)); % filter
Extra.Wf = zeros(size(COH));

for i = 1:size(COH_F,2)
    Yn = COH(:,i);
    n  = length(Yn);
    In = ones(n,1);
    valid = ones(n,1);

    %% Set Behavioral Model and Learning Procedure
    %  Create model
    Param = ay_create_state_space(1,0,1,0,1,1,0,1,1);
    Param.Ak = 1;
    Param.Ck = 1;
    % Set learning parameters
    Param = ay_set_learning_param(Param,Iter,0,1,1,1,1,0,1,2,1);
    % EM
    [rXSmt,rSSmt,Param,rXPos,rSPos,~,EYn]=ay_em([2 0],[],In,0,Yn,[],Param,valid);

    % Extract estimating data
    COH_F(:,i) = EYn + Param.S;
    Extra.S(i) = Param.S;
    Extra.Xs(:,i) = cell2mat(rXSmt);
    Extra.Ws(:,i) = cell2mat(rSSmt);
    Extra.Xf(:,i) = cell2mat(rXPos);
    Extra.Wf(:,i) = cell2mat(rSPos);
end




