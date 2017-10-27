function [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb]= ay_fft_filter(Yn,Iter)

%% Parameters
Ak  = eye(2,2);
Ck  = eye(2,2);
Wk  = [1 eps;eps 1];
Vk  = [1 eps;eps 1];
X0  = [0;0];
W0  = Wk;

 
K = size(Yn,1);
for iter=1:Iter
    %% data length
    K = size(Yn,1);
    %% variables
    % one-step mean and variance
    XPre = cell(K,1);
    SPre = cell(K,1);
    % filter mean and covariance
    XPos = cell(K,1);
    SPos = cell(K,1);
    % As
    As   = cell(K,1);
    % smoother mean and covariance
    XSmt = cell(K,1);
    SSmt = cell(K,1);
    % extra components
    Ckk   = cell(K,1);
    Ckk_1    = cell(K,1);
    Wkk_1    = cell(K,1);
    

    %% run filter step
    for k=1:K
        %% Run one step prediction
        if k == 1
            XPre{k} = Ak * X0 ;
            SPre{k} = Ak * W0 * Ak'+ Wk;
        else
             XPre{k} = Ak * XPos{k-1};
             SPre{k} = Ak * SPos{k-1}* Ak' + Wk;
        end
        %% Run filter
        Yp =  Yn(k,:)' - Ck * XPre{k};
        Sk =  Ck * SPre{k} * Ck' + Vk;
        Kk =  SPre{k}* Ck' * pinv(Sk);
        XPos{k} =  XPre{k} + Kk * Yp;
        SPos{k} = (I-Kk*Ck)*SPre{k};
    end
    
    %% run smoother step
    XSmt{end} = XPos{end};
    SSmt{end} = SPos{end};
    for k=K-1:-1:1
        % Ak, equation (A.10)
        As{k} = SPos{k} * Ak' *  SPre{k+1}^-1 ;
        % Smting function, equation (A.9)
        XSmt{k} = XPos{k} + As{k} * (XSmt{k+1}- XPre{k+1});
        % Variance update, equation (A.11)
        SSmt{k} = SPos{k} + As{k} * (SSmt{k+1}- SPre{k+1}) * As{k}';
    end
    % Kalman smoother for time 0
    As0   = W0 * Ak' *  SPre{1}^-1;
    XSmt0 = X0 + As0 * (XSmt{1}- XPre{1}) ;
    SSmt0 = W0 + As0 * (SSmt{1}- SPre{1}) * As0';
    
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk)
    % Ckk = E(Xk*Xk) prediction by smoothing 
    for k = 1:K
        % Wk update - Smting Xk*Xk
        Ckk{k} = SSmt{k} + XSmt{k} * XSmt{k}';
    end
    Ckk0 = SSmt0 + XSmt0 * XSmt0';
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk-1)
    % Ckk_1=E(Xk-1*Xk) prediction by smoothing - it is kept at index K
    % Wkk_1= Ckk_1 + Bias
    % Covariance for smoothed estimates in state space models - Biometrica
    % 1988- 601-602
    for k = 1:K
        % Wkk update - Smoothing Xk-1*Xk
        if k>1
            Wkk_1{k}   = As{k-1} * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt{k-1}'; 
        else
            Wkk_1{k}   = As0 * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt0'; 
        end
    end
    
    %% Parameter Update
    
    
    %% Parameters
    Param.Ak = Ak;
    Param.Ck = Ck;
    Param.Wk = Wk;
    Param.Vk = Vk;
    Param.W0 = W0;
    Param.X0 = X0;
end        

