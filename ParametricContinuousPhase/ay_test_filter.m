K  = 1000;
X  = zeros(K,2);
Y  = zeros(K,2);
Ws  = [0.2 0.03;0.03 0.1];
Sxx  = 0.1;
Syy  = 0.1; 
X0 = mvnrnd([0;0],Ws);
for k=1:K
    if k==1
        Mx =  X0;
        Sx =  Ws;
        X(k,:)=mvnrnd(Mx,Sx);
    else
        Mx =  X(k-1,:)';
        Sx =  Ws;
        X(k,:)=mvnrnd(Mx,Sx);
    end
    Y(k,1)= X(k,1) + sqrt(Sxx)*randn();
    Y(k,2)= X(k,2) + sqrt(Syy)*randn();
end
Iter= 100;
[rXSmt,rSSmt,Param,rXPos,rSPos] = ay_fft_filter([Y(:,1) Y(:,2)],Iter);

subplot(2,2,2)
K  = length(rXSmt);
Xm = zeros(K,1);
Xb = zeros(K,1);
for k=1:K
    Xm(k) = rXSmt{k}(1);
    Xb(k) = rSSmt{k}(1,1);
end
ay_plot_bound(1,(1:K),Xm',(Xm-2*sqrt(Xb))',(Xm+2*sqrt(Xb))');

subplot(2,2,4)
Xm = zeros(K,1);
Xb = zeros(K,1);
for k=1:K
    Xm(k) = rXSmt{k}(2);
    Xb(k) = rSSmt{k}(2,2);
end
ay_plot_bound(1,(1:K),Xm',(Xm-2*sqrt(Xb))',(Xm+2*sqrt(Xb))');
