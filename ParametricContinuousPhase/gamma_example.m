rho=[0.8 0.75 0.7 0.65 0.60 0.55 0.59 0.68 0.76 0.85]
%%
nd = 32;
U  = orth(randn(nd,nd));

for i = 1:10
    L(:,:,i) = diag(exp(-rho(i)*(0:nd-1)));
    S(:,:,i) = U*L(:,:,i)*U';
    M(:,i) = svd(U*L(:,:,i)*U');
    g(i) = M(1,i)/sum(M(:,i))
   % g(i) = L{i,1}/sum(L{i,:})

end

%% generate 10 sessions of data
df = 64;
Gs = [];

for i=1:30
    T = wishrnd(S(:,:,1),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,2),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,3),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,4),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,5),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,6),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,7),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,8),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,9),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

for i=1:30
    T = wishrnd(S(:,:,10),df);
    ls = svd(T);
    Gs = [Gs;ls(1)/sum(ls)];
end

true_g = [g(1)*ones(30,1);g(2)*ones(30,1);g(3)*ones(30,1);g(4)*ones(30,1);
          g(5)*ones(30,1);g(6)*ones(30,1);g(7)*ones(30,1);g(8)*ones(30,1);
          g(9)*ones(30,1);g(10)*ones(30,1)];
plot(Gs);hold on;plot(true_g);hold off

%%
In = ones(300,2);

%create model
Param = ay_create_state_space(1,0,2,0,1,1,0,1,1);
Param.Ak = 1;

%set learning parameters
Iter = 300;
Param = ay_set_learning_param(Param,Iter,0,1,1,1,1,0,1,2,1)

%
[rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb]=ay_em([2 0],[],In,0,Gs,[],Param,ones(300,1))

%visualization
figure(2)
K  = length(Gs);
xm = zeros(K,1);
xb = zeros(K,1);
for i=1:K
    temp=rXSmt{i};xm(i)=temp(1);
    temp=rSSmt{i};xb(i)=temp(1,1);
end
ay_plot_bound(1,(1:K),xm,(xm-2*sqrt(xb))',(xm+2*sqrt(xb))');
ylabel('x_k');
xlabel('Trial');
hold on 
plot(true_g)
hold on
plot(EYn+Param.S)
hold off

figure(3)
ml=[];
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');
