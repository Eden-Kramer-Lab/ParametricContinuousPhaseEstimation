s1    = 1;rho12=0.4;
s2    = 1;rho13=0.9;
s3    = 1;rho23=-0.8;
df    = 16;
% draw sample
for z=-9:9
    rho12 = z/10;
    res = [];
    %sigma = [s1 sqrt(s2*s1)*rho;sqrt(s2*s1)*rho s2];
    sigma = rand(4,4);
    sigma = 0.5*(sigma+sigma');
    sigma = sigma + 4*eye(4);
    for i=1:1e4
        out = wishrnd(sigma,df);
        s   = svd(out/df);
        res = [res;s(1)/(s(1)+s(2)+s(3)+s(4))];
    end
    subplot(2,1,1)
    hist(res,[0:0.01:1])
    xlim([0 1])
    title(['rho=' num2str(rho)])
    subplot(2,1,2)
    pf = hist(res,[0:0.01:1]);
    phat = betafit(res);
    of = betapdf([0:0.01:1],phat(1),phat(2));
    plot([0:0.01:1],of);
    hold on
    plot([0:0.01:1],pf/1e2);
    hold off
    pause()
end