function  p_in_range = ay_phase(rXSmt,rSSmt,ph_range)
sample = 1e4;
K = length(rXSmt);
p_in_range = zeros(K,1);
for k=1:K
    Mx = rXSmt{k};
    Sx = rSSmt{k};Sx = 0.5 * (Sx +Sx');
    xy = mvnrnd(Mx,Sx,sample);
    dp  = xy(:,1)./sqrt(xy(:,1).^2+xy(:,2).^2);
    ind = find(dp>=ph_range(1) & dp<=ph_range(2));
    p_in_range(k)=length(ind)/sample;
end