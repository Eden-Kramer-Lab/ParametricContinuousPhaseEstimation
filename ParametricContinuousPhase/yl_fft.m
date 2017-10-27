function [XF,filter_Y] = yl_fft(X,Ls,Iter)

Lt = length(X)/Ls;
XF = ones(Lt,Ls/2);
for i = 1:Lt
    ind_a = (i-1)*Ls+1;
    ind_b = i*Ls;
    temp1 = fft(X(ind_a:ind_b));
    XF(i,:) = temp1(1:Ls/2)/Ls;
end
% XF is a matrix of Lt* Fs

for i = 1:Ls/2
    Y = XF(:,i);
    Re_Y = real(Y);
    img_Y = imag(Y);
     
    disp(['Iteration' num2str(i)])
    
    [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb] = yl_filter(Re_Y,Iter);
     Real_XSmt{i} = rXSmt;
     Real_SSmt{i} = rSSmt;
     Real_rXPos{i} = rXPos;
     Real_rSPos{i} = rSPos;

     Real_Param{i} = Param;
     Real_EYn{i} = EYn;
     Real_rYn{i} = rYn;
    [rXSmt,rSSmt,Param,rXPos,rSPos,ML,EYn,EYb,rYn,rYb] = yl_filter(img_Y,Iter);
     img_XSmt{i} = rXSmt;
     img_SSmt{i} = rSSmt;
     img_rXPos{i} = rXPos;
     img_rSPos{i} = rSPos;

     img_Param{i} = Param;
     img_EYn{i} = EYn;
     img_rYn{i} = rYn;    
     
end
% Extra data
filter_Y =[];
for i = 1:500
    filter_Y(i,:) = complex(cell2mat(Real_rXPos{i}),cell2mat(img_rXPos{i}));
end
end