%% Load the data of filtering result
load('filter_Iter_300.mat');
load('filter_Iter_30.mat')
load('Subset_COH.mat');
% Extract data
for n = 1: 500
    for m = 1:200
    filter_X(m,n) = B_rXPos{1,n}{m};
    filter_Y(m,n) = B_EYn{1,n}(m)+B_Param{1,n}.S;
    end
end

for n = 1: 500
    for m = 1:200
    filter1_X(m,n) = B_rXPos{1,n}{m};
    filter1_Y(m,n) = B_EYn{1,n}(m)+B_Param{1,n}.S;
    end
end

% Visualize
figure()
subplot(1,3,1)
imagesc(COH)
colorbar
title('COH without filtering')

subplot(1,3,2)
imagesc(filter1_Y)
colorbar
title('Iter = 30, X+C0')

subplot(1,3,3)
imagesc(filter_Y)
colorbar
title('Iter = 300, X')