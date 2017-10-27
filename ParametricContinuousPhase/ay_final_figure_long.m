%% gerate data
ch_no   = 8;
sample_rate = 1000;
EEG_8 = generate_scenario2_signal(ch_no,sample_rate,160,5,3);


% EEG_1 = generate_scenario2_signal(4,1000,40,40,3);

%% no filter(fft) + no filter(global coherence)
fs_a= 20;
fs_b= 90;
Iter= 200;
wnd_len = 1000;
multiple_wnd = 8;
COH_8_1=ay_global_coherence(EEG_8,wnd_len,wnd_len*multiple_wnd,fs_a,fs_b);

figure()
imagesc(COH_8_1);

%% filter(fft) + no filter(global coherence)
wnd_len=1000;
COH2_test = ypl_global_coherence(EEG_8,wnd_len,wnd_len*multiple_wnd,size(EEG_8,1),Iter,fs_a,fs_b);

figure()
imagesc(COH2_test);

%% no filter(fft) + filter(global coherence)
[B_Param,B_EYn,B_rXPos,B_rSPos] = yl_COH_filter(COH2_test,Iter);
for n = 1: length(B_EYn)
    for m = 1: length(B_EYn{1})
        test_Y8(m,n) = B_EYn{1,n}(m)+B_Param{1,n}.S;
    end
end

for n = 1: length(B_EYn)
    for m = 1: length(B_EYn{1})
        test_Y8(m,n) = B_EYn{1,n}(m)+B_Param{1,n}.S;
    end
end
   
figure()
imagesc(test_Y8);


B_test_Y=[];
B__test_X = [];
%% filter(fft) +filter(global coherence)

[B1_Param,B1_EYn,B1_rXPos,B1_rSPos] = yl_COH_filter(COH_8_1,Iter);

for n = 1: length(B1_EYn)
    for m = 1: length(B1_EYn{1})
    test_Y8_1(m,n) = B1_EYn{1,n}(m)+B1_Param{1,n}.S;
    end
end


figure()
subplot(2,2,1)
imagesc(COH_8_1);
title('PNAS Paper Method');
xlabel('Frequency');
ylabel('Time Index');
xticklabels(39:20:79)
caxis([0,1])
colorbar

subplot(2,2,2)
imagesc(COH2_test);
title('Filtered FFT Only');
xlabel('Frequency');
ylabel('Time Index');
xticklabels(39:20:79)
caxis([0,1])
colorbar

subplot(2,2,3)
imagesc(test_Y8_1)
title('PNAS Method Followed By Gamma Filter');
xlabel('Frequency');
ylabel('Time Index');
xticklabels(39:20:79)
caxis([0,1])
colorbar

subplot(2,2,4)
imagesc(test_Y8);
title('Filtered FFT Followed By Gamma Filter')
xlabel('Frequency');
ylabel('Time Index');
xticklabels(39:20:79)
caxis([0,1])
colorbar



























;







