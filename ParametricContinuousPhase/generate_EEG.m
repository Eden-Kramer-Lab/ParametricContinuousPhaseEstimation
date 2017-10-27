COH2 = ypl_global_coherence(EEG_1_1,1000,8000,4,50);
COH3 = ypl_global_coherence(EEG_1_2,1000,8000,4,50);
COH4 = ypl_global_coherence(EEG_1_3,1000,8000,4,50);

COH2_test = ypl_global_coherence(EEG_1_1,1000,8000,4,200);
COH3_test = ypl_global_coherence(EEG_1_2,1000,8000,4,200);
COH4_test = ypl_global_coherence(EEG_1_3,1000,8000,4,200);




COH_2_1=ay_global_coherence(EEG_1_1,1000,8000);
COH_3_1=ay_global_coherence(EEG_1_2,1000,8000);
COH_4_1=ay_global_coherence(EEG_1_3,1000,8000);


figure()
imagesc(COH_2_1);

specgram(EEG_1_1(1,:));