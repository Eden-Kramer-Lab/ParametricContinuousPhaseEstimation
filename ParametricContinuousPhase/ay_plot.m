F  = COHBB;
G  = COHBB_F;

subplot(1,3,1)
imagesc(30:50,8:8:720,F);
xlabel('Frequency')
ylabel('Time (sec)')

subplot(1,3,2)
imagesc(30:50,8:8:720,G);
xlabel('Frequency')
ylabel('Time (sec)')

subplot(2,3,3)
plot(F(:,12),'LineWidth',2)
hold on
plot(G(:,12),'LineWidth',2)
hold off
ylabel('G_c')
axis tight
ylim([0 1])


subplot(2,3,6)
plot(F(:,15),'LineWidth',2)
hold on
plot(G(:,15),'LineWidth',2)
hold off
ylabel('G_c')
xlabel('Time Update')
axis tight
ylim([0 1])