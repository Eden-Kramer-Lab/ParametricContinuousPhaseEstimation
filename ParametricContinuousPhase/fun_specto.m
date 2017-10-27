function [] = fun_specto(data,discription)
[m,n] = size(data);
number_of_subplot = 4;
window=200;
noverlap=150;
f=1000;
fs=1000;

for i=1:m/number_of_subplot
    figure('units','normalized','outerposition',[0 0 1 1]);
    for j=1:number_of_subplot
        temp_data = data(4*(i-1)+j,:);
        subplot(2,2,j),spectrogram(temp_data,window,noverlap,f,fs),title(['spetcrogram of chanal #',num2str(4*(i-1)+j)]);
    end
    suptitle(discription)
end

end