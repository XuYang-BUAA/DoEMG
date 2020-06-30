clc;clear;
load('testdata.mat');
data = filter(filter5khz_bandpass,sEMG.data);
%data = sEMG.data;
[data_len,ch_num]=size(data);
x = 1:1:data_len;
mu_width = 75;
bias =10;
figure()
threshold = 0.1;
plot(x,data(:,:));
hold on
%%
[peaks, timings] =  findpeaks(abs(data(:,1)), ...
            'MinPeakHeight', threshold);
plot(timings,data(timings,1),'o')
