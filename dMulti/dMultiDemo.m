%==========================================================================
%              demo for decompositon of multi-channel sEMG                *
%                                                                         *
% Read sEMG Data                                                          *
% Format: Multi - Channel sEMG Data should be organized into a M * N      *
% matrix, which contain N channels and M points in each channel.          *
% x = [x_ch1_time1, ... , x_chN_time1;                                    *
%      x_ch1_time2, ... , x_chN_time2;                                    *
%              ..., ... ,         ...;                                    *
%      x_ch1_timeM, ... , x_chN_timeM]                                    *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================
clear all;clc;
close all;
%%
load('R05M1T25L25_trial1');
% load('data_tringle.mat');
%%
data = cell(40,10);
for i = 1:1:40
    for j = 1:2:19
        data{i,(j-1)/2+1}=EMG{i,j};
    end
end
for j=1:10
    for i=1:39
        data{i,j}=data{i,j}-data{i+1,j};
    end
end
for i = 1:1:39
    for j = 1:10
        sEMG.data((i-1)*10+j,:)=data{i,j};
    end
end

sEMG.dt = 1/2048;
sEMG.t0 = 0;
sEMG.ch = [39,10];
[~,N] = size(sEMG.data);
t = (sEMG.t0:N-1)' * sEMG.dt;

figure
for j=20
    for i=1:sEMG.ch(2)
        select_ch = [j,i];
        plot(t,sEMG.data((j-1)*sEMG.ch(2)+i,:)+300*(i-1))
        hold on
    end
end

figure
for i=1:sEMG.ch(1)
    select_ch = [i,9];
    plot(t,sEMG.data((i-1)*sEMG.ch(2)+select_ch(2),:)+300*i+100*40);
    hold on
end
% for i=1:39
%     select_ch = [i,10];
%     plot(t,sEMG.data{select_ch(1),select_ch(2)}+100*(i+41))
%     hold on
% end



































