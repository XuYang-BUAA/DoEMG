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
%%
sEMG.data = EMG;
sEMG.dt = 1/2048;
sEMG.t0 = 0;
sEMG.chn_num = [40,20];

[~,N] = size(sEMG.data{1,1});
t = (sEMG.t0:N-1)' * sEMG.dt;

for j=1:20
    for i=1:39
        sEMG.data{i,j}=sEMG.data{i,j}-sEMG.data{i+1,j};
    end
end
figure
for j=1
    for i=1:20
        select_ch = [j,i];
        plot(t,sEMG.data{select_ch(1),select_ch(2)}+100*i+300*j)
        hold on
    end
end
% for i=1:39
%     sEMG.data{i,9}=sEMG.data{i,9}-sEMG.data{i+1,9};
% end

% sEMG.data{20,8}=sEMG.data{20,8}-sEMG.data{20,9};
% sEMG.data{20,10}=sEMG.data{20,10}-sEMG.data{20,9};
% sEMG.data{20,9}=sEMG.data{20,9}-sEMG.data{19,9};


t = (sEMG.t0:N-1)' * sEMG.dt;

figure
for i=1:39
    select_ch = [i,9];
    plot(t,sEMG.data{select_ch(1),select_ch(2)}+100*i+100*40)
    hold on
end
% for i=1:39
%     select_ch = [i,10];
%     plot(t,sEMG.data{select_ch(1),select_ch(2)}+100*(i+41))
%     hold on
% end



































