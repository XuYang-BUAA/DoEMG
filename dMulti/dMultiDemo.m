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

select_ch = [20,10];
figure
plot(t,sEMG.data{select_ch(1),select_ch(2)})
hold on



































