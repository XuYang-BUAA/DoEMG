%==========================================================================
%                         demo for PD3                                    *
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
%    06/29/2020  : XuY Update.                                            *
%==========================================================================
%%
clear all;clc;
close all;
%%

% load('triangle.mat');
% sEMG.data = data(1,:)';
% sEMG.dt = 1/10240;
% sEMG.t0 = 0;
% sEMG.chn_num = 1;

load('decomp_lpy_cycle.mat');
sEMG = signals{1,1};
sEMG.data = signals{1,1}.data(600000:680000,:);
newdata = filter(filter5khz,sEMG.data);
sEMG.data = newdata;
sEMG.dt = signals{1,1}.dt;
sEMG.t0 = signals{1,1}.t0;
sEMG.chn_num = signals{1, 1}.chn_num;

sEMG.data = sEMG.data - repmat(mean(sEMG.data), length(sEMG.data), 1);
N = size(sEMG.data);
t = (sEMG.t0:N-1)' * sEMG.dt;



% Generate Time Label


figure()
plot(t,sEMG.data+2);
%plot(sEMG.data+2);
hold on;

plot(t,newdata);
hold on

%[result, resdata] = PD3(sEMG, 0.015);
result = PD3(sEMG, 0.02,0);

%%
plot_decomp_result(sEMG,result);
%%
fprintf('Finish!\n');
