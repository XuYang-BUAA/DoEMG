clear all;clc;
% =========================================================================
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
%    06/29/2020 XY : Update.                                              *
% =========================================================================

load('decomp_lpy_cycle.mat');
sEMG = signals{1,1};

sEMG.data = signals{1,1}.data(1:81920,:);
sEMG.dt = signals{1,1}.dt;
sEMG.t0 = signals{1,1}.t0;
sEMG.chn_num = signals{1, 1}.chn_num;

sEMG.data = sEMG.data - repmat(mean(sEMG.data), length(sEMG.data), 1);
N = size(sEMG.data);
% Generate Time Label
t = (sEMG.t0:N-1)' * sEMG.dt;

%[result, resdata] = PD3(sEMG, 0.015);
result = PD3(sEMG, 0.015,0);
%%
plot_decomp_result(sEMG,result);
%%
fprintf('Finish!\n');
