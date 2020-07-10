function [decomp_result] = dMulti(sEMG)
%==========================================================================
%              top of decompositon of multi-channel sEMG                  *
%                                                                         *
% INPUT:                                                                  *
%    sEMG            -- filtered multi-channel sEMG data                  *
%                                                                         *
% OUTPUT:                                                                 *
%    decomp_result   -- cell for the decompsition result                  *
%                                                                         *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================
% ========== Data information ========================================-----
    data = sEMG.data;
    % fdata = doFilter10KHz(data);
    fdata = data;
    t0 = sEMG.t0; dt = sEMG.dt;
    ch_num = sEMG.ch_num;
    [~,data_len] = size(sEMG.data{1,1});
    % Time series
    t = t0 + dt * (0:data_len-1);
    
    detectSpikes(sEMG,ch_num);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    