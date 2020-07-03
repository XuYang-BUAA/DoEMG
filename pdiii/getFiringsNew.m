function [ firings, resdata ] = getFirings(rdata, mu_template,thresh,spike_width)
% =========================================================================
%                          Get MU firing trains.                          *
%                                                                         *
%  INPUT:                                                                 *
%    rdata           -- sEMG data                                         *
%    spike_segment   -- spike segments                                    *
%    spike_train     -- timings of firing segments                        *
%    mu_template     -- mu templates                                      *
%                    -- members:                                          *
%                    -- shape      :template signal segments              *
%                    -- ch_num     :channel number                        *
%                    -- tmplt_mem  :memory of all templates               *
%                    -- mu_num     :MU number                             *
%                    -- mu_len     :length of MU segment                  *
%                    -- fire_times :firing times of each MU               *
%                    -- mu_id      :id index of MU templates              *
%                                                                         *
%  OUTPUT:                                                                *
%    firings         -- merged MU templates                               *
%    resdata         -- residual data                                     *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/2/2020 : XuY created                                               *
% =========================================================================
%%
% =========== Basic information & Settings=================================
    data_len = length(rdata);
    resdata = rdata;
    
    mu_num = mu_template.mu_num;
    mu_len = mu_template.mu_len;
    %mu_lenH = floor(mu_len/2);
    ch_num = mu_template.ch_num;
    mu_shapes = mu_template.shape;
    mu_ids = mu_template.mu_id;
    tmplt_mem = mu_template.tmplt_mem;
    
    
%%
%==========================================================================
    [mu_firings , res_data] = checkMaxMu(resdata,mu_template,thresh,spike_width);



























