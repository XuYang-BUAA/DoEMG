function spikes = detectSpikes(data,ch_num)
%==========================================================================
%                        detect spikes                                    *
%                                                                         *
% INPUT:                                                                  *
%    sEMG            -- filtered multi-channel sEMG data                  *
%                                                                         *
% OUTPUT:                                                                 *
%    spikes          -- spikes for clustering                             *
%                                                                         *
%                                                                         *
%                                                                         *
%  WARNINGS:                                                              *
%    None                                                                 *
%                                                                         *
%  HISTORY:                                                               *
%    07/08/2020  : XuY create.                                            *
%==========================================================================