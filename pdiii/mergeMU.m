function [ newMUPool ] = mergeMU( mus, prob )
% =========================================================================
%Merge the MUs which have a high probability being split from a single MU *
%      comparing template sj with sk                                      *
%      P = B(sj'*sk)/(|sj||sk|);                                          *
%      P is the estimated prbalility that sj,sk are from the same train   *
%      A = sj' * sk/(sj' * sj);                                           *
%      B = A for 0<=A<=1, B = 1/A for A>1, B=0 for A<0.                   *
%                                                                         *
%  INPUT:                                                                 *
%    mus             -- mu templates                                      *
%                    -- members:                                          *
%                    -- shape      :template signal segments              *
%                    -- ch_num     :channel number                        *
%                    -- tmplt_mem  :memory of all templates               *
%                    -- mu_num     :MU number                             *
%                    -- mu_len     :length of MU segment                  *
%                    -- fire_times :firing times of each MU               *
%                    -- mu_id      :id index of MU templates              *
%                                                                         *
%    prob            -- probabilty threshold for merging                  *
%                                                                         *
%  OUTPUT:                                                                *
%    newMUPool       -- merged MU templates                               *
%                                                                         *
%  WARNINGS:   none                                                       *
%                                                                         *
%  HISTORY:                                                               *
%    7/2/2020 : XuY update                                                *
% =========================================================================
    mu_num = mus.mu_num;
    templates = mus.shape;
    
    ss = sum(templates.^2);
    xss = sqrt(ss' * ss);
    [~,~,~,A] = caldiff(templates, templates);

    B = A'; B(B>1) = 1./B(B>1); B(B<0) = 0;
    % Calculate inclusion measure.
    Pin = B .* (templates' * templates) ./ xss;
%     Pin = calASR(templates,templates);
    P = Pin;
    P(1:mu_num+1:end) = 0;
    [r,c] = find(P>prob);
    G = graph(r,c,P(r+(c-1)*mu_num),cellstr(num2str((mus.mu_id)')));
    % figure; plot(G,'NodeLabel',G.Nodes.Name,'EdgeLabel',G.Edges.Weight);
%     title(sprintf('P = %f', prob));
    [bin, binsize] = conncomp(G);
    bin_idxs = find(binsize > 1);
    mu_del_mask = false(1,mu_num);
    for bin_idx = bin_idxs
        mu_inx = find(bin == bin_idx);
        merge_id = mus.mu_id(mu_inx);
        tmp_flag = sum(eq(mus.tmplt_mem(1,:), merge_id')) > 0;
        mus.tmplt_mem(1,tmp_flag) = merge_id(1);
        mus.mu_id(mu_inx) = merge_id(1);
        merge_num = length(mu_inx);
        fire_times = mus.fire_times(mu_inx);
        new_tmplt = sum(templates(:,mu_inx).*fire_times./sum(fire_times), 2);
        mus.shape(:,mu_inx) = repmat(new_tmplt, 1, merge_num);
        mus.fire_times(:,mu_inx) = sum(fire_times);
        mus.mu_num = mus.mu_num - sum(mu_inx) + 1;
        mu_del_mask(mu_inx(2:end)) = true;
    end
    mus.shape(:,mu_del_mask) = [];
    mus.fire_times(:,mu_del_mask) = [];
    mus.mu_id(mu_del_mask) = [];
    mus.mu_num = mu_num - sum(mu_del_mask);
    newMUPool = mus;
end

%--------------------------------------------------------------------------
function [cost] = cal_cost(pmat, group)
% Calculate merge cost
    [all_num,~] = size(pmat);
    in_num = length(group); apair = in_num*(in_num);
    ex_num = all_num - in_num; xpair = 2*(in_num*ex_num);
    Pex = sum(sum(1 - pmat(group,group)));
    Pin = sum(sum([pmat(group,:);pmat(:,group)']));
    cost = Pex/apair + (Pin+2*sum(sum(pmat(group,group))))/xpair;
end

%------------------------------EOF-----------------------------------------
