function [ind_timewin_min, ind_timewin_max, timewin, pos_cluster_electrode, neg_cluster_electrode, p_pos, p_neg] = find_cluster(stat, alpha)
%shamelessly copied from the
% FT_CLUSTERPLOT fieldtrip function

% first squeeze all values from the statmatrix
stat.prob = squeeze(stat.prob);
stat.posclusterslabelmat = squeeze(stat.posclusterslabelmat);
stat.negclusterslabelmat = squeeze(stat.negclusterslabelmat);
stat.cirange = squeeze(stat.cirange);
stat.mask = squeeze(stat.mask);
stat.stat = squeeze(stat.stat);
stat.ref = squeeze(stat.ref);



cfg.alpha = alpha;

% this determines the labels in the figure
hastime = isfield(stat, 'time');
hasfreq = isfield(stat, 'freq');

% use the vector time, even though the 2nd dimension might be freq
if hastime
    time = stat.time;
elseif hasfreq
    time = stat.freq;
end

if issubfield(stat, 'cfg.correcttail') && ((strcmp(stat.cfg.correcttail, 'alpha') || strcmp(stat.cfg.correcttail, 'prob')) && (stat.cfg.tail == 0));
    if ~(cfg.alpha >= stat.cfg.alpha)
        ft_warning(['the pvalue you plot: cfg.alpha = ' num2str(cfg.alpha) ' is higher than the correcttail option you tested: stat.cfg.alpha = ' num2str(stat.cfg.alpha)]);
    end
end

% find significant clusters
haspos = isfield(stat, 'posclusters');
hasneg = isfield(stat, 'negclusters');

if haspos == 0 && hasneg == 0
    fprintf('%s\n', 'no clusters exceeded the nominal threshold in the data; nothing to plot')
else
    if haspos
        probpos = [stat.posclusters.prob];
        sigpos  = probpos < cfg.alpha;
    else
        sigpos  = [];
    end
    if hasneg
        probneg = [stat.negclusters.prob];
        signeg  = probneg < cfg.alpha;
    else
        signeg  = [];
    end
    sigpos  = find(sigpos == 1);
    signeg  = find(signeg == 1);
    Nsigpos = length(sigpos);
    Nsigneg = length(signeg);
    Nsigall = Nsigpos + Nsigneg;
    
    if Nsigall == 0
        ft_error('no clusters present with a p-value lower than the specified alpha, nothing to plot')
    end
    
    % make clusterslabel matrix per significant cluster
    if haspos
        posCLM = stat.posclusterslabelmat;
        sigposCLM = zeros(size(posCLM));
        for iPos = 1:Nsigpos
            sigposCLM(:,:,iPos) = (posCLM == sigpos(iPos));
            %hlsignpos(iPos) = prob2hlsign(probpos(iPos), cfg.highlightsymbolseries);
        end
    else
        sigposCLM = [];
        probpos = [];
    end
    
    if hasneg
        negCLM = stat.negclusterslabelmat;
        signegCLM = zeros(size(negCLM));
        for iNeg = 1:Nsigneg
            signegCLM(:,:,iNeg) = (negCLM == signeg(iNeg));
            %hlsignneg(iNeg) = prob2hlsign(probneg(iNeg), cfg.highlightsymbolseries);
        end
    else % no negative clusters
        signegCLM = [];
        probneg = [];
    end
    
    fprintf('There are %d clusters smaller than alpha (%g)\n', Nsigall, cfg.alpha);
    
    
    % define time or freq window per cluster
    for iPos = 1:Nsigpos
        possum_perclus = sum(sigposCLM(:,:,iPos),1); %sum over chans for each time- or freq-point
        ind_min = find(possum_perclus~=0, 1 );
        ind_max = find(possum_perclus~=0, 1, 'last' );
        time_perclus = [time(ind_min) time(ind_max)];
    end
    for iNeg = 1:Nsigneg
        negsum_perclus = sum(signegCLM(:,:,iNeg),1);
        ind_min = find(negsum_perclus~=0, 1 );
        ind_max = find(negsum_perclus~=0, 1, 'last' );
        time_perclus = [time(ind_min) time(ind_max)];
    end
    
    % define time- or freq-window containing all significant clusters
    possum = sum(sum(sigposCLM,3),1); %sum over Chans for timevector
    negsum = sum(sum(signegCLM,3),1);
    
    if haspos && hasneg
        allsum = possum + negsum;
    elseif haspos
        allsum = possum;
    else
        allsum = negsum;
    end
    
    % first and last time points of any cluster
    ind_timewin_min = find(allsum~=0, 1);
    ind_timewin_max = find(allsum~=0, 1, 'last');
    timewin = time(ind_timewin_min:ind_timewin_max);
    
    % find all electrodes that are part of a cluster
    pos_cluster_ind = find(sum((sigposCLM~=0),2) ~= 0);
    pos_cluster_electrode = stat.label(pos_cluster_ind);
    neg_cluster_ind = find(sum((signegCLM~=0),2) ~= 0);
    neg_cluster_electrode = stat.label(neg_cluster_ind);
    
    % find significance value of positive or negative cluster
    if ~isempty(sigpos) 
        p_pos = stat.posclusters(sigpos).prob;
    else
        p_pos = [];
    end
    if ~isempty(signeg)
        p_neg = stat.negclusters(signeg).prob;
    else
        p_neg = [];
    end
        
    
end