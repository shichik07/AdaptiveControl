function NEW_TBL =  theta_mat_to_R(freq_mat, cfg)

% group indices
grp_CO= cfg.grp_idx == 0;
grp_PD= cfg.grp_idx == 1;


% load electrodes of significant cluster
 save_folder = 'E:\AdaptiveControl\Data\TF_for_R';
% get electrode load location
if cfg.Ana_type == "LWPC"
    if cfg.lck =="SL"
        save_name_CO = ['Cluster_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_PD_' + cfg.Ana_type + '_Theta_electrodes.csv'];
    elseif cfg.lck == "RL"
        save_name_CO = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_RL_PD_' + cfg.Ana_type + '_Theta_electrodes.csv'];
    end
% for the ISPC effect we only load the CO cluster, because there was no significant cluster in the PD group 
elseif cfg.Ana_type == "ISPC"
    if cfg.lck =="SL"
        save_name_CO = ['Cluster_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
    elseif cfg.lck == "RL"
        save_name_CO = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
    end
end

% load the data
save_loc_PD = fullfile(save_folder,save_name_PD);
cluster_PD = readtable(save_loc_PD);
elec_idx_PD = find(ismember(cfg.chanlocs, table2cell(cluster_PD)'));
save_loc_CO = fullfile(save_folder,save_name_CO);
cluster_CO = readtable(save_loc_CO);
elec_idx_CO = find(ismember(cfg.chanlocs, table2cell(cluster_CO)'));


% get correct electrode indices

% get table
NEW_TBL1 = table;

% get the PD data
freq_mat_PD              = double(squeeze(mean(freq_mat(elec_idx_PD, cfg.freq_idx, :, grp_PD),1))); % average over electrodes
freq_mat_PD              = double(squeeze(mean(freq_mat_PD,1))); % average over frequencies
NEW_TBL1.SE              = transpose(std(freq_mat_PD')/ sqrt(size(freq_mat_PD,2))); % Get SE
NEW_TBL1.Pw              = squeeze(mean(freq_mat_PD,2)); %  averageover participant
NEW_TBL1.Condition       = repmat({cfg.var_n}, length(freq_mat_PD), 1);
NEW_TBL1.Group           = repmat({"PD"}, length(freq_mat_PD),1);
NEW_TBL1.Lock            = repmat({cfg.lck}, length(freq_mat_PD),1);
NEW_TBL1.Ana_type        = repmat({cfg.Ana_type}, length(freq_mat_PD),1);
NEW_TBL1.Time            = transpose(cfg.Time);

NEW_TBL2 = table;
% get the CO data
freq_mat_CO              = double(squeeze(mean(freq_mat(elec_idx_CO, cfg.freq_idx, :, grp_CO),1))); % average over electrodes
freq_mat_CO              = double(squeeze(mean(freq_mat_CO,1))); % average over frequencies
NEW_TBL2.SE              = transpose(std(freq_mat_CO')/ sqrt(size(freq_mat_CO,2))); % Get SE
NEW_TBL2.Pw              = squeeze(mean(freq_mat_CO,2)); %  averageover participant
NEW_TBL2.Condition       = repmat({cfg.var_n}, length(freq_mat_CO), 1);
NEW_TBL2.Group           = repmat({"CO"}, length(freq_mat_CO),1);
NEW_TBL2.Lock            = repmat({cfg.lck}, length(freq_mat_CO),1);
NEW_TBL2.Ana_type        = repmat({cfg.Ana_type}, length(freq_mat_PD),1);
NEW_TBL2.Time            = transpose(cfg.Time);

% Finally join tables
NEW_TBL = [NEW_TBL1; NEW_TBL2];
