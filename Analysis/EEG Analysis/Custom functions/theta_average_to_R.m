function TBL_brain_behav =  theta_average_to_R(MC, MI,  cfg)

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
        save_name_CO_t = ['Cluster_CO_' + cfg.Ana_type + '_Theta_time.csv'];
        save_name_PD_t = ['Cluster_PD_' + cfg.Ana_type + '_Theta_time.csv'];
    elseif cfg.lck == "RL"
        save_name_CO = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_RL_PD_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_CO_t = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_time.csv'];
        save_name_PD_t = ['Cluster_RL_PD_' + cfg.Ana_type + '_Theta_time.csv'];
    end
% for the ISPC effect we only load the CO cluster, because there was no significant cluster in the PD group 
elseif cfg.Ana_type == "ISPC"
    if cfg.lck =="SL"
        save_name_CO = ['Cluster_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_CO_t = ['Cluster_CO_' + cfg.Ana_type + '_Theta_time.csv'];
        save_name_PD_t = ['Cluster_CO_' + cfg.Ana_type + '_Theta_time.csv'];
    elseif cfg.lck == "RL"
        save_name_CO = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_PD = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_electrodes.csv'];
        save_name_CO_t = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_time.csv'];
        save_name_PD_t = ['Cluster_RL_CO_' + cfg.Ana_type + '_Theta_time.csv'];
    end
end

% load the electrode cluster data
save_loc_PD = fullfile(save_folder,save_name_PD);
cluster_PD = readtable(save_loc_PD);
elec_idx_PD = find(ismember(cfg.chanlocs, table2cell(cluster_PD)'));
save_loc_CO = fullfile(save_folder,save_name_CO);
cluster_CO = readtable(save_loc_CO);
elec_idx_CO = find(ismember(cfg.chanlocs, table2cell(cluster_CO)'));

% load the time cluster data
save_loc_PD_t = fullfile(save_folder,save_name_PD_t);
t_PD = readtable(save_loc_PD_t);
t_idx_PD = find(ismember(cfg.Time', table2array(t_PD)));
save_loc_CO_t = fullfile(save_folder,save_name_CO_t);
t_CO = readtable(save_loc_CO_t);
t_idx_CO = find(ismember(cfg.Time', table2array(t_CO)));

% get correct electrode indices

% get table
TBL_brain_behav1 = table;
%% Because we are interested in Brain - Behavior Correlations, we als want a table that summarizes the activity for every participant
% in the period in question

% get the PD data
freq_mat_MI_PD              = double(squeeze(mean(MI(elec_idx_PD, cfg.freq_idx,t_idx_PD, grp_PD),1))); % average over electrodes
freq_mat_MC_PD              = double(squeeze(mean(MC(elec_idx_PD, cfg.freq_idx,t_idx_PD, grp_PD),1))); % average over electrodes
freq_mat_MI_PD              = double(squeeze(mean(freq_mat_MI_PD ,1))); % average over frequencies
freq_mat_MC_PD              = double(squeeze(mean(freq_mat_MC_PD ,1))); % average over frequencies
TBL_brain_behav1.Pw            = [double(squeeze(mean(freq_mat_MI_PD ,1))) - double(squeeze(mean(freq_mat_MC_PD ,1)))]'; % average over time and subtract MC from MI
TBL_brain_behav1.Group           = repmat("PD", sum(grp_PD), 1);
TBL_brain_behav1.Condition       = repmat({cfg.var_n}, sum(grp_PD),1);
TBL_brain_behav1.Lock            = repmat({cfg.lck}, sum(grp_PD),1);
TBL_brain_behav1.Ana_type        = repmat({cfg.Ana_type}, sum(grp_PD),1);
TBL_brain_behav1.Part            = cfg.part_name(grp_PD);

TBL_brain_behav2 = table;
% get the CO data
freq_mat_MI_CO              = double(squeeze(mean(MI(elec_idx_CO, cfg.freq_idx,t_idx_CO, grp_CO),1))); % average over electrodes
freq_mat_MC_CO              = double(squeeze(mean(MC(elec_idx_CO, cfg.freq_idx,t_idx_CO, grp_CO),1))); % average over electrodes
freq_mat_MI_CO              = double(squeeze(mean(freq_mat_MI_CO ,1))); % average over frequencies
freq_mat_MC_CO              = double(squeeze(mean(freq_mat_MC_CO ,1))); % average over frequencies
TBL_brain_behav2.Pw            = [double(squeeze(mean(freq_mat_MI_CO ,1))) - double(squeeze(mean(freq_mat_MC_CO ,1)))]'; % average over time and subtract MC from MI
TBL_brain_behav2.Group           = repmat("CO", sum(grp_CO), 1);
TBL_brain_behav2.Condition       = repmat({cfg.var_n}, sum(grp_CO),1);
TBL_brain_behav2.Lock            = repmat({cfg.lck}, sum(grp_CO),1);
TBL_brain_behav2.Ana_type        = repmat({cfg.Ana_type}, sum(grp_CO),1);
TBL_brain_behav2.Part            = cfg.part_name(grp_CO);


TBL_brain_behav = [TBL_brain_behav1; TBL_brain_behav2];
