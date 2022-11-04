%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Statistical Analysis using Fieldtrip           %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Explore beta distribution
clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
LIMO                         = 'C:\Program Files\MATLAB\LIMO';
cm                           = 'C:\Program Files\MATLAB\FileExchangeAddOns'; % colormaps
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_9'})); %remove PD 9 bad data
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

addpath(LIMO)
addpath(eegl)
addpath(ft)
addpath(cm)
addpath(dirs.functions)



%% Load example data file for some useful information
examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_npl_bslC"));
data                        = load(examp, 'TF_non_phase');
Infos.chanlocs                    = data.TF_non_phase.chanlocs;
Infos.freqs_base                       = data.TF_non_phase.Frequencies;
Infos.time                        = data.TF_non_phase.Time;

%participants with too little data for a GLM
LWPC_missing = ["sub-PD_11", "sub-PD_21"];
ISPC_missing = ["sub-CO_6", "sub-PD_15", "sub-PD_21"];

% get layout
cfg           = [];
cfg.layout =  'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis\EasyCap128.mat';
cfg.feedback  = 'yes';
[EasyCap128, cfg] = ft_prepare_layout(cfg)

for eff = 1:2
    if eff == 1
        eff_name = 'ISPC';
        umgsprl  = 'Reactive';
        missing = ISPC_missing;
    else
        eff_name = 'LWPC';
        umgsprl  = 'Proactive';
        missing = LWPC_missing;
    end
    for freq_type = 1:2
        if freq_type == 1
            freq_name = 'Theta';
            phase_type = 'pl';
            frequencies2ana         = [4, 8]; % we are interested 4 to 8 Hz
        else
            freq_name = 'Alpha';
            phase_type = 'pl';
            frequencies2ana         = [8, 13]; % we are interested 4 to 8 Hz
        end
        % load test data set
        load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_',eff_name ,'_', phase_type,'.mat']));
        load(load_loc, 'Betas');
        
        Betas = reshape_3d_limo_fixed(Betas, Infos.freqs_base, Infos.time);
        
        Beta_s = size(Betas);
        
        
        
        %initialze data matrices
        LW.Betas_intercept              = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Congruency             = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Block                  = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Interaction            = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Group                        = NaN(Part_N,1);
        
        
        % Load the weight matrices and get the data
        for sub = 2:Part_N
            if any(matches(missing, Participant_IDs{sub}))
                continue
            else
                load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_',eff_name ,'_', phase_type,'.mat']));
                load(load_loc, 'Betas');
                
                Betas                           = reshape_3d_limo_fixed(Betas, Infos.freqs_base, Infos.time);
                LW.Betas_intercept(:,:,:,sub)   = Betas(:,:,:, 4);
                LW.Betas_Congruency(:,:,:,sub)  = Betas(:,:,:, 1);
                LW.Betas_Block(:,:,:,sub)       = Betas(:,:,:, 2);
                LW.Betas_Interaction(:,:,:,sub) = Betas(:,:,:, 3);
                
                % get index variable
                if Participant_IDs{sub}(5:6) == 'CO'
                    LW.Group(sub) = 0;
                elseif Participant_IDs{sub}(5:6) == 'PD'
                    LW.Group(sub) = 1;
                end
                clear Betas
            end
        end
        
        % calculate marginal means
        MI_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + ( 0.25)*LW.Betas_Interaction;
        MI_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + (-0.25)*LW.Betas_Interaction;
        MC_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + (-0.25)*LW.Betas_Interaction;
        MC_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + ( 0.25)*LW.Betas_Interaction;
        
        % Caluclate Congruency Effects per block
        MI_Con = MI_I - MI_C;
        MC_Con = MC_I - MC_C;
        
        M_Stroop_dif = MI_Con - MC_Con;
        
        % Average over 4-8Hz to make it one frequency
        frequencies2ana         = Infos.freqs_base >= frequencies2ana(1) & Infos.freqs_base <= frequencies2ana(2);
        Infos.freqs             = Infos.freqs_base; %so we do not always add to it
        Infos.freqs(end+1)      = 1; % we add another frequency as a standin for our averaged theta data
        
        MI_Con(:, end+1, :, :)  = squeeze(mean(MI_Con(:, frequencies2ana, :, :),2));
        MC_Con(:, end+1, :, :)  = squeeze(mean(MC_Con(:, frequencies2ana, :, :),2));
        % MC_C(:, end+1, :, :)  = squeeze(mean(MC_C(:, frequencies2ana, :, :),2));
        % MC_I(:, end+1, :, :)  = squeeze(mean(MC_I(:, frequencies2ana, :, :),2));
        
        
        
        % Get the data into fieldtrip format
        Group = LW.Group == 0;
        LWPC_MICon_CO = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
        LWPC_MCCon_CO = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);
        Group = LW.Group == 1;
        LWPC_MICon_PD = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
        LWPC_MCCon_PD = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);
        
        % %Test
        % Group = LW.Group == 1;
        % LWPC_MI_PD = get_pseudo_FT_struct_cluster(MC_I(:,:,:,Group), Infos);
        % LWPC_MC_PD = get_pseudo_FT_struct_cluster(MC_C(:,:,:,Group), Infos);
        
        
        
        
        
        
        %% Fieldtrip Cluster-based Permutation Test statistics Loop through group
        for grp = 1:2
            if grp == 1
                varg1 = LWPC_MICon_CO;
                varg2 = LWPC_MCCon_CO;
                groupname = 'CO'
            else
                varg1 = LWPC_MICon_PD;
                varg2 = LWPC_MCCon_PD;
                groupname = 'PD'
            end
            
            
            cfg = [];
            cfg.channel          = {'EEG'};
            cfg.latency          = [0 0.8];
            cfg.frequency        = 1; % one the standin for our theta averaged data
            cfg.method           = 'montecarlo';
            cfg.statistic        = 'ft_statfun_depsamplesT';
            cfg.correctm         = 'cluster';
            cfg.clusteralpha     = 0.05;
            cfg.clusterstatistic = 'maxsum';
            cfg.minnbchan        = 2;
            cfg.tail             = 0;
            cfg.clustertail      = 0;
            cfg.alpha            = 0.05;
            cfg.numrandomization = 10000;
            % specifies with which sensors other sensors can form clusters
            cfg_neighb.method    = 'distance';
            cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, varg1);
            
            
            subj = size(varg1.powspctrm, 1);
            design = zeros(2,2*subj);
            for i = 1:subj
                design(1,i) = i;
            end
            for i = 1:subj
                design(1,subj+i) = i;
            end
            design(2,1:subj)        = 1;
            design(2,subj+1:2*subj) = 2;
            
            cfg.design   = design;
            cfg.uvar     = 1;
            cfg.ivar     = 2;
            
            [stat] = ft_freqstatistics(cfg,varg1, varg2)
            
            try
                % plot results of the analysis
                cfg = [];
                cfg.highlightsymbolseries = ['*','*'];
                cfg.markersize = [20 20];
                cfg.alpha  = 0.025;
                cfg.subplotsize = [2 5];
                cfg.highlightsizeseries = [6 6 6 6 6];
                cfg.highlightfontsizeseries = [200 200 200 200 200];
                cfg.highlightfontsize = [20]
                cfg.parameter = 'stat';
                cfg.zlim   = [-4 4];
                cfg.toi    = [0.3:0.1:0.8];
                cfg.layout = 'easycapM15.mat';
                f1 = figure
                plt = ft_clusterplot_mod(cfg, stat);
                
                colormap(flipud(brewermap([], "Rdbu")));
                h = colorbar;
                set(h,'Position',[0.93 0.05 0.03 0.9])
                h.FontSize = 10;
                h.FontWeight = 'bold';
                sgt = sgtitle([freq_name, ' Cluster ', eff_name,' Effect: ', groupname ,' Group']);
                sgt.FontSize = 15;
                sgt.FontWeight = 'bold'
                
                f = gcf;
                
                save_n = ['Cluster_', groupname ,'_', eff_name,'_', freq_name,'_interaction_StimL_pl.png'];
                save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\fin_dataset\SL_npl';
                save_loc = fullfile(save_fl,save_n);
                exportgraphics(f, save_loc, 'Resolution', 300)
                
                close(f1)
                
                
                %% Next let us plot the cluster over the active period only 
                % first we need to identiffy electrodes and timepoints tha
                % are significant. First for the negative then the positive
                % clusters.
                % note, really crappy function - be careful if you have
                % more than one cluster
                
                [min_ind, max_ind, timewin_cl, pos_ele, neg_ele, p_pos, p_neg] = find_cluster(stat, 0.025);
                %times2plot = [0.6];
                %%get index of time
                t_ind = Infos.time >= timewin_cl(1) & Infos.time <= timewin_cl(end);
                %t_ind = Infos.time >= times2plot(1) & Infos.time <= times2plot(2);
                %[~,t_ind]  = (min(abs(Infos.time - times2plot)));
                %index group to plot
                grp_ind = LW.Group == grp-1;
                
                % time indices to plot
                freqdata2plot = double(squeeze(mean(M_Stroop_dif(:, frequencies2ana, t_ind, grp_ind),4))); % average over participants
                freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
                freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
                
               
                
                if (~isempty(pos_ele) + ~isempty(neg_ele)) == 2
                    warning('We have both postive and negative clusters')
                    break
                end
                
                if ~isempty(pos_ele) == 1
                    % find index of labeled electrodes
                    List_one = {Infos.chanlocs.labels};
                    List_two = pos_ele';
                    idx=find(ismember(List_one,List_two))
                    % Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-0.3 0.3], 'electrodes', 'on', 'emarker', {'.','k',{},4}, 'emarker2', {idx,'+','k', 6, 2});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    txt =  num2str(timewin_cl(1)) + "s to " + num2str(timewin_cl(end)) + "s; p = " + num2str(p_pos);
                    text(-0.3, -0.6,txt,'FontSize',12, 'FontWeight', 'bold')
                
                elseif  ~isempty(neg_ele) == 1
                    % find index of labeled electrodes
                    List_one = {Infos.chanlocs.labels};
                    List_two = neg_ele';
                    idx=find(ismember(List_one,List_two))
                    % Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-0.3 0.3], 'electrodes', 'on', 'emarker', {'.','k',{},4}, 'emarker2', {idx,'*','k', 6, 1});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    txt =  num2str(timewin_cl(1)) + "s to " + num2str(timewin_cl(end)) + "s; p = " + num2str(p_neg);
                    text(-0.3, -0.6,txt,'FontSize',12, 'FontWeight', 'bold')
                    
                elseif  isempty(neg_ele) == 1 & isempty(pos_ele) == 1
                    t_ind = Infos.time >= 0.3 & Infos.time <= 0.7;
                    %get data
                    % time indices to plot
                    freqdata2plot = double(squeeze(mean(M_Stroop_dif(:, frequencies2ana, t_ind, grp_ind),4))); % average over participants
                    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
                    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
                    
                    % Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-0.3 0.3], 'electrodes', 'on', 'emarker', {'.','k',{},4});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    txt = "0.3s to 0.7s; p > 0.05 " ;
                    text(-0.3, -0.6,txt,'FontSize',12, 'FontWeight', 'bold')
                
                end
                
             
                set(h,'Position',[0.9 0.05 0.03 0.9])
                h.FontSize = 10;
                h.FontWeight = 'bold';
                sgt = sgtitle([freq_name, ' Cluster ', umgsprl,' Control: ', groupname ,' Group']);
                sgt.FontSize = 15;
                sgt.FontWeight = 'bold';
               
                f = gcf;
                
                save_n = ['Cluster_', groupname ,'_', eff_name,'_', freq_name,'_interaction_summary_StimL_pl.png'];
                save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\fin_dataset\SL_npl';
                save_loc = fullfile(save_fl,save_n);
                exportgraphics(f, save_loc, 'Resolution', 300)
                close(f2)
                
                
                
               
    
                
            catch
                warning('no clusters present with a p-value lower than the specified alpha, nothing to plot');
                
            end
            
        end
    end
end

% Get cluster
cluster = squeeze(stat.negclusterslabelmat) == 1;
ind = zeros(size(cluster,1),1);
for i = 1:length(ind)
    if any(cluster(i,:) == 1)
        ind(i) = 1;
    end
end
ind = find(ind);
Labels = length(stat.label(ind));

start = 0;
stop = 0;
for ii = 1:size(cluster,2)
    if any(cluster(:,ii) == 1)
        if start == 0
            start = ii;
        end
    elseif start ~= 0
        if stop == 0
            stop = ii;
        end
    end
end
    
CO_start = start;
CO_stop = stop;


PD_start = start;
PD_stop = stop;



CO_MCCon = squeeze(LWPC_MCCon_CO.powspctrm(:, ind_CO, 21, :));
CO_MICon = squeeze(LWPC_MICon_CO.powspctrm(:, ind_CO, 21, :));
CO_ClusterSize = length(ind_CO);
CO_sigDif = [stat.time(CO_start), stat.time(CO_stop)];

PD_MCCon = squeeze(LWPC_MCCon_PD.powspctrm(:, ind_PD, 21, :));
PD_MICon = squeeze(LWPC_MICon_PD.powspctrm(:, ind_PD, 21, :));
PD_ClusterSize = length(ind_PD);
PD_sigDif = [stat.time(PD_start), stat.time(PD_stop)];% Cluster looked only between 0 and 0.8 s, so we need to recalculate time


CO_MCCon = squeeze(mean(CO_MCCon, 2));
CO_MICon = squeeze(mean(CO_MICon, 2));
PD_MCCon = squeeze(mean(PD_MCCon, 2));
PD_MICon = squeeze(mean(PD_MICon, 2));

Time = LWPC_MICon_PD.time;

save_n = 'Cluster_data_theta.mat';
save_loc= fullfile('E:\AdaptiveControl\Data', 'TF_for_R', save_n);
save(save_loc,'Time', 'PD_MCCon', 'PD_MICon', 'CO_MCCon', 'CO_MICon', ...
    'CO_ClusterSize', 'PD_ClusterSize','CO_sigDif', 'PD_sigDif');