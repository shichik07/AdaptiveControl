%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Statistical Analysis using Fieldtrip           %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
LIMO                         = 'C:\Program Files\MATLAB\LIMO';
cm                           = 'C:\Program Files\MATLAB\FileExchangeAddOns'; % colormaps
% set directories
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% set new home directory for baseline correction
dirs.home                    = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData';
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData';

% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

addpath(LIMO)
addpath(eegl)
addpath(ft)
addpath(cm)
addpath(dirs.functions)





%% Load example data file for some useful information
examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_pl_RL_bslC_RL.mat"));
data                        = load(examp);
InfosRL.chanlocs                    = data.TF_phase.chanlocs;
InfosRL.freqs_base                  = data.TF_phase.Frequencies;
InfosRL.time                        = data.TF_phase.Time;
scale = 0.6; % scale for the single plots in dB

examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_pl_bslC"));
data                        = load(examp);
InfosSL.chanlocs                    = data.TF_phase.chanlocs;
InfosSL.freqs_base                       = data.TF_phase.Frequencies;
InfosSL.time                        = data.TF_phase.Time;




% participants with too little data
LWPC_missing = ["sub-PD_11", "sub-PD_21"];
ISPC_missing = ["sub-CO_6", "sub-PD_15", "sub-PD_21"];

% get layout
cfg           = [];
cfg.layout =  'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis\EasyCap128.mat';
cfg.feedback  = 'yes';
[EasyCap128, cfg] = ft_prepare_layout(cfg)

lock = ['SL', 'RL']; %response or stimulus locked

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
    
    freq_name = 'Theta';
    phase_type = 'pl';
    frequencies2ana         = [4, 8]; % we are interested 4 to 8 Hz
    
    for lc = 1: length(lock)
        
        if lc == 1
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_',eff_name ,'_', phase_type,'.mat']));
            load(load_loc, 'Betas');
            Infos = InfosSL;
            Order = [1, 2];
        elseif lc == 2
            % load test data set
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_',eff_name ,'_', phase_type,'_RL.mat']));
            load(load_loc, 'Betas');
            Infos = InfosRL;
            Order = [3, 4];
        end
        
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
                if lc == 1
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_',eff_name ,'_', phase_type,'.mat']));
                    load(load_loc, 'Betas');
                    
                elseif lc == 2
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_',eff_name ,'_', phase_type,'_RL.mat']));
                    load(load_loc, 'Betas');
                end
                
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
        MI_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + ( 0.5)*LW.Betas_Interaction;
        MI_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + (-0.5)*LW.Betas_Interaction;
        MC_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + (-0.5)*LW.Betas_Interaction;
        MC_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + ( 0.5)*LW.Betas_Interaction;
        
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
        
        
        % Get the data into fieldtrip format
        Group = LW.Group == 0;
        LWPC_MICon_CO = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
        LWPC_MCCon_CO = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);
        Group = LW.Group == 1;
        LWPC_MICon_PD = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
        LWPC_MCCon_PD = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);
        
        
        
        
        
        %% Fieldtrip Cluster-based Permutation Test statistics Loop through group
        for grp = 1:2
            if grp == 1
                varg1 = LWPC_MICon_CO;
                varg2 = LWPC_MCCon_CO;
                groupname = 'CO';
                pos = Order(1); 
            else
                varg1 = LWPC_MICon_PD;
                varg2 = LWPC_MCCon_PD;
                groupname = 'PD';
                pos = Order(2); 
            end
            
            
            cfg = [];
            
            if lc == 1
                cfg.latency          = [0 0.8];
            elseif lc == 2
                cfg.latency    = [-0.6 0.2];
            end
            
            cfg.channel          = {'EEG'};
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
            cfg.numrandomization = 1000;
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
                cfg = [];
                
                if lc == 1
                    cfg.toi    = [0.3:0.1:0.8];
                elseif lc == 2
                    cfg.toi    = [-0.600:0.075:0.100];
                end
                    
                cfg.highlightsymbolseries = ['*','*'];
                cfg.markersize = [20 20];
                cfg.alpha  = 0.025;
                cfg.subplotsize = [2 5];
                cfg.highlightsizeseries = [6 6 6 6 6];
                cfg.highlightfontsizeseries = [200 200 200 200 200];
                cfg.highlightfontsize = [20]
                cfg.parameter = 'stat';
                cfg.zlim   = [-4 4];
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
                
                close(f1)
                %% Next let us plot the cluster over the active period only
                % first we need to identiffy electrodes and timepoints tha
                % are significant. First for the negative then the positive
                % clusters.
                % note, really crappy function - be careful if you have
                % more than one cluster
                
                [min_ind, max_ind, timewin_cl, pos_ele, neg_ele, p_pos, p_neg] = find_cluster(stat, 0.025);
                times2plot = [0.6];
                %%get index of time
                t_ind = Infos.time >= timewin_cl(1) & Infos.time <= timewin_cl(end);
                %[~,t_ind]  = (min(abs(Infos.time - times2plot)));
                %index group to plot
                grp_ind = LW.Group == grp-1;
                
                % save the negative cluster, so that we can use it
                % elsewehere as a table
                Cluster_elec = table(neg_ele);
                Cluster_time = table(timewin_cl'); % transpose to get it into a column
                save_name_elec = ['Cluster_RL_', groupname ,'_', eff_name,'_', freq_name,'_electrodes.csv'];
                save_name_time = ['Cluster_RL_', groupname ,'_', eff_name,'_', freq_name,'_time.csv'];
                save_folder = 'E:\AdaptiveControl\Data\TF_for_R';
                save_loca_elec = fullfile(save_folder,save_name_elec);
                save_loca_time = fullfile(save_folder,save_name_time);
                
                % time indices to plot
                freqdata2plot = double(squeeze(mean(M_Stroop_dif(:, frequencies2ana, t_ind, grp_ind),4))); % average over participants
                freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
                freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
                
                
                
                if (~isempty(pos_ele) + ~isempty(neg_ele)) == 2
                    warning('We have both postive and negative clusters')
                    break
                end
                
               
                plt1 = subplot(1,4,pos);
                
                if ~isempty(pos_ele) == 1
                    % find index of labeled electrodes
                    List_one = {Infos.chanlocs.labels};
                    List_two = pos_ele';
                    idx=find(ismember(List_one,List_two))
                    %Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-scale scale], 'electrodes', 'on', 'emarker', {'.','k',{},4}, 'emarker2', {idx,'+','k', 6, 2});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    
                elseif  ~isempty(neg_ele) == 1
                    % find index of labeled electrodes
                    List_one = {Infos.chanlocs.labels};
                    List_two = neg_ele';
                    idx=find(ismember(List_one,List_two))
                    % Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-scale scale], 'electrodes', 'on', 'emarker', {'.','k',{},4}, 'emarker2', {idx,'*','k', 6, 1});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    
                elseif  isempty(neg_ele) == 1 & isempty(pos_ele) == 1
                    t_ind = Infos.time >= -0.6 & Infos.time <= -0.2;
                    % time indices to plot
                    freqdata2plot = double(squeeze(mean(M_Stroop_dif(:, frequencies2ana, t_ind, grp_ind),4))); % average over participants
                    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
                    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
                    
                    % Plot Topo
                    f2 = figure
                    topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-scale scale], 'electrodes', 'on', 'emarker', {'.','k',{},4});
                    colormap(flipud(brewermap([], "Rdbu")));
                    h = colorbar;
                    
                    
                end
                
                
            catch
                warning('no clusters present with a p-value lower than the specified alpha, nothing to plot');
                plt1 = subplot(1,4,pos);
                % But that doesnt mean we will not plot for the PD partcipants! We load the cluster
                % electrodes of the CO Participants and get ourslef the time
                % window to
                t_ind = Infos.time >= -0.6 & Infos.time <= -0.2;
                freqdata2plot = double(squeeze(mean(M_Stroop_dif(:, frequencies2ana, t_ind, grp_ind),4))); % average over participants
                freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
                freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
                
                List_one = {Infos.chanlocs.labels};
                List_two = pos_ele';
                idx=find(ismember(List_one,List_two))
                % Plot Topo
                topoplot(freqdata2plot,Infos.chanlocs,'maplimits',[-scale scale], 'electrodes', 'on', 'emarker', {'.','k',{},4}, 'emarker2', {idx,'+','k', 6, 2});
            end
            
        end
    end
    % get colorbar
        colormap(flipud(brewermap([], "Rdbu")));
        h = colorbar;
        set(h,'Position',[0.9 0.05 0.03 0.9])

        f = gcf;
        set(gcf, 'Units', 'centimeters', 'Position', [0 0 16 8]) % set to have size 16 cm width
        save_n = ['BlockCongruence_Cluster_' + eff_name + '_' + lock + '_' + phase_type + 'new.png'];
        
        save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\fin_dataset\RL_pl';
        save_loc = fullfile(save_fl,save_n);
        exportgraphics(f, save_loc, 'Resolution', 300)
        saveas(f, fullfile(save_fl,save_fign))
end