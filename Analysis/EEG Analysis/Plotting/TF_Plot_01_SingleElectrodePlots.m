%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Single Electrode Plotting                      %
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
examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_pl_bslC"));
data                        = load(examp);
Infos.chanlocs                    = data.TF_phase.chanlocs;
channel_l = {Infos.chanlocs(:).labels};
Infos.freqs_base                       = data.TF_phase.Frequencies;
freqs = Infos.freqs_base        ;
Infos.time                        = data.TF_phase.Time;
Time_SL = Infos.time;
cont_n = 15; % number of contours for contourblock

%% load RL Time
examp                       = fullfile(dirs.eegsave, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_pl_RL_bslC_RL.mat"));
data                        = load(examp);
Infos.time                        = data.TF_phase.Time;
Time_RL = Infos.time;
%participants with too little data for a GLM
LWPC_missing = ["sub-PD_11", "sub-PD_21"];
ISPC_missing = ["sub-CO_6", "sub-PD_15", "sub-PD_21"];

ROI_elect = ["FCz"]; %, "Cz", "CPz"];
freq_low = 2; %2Hz
freq_up = 30; %30Hz
max_scale = 0.8;
min_scale = -0.8;

locks = ["SL", "RL"];


p_type = 2
phase_type = 'pl';

for eff = 1:2
    if eff == 1
        eff_name = "ISPC";
        umgsprl  = "Reactive";
        missing = ISPC_missing;
    else
        eff_name = "LWPC";
        umgsprl  = "Proactive";
        missing = LWPC_missing;
    end
    clf
    for lc = 1:length(locks)
        lock = locks(lc);
        if lc == 1
            Time = Time_SL;
            order = [1, 2, 5, 6];
        elseif lc == 2
            Time = Time_RL;
            order = [3, 4, 7, 8];
        end
        % load test data set
        if lc == 1
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
        elseif lc == 2
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_' + eff_name + '_' + phase_type + '_RL.mat']));
        end
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
                if lc == 1
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
                elseif lc == 2
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_' + eff_name + '_' + phase_type + '_RL.mat']));
                end
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
        MI_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + ( 0.5)*LW.Betas_Interaction;
        MI_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + (-0.5)*LW.Betas_Interaction;
        MC_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + (-0.5)*LW.Betas_Interaction;
        MC_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + ( 0.5)*LW.Betas_Interaction;
        
        % Caluclate Congruency Effects per block
        MI_Con = MI_I - MI_C;
        MC_Con = MC_I - MC_C;
        
        % Calculate Congruency effects
        Congruency_effect = (MI_Con + MC_Con)/.2;
        
        % Interaction
        M_Inter = MI_Con - MC_Con;
        
        
        %% Plot single electrode Power spectrum
        for elect = 1:length(ROI_elect)
            %electrode of interest
            ch_t = find(strcmp(channel_l, ROI_elect(elect)));
            
            %% Cognruency Effect
            
            %limits congruency effect
            max_scale_con = 10;
            min_scale_con = -10;
            
            
            %% Plot Congruence effects per Block
            %limits congruency effect
            
            max_scale_blc = 1.2;
            min_scale_blc = -1.2;
            
            
            
            % Plot Congruence
            fig3 = figure(3)
            plt1 = subplot(2,4,order(1));
            
            
            %Plot CO group
            Group = LW.Group == 0;
            data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('HC - Group: MI Block')
            
            % slightly shift the plots
            pos = get(plt1, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt1, 'Position', posnew)
            
            % Plot PD Group
            
            plt2 = subplot(2,4,order(2));
            Group = LW.Group == 1;
            
            data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('PD - Group: MI Block')
            
            % slightly shift the plots
            pos = get(plt2, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt2, 'Position', posnew)
            
            plt3 = subplot(2,4,order(3));
            
            %Plot CO group
            Group = LW.Group == 0;
            data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('HC - Group: MC Block')
            
            % slightly shift the plots
            pos = get(plt3, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt3, 'Position', posnew)
            
            plt4 = subplot(2,4,order(4));
            
            % Plot PD Group
            Group = LW.Group == 1;
            data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('PD - Group: MC Block')
            
            % slightly shift the plots
            pos = get(plt4, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt4, 'Position', posnew)
            
        end
        
    end
    % get colorbar
    colormap(flipud(brewermap([], "Rdbu")));
    h = colorbar;
    set(h,'Position',[0.9 0.05 0.03 0.9])
    h.FontSize = 10;
    h.FontWeight = 'bold';
    t_input = 'Conflict ' + eff_name + ' Effect: ' + ROI_elect(elect);
    sgt = suptitle([ t_input]);%sgtitle([ t_input]);
    sgt.FontSize = 15;
    sgt.FontWeight = 'bold';
    
    % labels
    han=axes(fig3,'visible','off');
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Frequency in Hz', 'Position',[-0.12   0.46]);
    xlabel(han,'Time in s', 'Position',[0.43   -0.06]);
    han.FontWeight = 'bold';
    
    
    
    f = gcf;
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 16 8]) % set to have size 16 cm width
    save_n = ['BlockCongruence_difference_' + ROI_elect(elect)+ '_' + eff_name + '_' + lock + '_' + phase_type + 'new.png'];
    save_fign = ['BlockCongruence_difference_' + ROI_elect(elect)+ '_' + eff_name + '_' + lock + '_' + phase_type + 'new.fig'];
    
    save_fl = fullfile('C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\SingleElectrodes\', eff_name, lock);
    save_loc = fullfile(save_fl,save_n);
    exportgraphics(f, save_loc, 'Resolution', 300)
    saveas(f, fullfile(save_fl,save_fign))
end
%% Next let us extract the Infos for the individual timelines by electrodes for the pl data only
frequencies2ana         = [4, 8]; % we are interested 4 to 8 Hz
%elec2ana = {'FCz', 'FC1', 'FC2', 'FCC1h', 'FCC2h', 'FFC1h', 'FFC2h'};
elec2ana = {'Cz'}%, 'C1', 'C2'};



for lc = 1:length(locks)
    lock = locks(lc);
    if lc == 1
        
        Time = Time_SL;
    elseif lc == 2
       
        Time = Time_RL;
    end
    
    phase_type = 'pl';
    
    for eff = 1:2
        if eff == 1
            eff_name = "ISPC";
            umgsprl  = "Reactive";
            missing = ISPC_missing;
        else
            eff_name = "LWPC";
            umgsprl  = "Proactive";
            missing = LWPC_missing;
        end
        
        % load test data set
        if lc == 1
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
        elseif lc == 2
            load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_' + eff_name + '_' + phase_type + '_RL.mat']));
        end
        load(load_loc, 'Betas');
        
        Betas = reshape_3d_limo_fixed(Betas, Infos.freqs_base, Infos.time);
        Beta_s = size(Betas);
        
        %initialze data matrices
        LW.Betas_intercept              = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Congruency             = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Block                  = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Betas_Interaction            = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
        LW.Group                        = NaN(Part_N,1);
        LW.Participant                  = cell(Part_N,1);
        
        % Load the weight matrices and get the data
        for sub = 2:Part_N
            if any(matches(missing, Participant_IDs{sub}))
                continue
            else
                if lc == 1
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
                elseif lc == 2
                    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_' + eff_name + '_' + phase_type + '_RL.mat']));
                end
                load(load_loc, 'Betas');
                
                Betas                           = reshape_3d_limo_fixed(Betas, Infos.freqs_base, Infos.time);
                LW.Betas_intercept(:,:,:,sub)   = Betas(:,:,:, 4);
                LW.Betas_Congruency(:,:,:,sub)  = Betas(:,:,:, 1);
                LW.Betas_Block(:,:,:,sub)       = Betas(:,:,:, 2);
                LW.Betas_Interaction(:,:,:,sub) = Betas(:,:,:, 3);
                LW.Participant{sub}             = Participant_IDs{sub};
                
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
        
        % Average over 4-8Hz to make it one frequency
        freqsana         = Infos.freqs_base >= frequencies2ana(1) & Infos.freqs_base <= frequencies2ana(2);
        
        
        %plot all conditions
        % get electrodes
        List_one = {Infos.chanlocs.labels};
        idx=find(ismember(List_one,elec2ana));
        
        
        % prepare information to be contained in tables
        cfg = [];
        cfg.Time = Time; % time
        cfg.lck = lock; % response locked or stimulus locked
        cfg.grp_idx = LW.Group; % index info for groups
        cfg.Ana_type = eff_name; % effect we want to anlyze
        cfg.freq_idx = freqsana; % indexes we want to summarize
        cfg.elec_idx = idx; % index of electrodes we preselected (now legacy)
        cfg.chanlocs =  {Infos.chanlocs.labels}; % channel labels
        
        % Get the data into table form
        cfg.var_n = "MI_C";
        MI_C_T = theta_mat_to_R(MI_C , cfg);
        cfg.var_n = "MI_I";
        MI_I_T = theta_mat_to_R(MI_I , cfg);
        cfg.var_n = "MC_C";
        MC_C_T = theta_mat_to_R(MC_C , cfg);
        cfg.var_n = "MC_I";
        MC_I_T = theta_mat_to_R(MC_I , cfg);
        
        % concatenate tables
        Theta_table = [MI_C_T; MI_I_T; MC_C_T; MC_I_T];
        
        % save table
        t_save_n = ["Theta_data" + "_" + eff_name + "_" + lock + ".csv"];
        save_fl = 'E:\AdaptiveControl\Data\TF_for_R';
        save_loc = fullfile(save_fl,t_save_n);
        writetable(Theta_table,save_loc)
        
        % Now do the same to get the tables that you want to use for the
        % brain behavior correlations
        
        % prepare information to be contained in tables
        cfg             = [];
        cfg.Time        = Time; % time
        cfg.lck         = lock; % response locked or stimulus locked
        cfg.grp_idx     = LW.Group; % index info for groups
        cfg.Ana_type    = eff_name; % effect we want to anlyze
        cfg.freq_idx    = freqsana; % indexes we want to summarize
        cfg.chanlocs    = {Infos.chanlocs.labels}; % channel labels
        cfg.part_name   = LW.Participant; % participant name info
        
        cfg.var_n       = "Congruent";
        Diff_con        = theta_average_to_R(MI_C, MC_C, cfg);
        cfg.var_n       = "Incongruent";
        Diff_incon      = theta_average_to_R(MI_I, MC_I, cfg);
        
        % concatenate tables
        Theta_diff_table = [Diff_con; Diff_incon];
        
        % save table
        t_save_n = ["Theta_cluster_averaged" + "_" + eff_name + "_" + lock + ".csv"];
        save_fl = 'E:\AdaptiveControl\Data\TF_for_R';
        save_loc = fullfile(save_fl,t_save_n);
        writetable(Theta_diff_table,save_loc)
        
    end
end

