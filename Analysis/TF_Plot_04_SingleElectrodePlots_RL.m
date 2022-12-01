clear all; close all; clc
dbstop if error

spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
LIMO                         = 'C:\Program Files\MATLAB\LIMO';
cm                           = 'C:\Program Files\MATLAB\FileExchangeAddOns'; % colormaps
% set directories
dirs.home                    = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
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
channel_l = {Infos.chanlocs(:).labels};
Infos.freqs_base                       = data.TF_non_phase.Frequencies;
freqs = Infos.freqs_base        ;
Infos.time                        = data.TF_non_phase.Time;
Time = Infos.time;

%participants with too little data for a GLM
LWPC_missing = ["sub-PD_11", "sub-PD_21"];
ISPC_missing = ["sub-CO_6", "sub-PD_15", "sub-PD_21"];

% get layout
cfg           = [];
cfg.layout =  'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis\EasyCap128.mat';
cfg.feedback  = 'yes';
[EasyCap128, cfg] = ft_prepare_layout(cfg)


ROI_elect = ["FCz", "Cz", "CPz"];
freq_low = 2; %2Hz
freq_up = 30; %30Hz
max_scale = 0.8;
min_scale = -0.8;

lock = "RL";


for p_type = 1:2
    if p_type == 1
        phase_type = 'npl';
    elseif p_type == 2
        phase_type = 'pl';
    end
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
        load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
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
                load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, ['_GLM_' + eff_name + '_' + phase_type + '.mat']));
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
        
        % Calculate Congruency effects
        Congruency_effect = (MI_Con + MC_Con)/.2;
        
        % Interaction 
        M_Inter = MC_Con - MI_Con;

        
        
        
        %% Plot single electrode Power spectrum
        for grp = 1:2
            if grp == 1
                groupname = 'CO'
            else
                groupname = 'PD'
            end
            for elect = 1:length(ROI_elect)
                %electrode of interest
                ch_t = find(strcmp(channel_l, ROI_elect(elect)));
            
                %% Cognruency Effect
               
                %limits congruency effect
                max_scale_con = 10;
                min_scale_con = -10;
                
                % Plot Congruence
                figure(1)
                plt = tiledlayout(1,2);
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = squeeze(mean(Congruency_effect(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_con, max_scale_con],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group')
                
                % Plot PD Group
                Group = LW.Group == 1;
                nexttile
                data_C = squeeze(mean(Congruency_effect(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_con, max_scale_con],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('PD - Group')
                
                colormap(flipud(brewermap([], "Rdbu")));
                h = colorbar;
                h.FontSize = 10;
                h.FontWeight = 'bold';
                t_input = 'Conflict ' + umgsprl + ' Control: ' + ROI_elect(elect) + ' (' + lock + ')' ;
                sgt = sgtitle([ t_input]);
                sgt.FontSize = 15;
                sgt.FontWeight = 'bold'
                
                f = gcf;
                save_n = ['Congruency_effect_' + ROI_elect(elect)+ '_' + eff_name + '_' + lock + '.png'];
                save_fl = fullfile('C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\SingleElectrodes\', eff_name, lock);
                save_loc = fullfile(save_fl,save_n);
                exportgraphics(f, save_loc, 'Resolution', 300)
                
                %% Plot Interaction
                %limits congruency effect
                max_scale_int = 0.50;
                min_scale_int = -0.50;
                
                % Plot Congruence
                figure(2)
                plt = tiledlayout(1,2);
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_int, max_scale_int],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group')
                
                % Plot PD Group
                Group = LW.Group == 1;
                nexttile
                data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_int, max_scale_int],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('PD - Group')
                
                colormap(flipud(brewermap([], "Rdbu")));
                h = colorbar;
                h.FontSize = 10;
                h.FontWeight = 'bold';
                t_input = 'Conflict interaction ' + umgsprl + ' Control: ' + ROI_elect(elect) + ' (' + lock + ')' ;
                sgt = sgtitle([ t_input]);
                sgt.FontSize = 15;
                sgt.FontWeight = 'bold'
                
                f = gcf;
                save_n = ['Interaction_effect_' + ROI_elect(elect)+ '_' + eff_name + '_' + lock + '.png'];
                save_fl = fullfile('C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\SingleElectrodes\', eff_name, lock);
                save_loc = fullfile(save_fl,save_n);
                exportgraphics(f, save_loc, 'Resolution', 300)
                
                %% Plot Congruence effects per Block
                %limits congruency effect
                max_scale_blc = 1.2;
                min_scale_blc = -1.2;
                
                % Plot Congruence
                figure(3)
                plt = tiledlayout(2,2);
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group: MI Block')
                
                % Plot PD Group
                Group = LW.Group == 1;
                nexttile
                data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('PD - Group: MI Block')
                
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group: MC Block')
                
                % Plot PD Group
                Group = LW.Group == 1;
                nexttile
                data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('PD - Group: MC Block')
                
                colormap(flipud(brewermap([], "Rdbu")));
                h = colorbar;
                h.FontSize = 10;
                h.FontWeight = 'bold';
                t_input = 'Conflict Block Difference ' + umgsprl + ' Control: ' + ROI_elect(elect) + ' (' + lock + ')' ;
                sgt = sgtitle([ t_input]);
                sgt.FontSize = 15;
                sgt.FontWeight = 'bold'
                
                f = gcf;
                save_n = ['BlockCongruence_difference_' + ROI_elect(elect)+ '_' + eff_name + '_' + lock + '.png'];
                save_fl = fullfile('C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\SingleElectrodes\', eff_name, lock);
                save_loc = fullfile(save_fl,save_n);
                exportgraphics(f, save_loc, 'Resolution', 300)
                
                
            end
            
        end
    end
end
