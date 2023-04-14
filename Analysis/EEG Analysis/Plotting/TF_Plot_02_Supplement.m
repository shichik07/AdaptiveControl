%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Single Electrode Plotting by Condition         %
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
freq2ana         = [4, 8]; % we are interested 4 to 8 Hz

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
    
    for block_type = 1:2
        if block_type == 1
            blck = 'MC';
        else
            blck = 'MI';
        end
        clf
        fig3 = figure(3)
        for lc = 1:length(locks)
            lock = locks(lc);
            if lc == 1
                %dirs.home = 'E:\AdaptiveControl\Data\FrequencyData\';
                Time = Time_SL;
                order = [1, 2, 5, 6];
                theta_t = [0.3, 0.7];
                order = [1, 2, 5, 6];
                order2 = [9, 10, 13, 14];
            elseif lc == 2
                %dirs.home = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData';
                Time = Time_RL;
                order = [3, 4, 7, 8];
                theta_t = [-0.6, -0.2];
                order = [3, 4, 7, 8];
                order2 = [11, 12, 15, 16];
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
                    LW.Betas_intercept(:,:,:,sub)   = Betas(:,:,:, 4);%4);
                    LW.Betas_Congruency(:,:,:,sub)  = Betas(:,:,:, 1);%1);
                    LW.Betas_Block(:,:,:,sub)       = Betas(:,:,:, 2);%2);
                    LW.Betas_Interaction(:,:,:,sub) = Betas(:,:,:, 3);%3);
                    
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
            
            if blck == 'MC'
                data_In = MC_I;
                data_Con = MC_C;
            else
                data_In = MI_I;
                data_Con = MI_C;
            end
            
            elect = 1
            %electrode of interest
            ch_t = find(strcmp(channel_l, ROI_elect(elect)));
            
            %limits congruency effect
            
            max_scale_blc = -1.5;
            min_scale_blc = -4.5;
            
            %% Theta data
            % Average over 4-8Hz to make it one frequency
            frequencies2ana         = Infos.freqs_base >= freq2ana(1) & Infos.freqs_base <= freq2ana(2);
            t_ind = Time >= theta_t(1) & Time <= theta_t(2);
            Group_CO = LW.Group == 0;
            Group_PD = LW.Group == 1;
            %get data
            
            
            %CO _theta data
            theta_inc_CO = double(squeeze(mean(data_In(:, frequencies2ana, t_ind, Group_CO),4))); % average over participants
            theta_inc_CO = double(squeeze(mean(theta_inc_CO,3))); %average over time
            theta_inc_CO = double(squeeze(mean(theta_inc_CO,2))); %average over freqeuncies
            
            theta_con_CO = double(squeeze(mean(data_Con(:, frequencies2ana, t_ind, Group_CO),4))); % average over participants
            theta_con_CO = double(squeeze(mean(theta_con_CO,3))); %average over time
            theta_con_CO = double(squeeze(mean(theta_con_CO,2))); %average over freqeuncies
            theta_inc_PD = theta_inc_CO -  theta_con_CO;
            %PD _theta data
            freqdata2plot = double(squeeze(mean(data_In(:, frequencies2ana, t_ind, Group_PD),4))); % average over participants
            freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
            theta_inc_PD = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
            
            freqdata2plot = double(squeeze(mean(data_Con(:, frequencies2ana, t_ind, Group_PD),4))); % average over participants
            freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
            theta_con_PD = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
            
            % Plot the topgraphies
            
            
            plt1 = subplot(4,4,order2(1));
            topoplot(theta_inc_CO,Infos.chanlocs,'maplimits',[min_scale_blc max_scale_blc], 'electrodes', 'on', 'emarker', {'.','k',{},4});
            colormap(flipud(brewermap([], "Rdbu")));
            title('HC: Incongruent')
            
            pos = get(plt1, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt1, 'Position', posnew)
            
            
            plt2 = subplot(4,4,order2(3));
            topoplot(theta_con_CO,Infos.chanlocs,'maplimits',[min_scale_blc max_scale_blc], 'electrodes', 'on', 'emarker', {'.','k',{},4});
            colormap(flipud(brewermap([], "Rdbu")));
            title('HC: Congruent')
            
            pos = get(plt2, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt2, 'Position', posnew)
            
            
            plt3 = subplot(4,4,order2(2));
            topoplot(theta_inc_PD,Infos.chanlocs,'maplimits',[min_scale_blc max_scale_blc], 'electrodes', 'on', 'emarker', {'.','k',{},4});
            colormap(flipud(brewermap([], "Rdbu")));
            title('PD: Incongruent')
            
            pos = get(plt3, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt3, 'Position', posnew)
            
            
            plt4 = subplot(4,4,order2(4));
            topoplot(theta_con_PD,Infos.chanlocs,'maplimits',[min_scale_blc max_scale_blc], 'electrodes', 'on', 'emarker', {'.','k',{},4});
            colormap(flipud(brewermap([], "Rdbu")));
            title('PD: Congruent')
            
            pos = get(plt4, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt4, 'Position', posnew)
            
            
            
            %% Plot time frequency spectra
            
            plt5 = subplot(4,4,order(1));
            
            
            %Plot CO group
            Group = LW.Group == 0;
            data_C = squeeze(mean(data_In(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('HC: Incongruent')
            
            % slightly shift the plots
            pos = get(plt5, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt5, 'Position', posnew)
            
            
            %Plot CO group
            plt6 = subplot(4,4,order(3));
            Group = LW.Group == 0;
            data_C = squeeze(mean(data_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('HC: Congruent')
            
            % slightly shift the plots
            pos = get(plt6, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt6, 'Position', posnew)
            
            % Plot PD Group
            plt7 = subplot(4,4,order(2));
            Group = LW.Group == 1;
            
            data_C = squeeze(mean(data_In(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('PD: Incongruent')
            
            % slightly shift the plots
            pos = get(plt7, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt7, 'Position', posnew)
            
            plt8 = subplot(4,4,order(4));
            
            %Plot PD group
            Group = LW.Group == 1;
            data_C = squeeze(mean(data_Con(ch_t,:,:,Group),4));
            contourf(Time,freqs,data_C,cont_n,'linecolor','none')
            set(gca,'clim',[min_scale_blc, max_scale_blc],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
            title('PD: Congruent')
            
            % slightly shift the plots
            pos = get(plt8, 'Position')
            posnew = pos;
            posnew(1) = posnew(1) - 0.04;
            set(plt8, 'Position', posnew)
            
            
        end
        % get colorbar
        colormap(flipud(brewermap([], "Rdbu")));
        h = colorbar;
        set(h,'Position',[0.9 0.05 0.03 0.9])
        h.FontSize = 10;
        h.FontWeight = 'bold';
        t_input = 'Conflict ' + eff_name + ' ' + blck + ' Proportion: ' + ROI_elect(elect);
        sgt = suptitle([ t_input]);%sgtitle([ t_input]);
        sgt.FontSize = 15;
        sgt.FontWeight = 'bold';
        
        
        f = gcf;
        set(gcf, 'Units', 'centimeters', 'Position', [0 0 16 16]) % set to have size 16 cm width
        save_n = ['SingleConditions_' + ROI_elect(elect)+ '_' + eff_name + '_' + blck + '_' + phase_type + '.png'];
        
        save_fl = fullfile('C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\SingleElectrodes\', eff_name, lock);
        save_loc = fullfile(save_fl,save_n);
        exportgraphics(f, save_loc, 'Resolution', 300)
    end
    
end

