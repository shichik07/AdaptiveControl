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

%% Variables for permutation
n_perm = 1000; % number of permutations
cluster_pval = 0.05;
voxel_pval = 0.05;

%% Load example data file for some useful information
examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_npl_bslC"));
data                        = load(examp, 'TF_non_phase');
Infos.chanlocs                    = data.TF_non_phase.chanlocs;
channel_l = {Infos.chanlocs(:).labels};
Infos.freqs_base                       = data.TF_non_phase.Frequencies;
freqs = Infos.freqs_base        ;
Infos.time                        = data.TF_non_phase.Time;
Time_SL = Infos.time;

%% load RL Time
examp                       = fullfile(dirs.eegsave, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_npl_RL_bslC_RL.mat"));
data                        = load(examp, 'TF_non_phase');
Infos.time                        = data.TF_non_phase.Time;
Time_RL = Infos.time;
%participants with too little data for a GLM
LWPC_missing = ["sub-PD_11", "sub-PD_21"];
ISPC_missing = ["sub-CO_6", "sub-PD_15", "sub-PD_21"];

ROI_elect = ["FCz", "Cz", "CPz"];
freq_low = 2; %2Hz
freq_up = 30; %30Hz
max_scale = 0.8;
min_scale = -0.8;

locks = ["SL", "RL"];

for lc = 1:length(locks)
    lock = locks(lc);
    if lc == 1
        dirs.home = 'E:\AdaptiveControl\Data\FrequencyData\';
        Time = Time_SL;
    elseif lc == 2
        dirs.home = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData';
        Time = Time_RL;
    end
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
            
            % try it for one electrode
            for grp = 1:2
                % we want to permute for each time/frequency point between
                % our two conditions to calculate the interaction. For test
                % purposes we start with the Stroop Effect in the MC
                % condition. First we pick an electrode (Cz makes sense)
                ch_t = find(strcmp(channel_l, "Cz"));
                Group = LW.Group == grp - 1;
                temp_power_C1 = squeeze(MC_Con(ch_t,:,:,Group));
                temp_power_C2 = squeeze(MI_Con(ch_t,:,:,Group));
                
                %create null hypotheses matrices
                permuted_maxvals = zeros(n_perm,2,length(freqs));
                permuted_vals    = zeros(n_perm,length(freqs),numel(Time));
                max_clust_info   = zeros(n_perm,1);
                
                % get the real t_values
                tnum   = squeeze(mean(temp_power_C2,3)-mean(temp_power_C1,3));
                tdenom = sqrt( (std(temp_power_C2,0,3).^2)./sum(Group) + (std(temp_power_C1,0,3).^2)./sum(Group));
                realt   = tnum./tdenom;
                %tmap(-2 < tmap & tmap < 2) = 0;
                
                % we will perform a paired permutation test - flip the
                % label randomly for each participant
                
                for i = 1:n_perm %perform the analysis for each group seprately
                    switches = sign(randn(sum(Group),1)); % participants for which values will be switched
                    temp_store = temp_power_C2; % store a temporary array
                    % replace values fo C2
                    temp_power_C2(:,:,(switches > 0)) = temp_power_C1(:,:,(switches > 0));
                    % now do the same C1
                    temp_power_C1(:,:,(switches > 0)) = temp_store(:,:,(switches > 0));
                    
                    % compute t-map of null hypothesis
                    tnum   = squeeze(mean(temp_power_C2,3)-mean(temp_power_C1,3));
                    tdenom = sqrt( (std(temp_power_C2,0,3).^2)./sum(Group) + (std(temp_power_C1,0,3).^2)./sum(Group));
                    tmap   = tnum./tdenom;
                    
                    % save permuted values
                    permuted_vals(i,:,:) = tmap;
                    
                    % now for the clustercorrection part I am using some of
                    % cohens script, which to be honest, I do not quite
                    % understand
                    
                    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
                    % note that here, clusters were obtained by parametrically thresholding
                    % the t-maps
                    tmap(abs(tmap)<tinv(1-voxel_pval,sum(Group)-1))=0;
                    
                    % get number of elements in largest supra-threshold cluster
                    clustinfo = bwconncomp(tmap);
                    max_clust_info(i) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
                end
                % now compute Z-map
                zmap = (realt-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));
                %limits congruency effect
                max_scale_con = 3;
                min_scale_con = -3;
                
                % Plot Congruence
                figure(1)
                plt = tiledlayout(1,2);
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = tmap;%squeeze(mean(Congruency_effect(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_con, max_scale_con],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group')

                
            end
        end
            
            ch_t = find(strcmp(channel_l, "Cz"));
                
                %% Cognruency Effect
                
                %limits congruency effect
                max_scale_con = 3;
                min_scale_con = -3;
                
                % Plot Congruence
                figure(1)
                plt = tiledlayout(1,2);
                nexttile
                
                %Plot CO group
                Group = LW.Group == 0;
                data_C = tmap;%squeeze(mean(Congruency_effect(ch_t,:,:,Group),4));
                contourf(Time,freqs,data_C,40,'linecolor','none')
                set(gca,'clim',[min_scale_con, max_scale_con],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
                title('HC - Group')
                
            
        end
    end
end