%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Convert EEGlab data, 2 fieldtrip and export 2 spm %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
% set directories
dirs.home                    = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData';
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';
addpath(dirs.functions)

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

addpath(eegl)
addpath(dirs.functions)

ISPC_only = 0; % when only doing the ISPC analysis

% Baseline Parameter
% Load one dataset to get parameters for wavelet anaylsis
sub                          = 3;
fileID                       = strcat(Participant_IDs{sub}, '_epoched_freq_SLplc_ResponseL.set'); %get file ID
folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file


t = -1000/1000: 1/(EEG.srate/2) : 200/1000; % divide sampling rate because we downsample data 
full_trial_bs_idx = dsearchn(t', -1000/1000): dsearchn(t', 0/1000); 
BS = "post"; %"post" select whether your want the pre or post baseline correction

    
if BS == "post"
    bsl_idx = full_trial_bs_idx;
end

      

for sub = 2:Part_N
    % get file IDs
    folderID                    = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    folderID_t1                 = dir(fullfile(folderID, '**', '*_pl_RL.mat'));
    folderID_t2                 = dir(fullfile(folderID, '**', '*_npl_RL.mat'));
    folderID_t1                 = folderID_t1(~contains({folderID_t1(:).name},{'GLM'})); % not sure where this mistake came from
    folderID_t2                 = folderID_t2(~contains({folderID_t2(:).name},{'GLM'}));
    folderID_t1                 = {folderID_t1(:).name}; %phase_locked files
    folderID_t2                 = {folderID_t2(:).name}; %non_phase locked files
    ISPC_t2                     = folderID_t2(~contains(folderID_t2(:),{'LWPC'})); %ISPC t2 data only for baseline correction
    ISPC_t1                     = folderID_t1(~contains(folderID_t1(:),{'LWPC'})); %%ISPC t1 data only for baseline correction
    
    for p_type = 1%:2 % do this for the phase locked and non-phase locked data
        % load both the phase locked and non_phase locked data to see if it has
        % been baseline corrected
        
        fprintf('Perform baseline correction for participant %s. \n',Participant_IDs{sub})
        save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t1{1});
        load(save_loc, 'TF_phase');
        save_loc2 = fullfile(dirs.home, Participant_IDs{sub},  folderID_t2{1});
        load(save_loc2, 'TF_non_phase');
        
        % now check if a baseline field already exists - if so break
        
        if isfield(TF_non_phase, 'baseline') | isfield(TF_phase, 'baseline')
            fprintf('Baseline_correction was performed already in participant %s. \n',Participant_IDs{sub})
            break
        end
        
        %remove the loaded files to not clutter your workspace
        clear TF_phase TF_non_phase
        
        
        % Load all the data to get the baseline across conditions
        if p_type == 1 % for phase_locked frequency data
            data = [];
            for con = 1:length(folderID_t1)
                save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t1{con});
                clear TF_phase % remove data to save up space
                load(save_loc, 'TF_phase');
                data = cat(4, data, TF_phase.power);
            end
        else % for non_phase locked frequency data
            data = [];
            for con = 1:length(folderID_t2)
                save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t2{con});
                clear TF_non_phase % remove data to save up space
                load(save_loc, 'TF_non_phase');
                data = cat(4, data, TF_non_phase.power);
            end
        end
        
        % calculate the baseline power
        base_power = squeeze(mean(data(:,:,bsl_idx,:),4)); % mean over trials
        base_power = squeeze(mean(base_power,3)); % mean over time
        
        % make room in your workspace
        clear data
        
        % perform baseline correction for each dataset separately
        if ISPC_only == 1;
            corrected_con_t1 = ISPC_t1;
            corrected_con_t2 = ISPC_t2;
        else
            corrected_con_t1 = folderID_t1;
            corrected_con_t2 = folderID_t2;
        end
        
            
        
        % Load all the data to get the baseline across conditions
        if p_type == 1 % for phase_locked frequency data
            data = [];
            for con = 1:length(corrected_con_t1)
                fprintf('Perform Baselin Correction for participant %s and condition %s. \n',Participant_IDs{sub}, corrected_con_t1{con}(11:19))
                load_loc = fullfile(dirs.home, Participant_IDs{sub},  corrected_con_t1{con});
                clear TF_phase % remove data to save up space
                load(load_loc, 'TF_phase');
                
                for freq = 1:size(base_power,2)
                    for chann = 1:size(base_power,1)
                        TF_phase.power(chann,freq,:,:) = 10*log10(TF_phase.power(chann,freq,:,:)./base_power(chann, freq));
                    end
                end
                
                % add indicator that data is baseline corrected
                TF_phase.baseline = 1;
                %save data on different hard drive
                new_save = strcat(corrected_con_t1{con}(1:end-4),'_bslC_RL.mat');
                save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, new_save);
                %check if directory exists
                if ~exist(fullfile(dirs.eegsave, Participant_IDs{sub})) mkdir(fullfile(dirs.eegsave, Participant_IDs{sub})); end
                
                save(save_loc,'TF_phase','-v7.3');
            end
        else % for non_phase locked frequency data
            data = [];
            for con = 1:length(corrected_con_t2)
                fprintf('Perform Baselin Correction for participant %s and condition %s. \n',Participant_IDs{sub}, corrected_con_t2{con}(11:19))
                load_loc = fullfile(dirs.eegsave, Participant_IDs{sub},  corrected_con_t2{con});
                clear TF_non_phase % remove data to save up space
                load(load_loc, 'TF_non_phase');
                
                for freq = 1:size(base_power,2)
                    for chann = 1:size(base_power,1)
                        TF_non_phase.power(chann,freq,:,:) = 10*log10(TF_non_phase.power(chann,freq,:,:)./base_power(chann, freq));
                    end
                end
                
                % add indicator that data is baseline corrected
                TF_non_phase.baseline = 1;
                %save data on different hard drive
                new_save = strcat(corrected_con_t2{con}(1:end-4),'_bslC_RL.mat');
                save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, new_save);
                %check if directory exists
                if ~exist(fullfile(dirs.eegsave, Participant_IDs{sub})) mkdir(fullfile(dirs.eegsave, Participant_IDs{sub})); end
                
                save(save_loc,'TF_non_phase','-v7.3'); 
            end
        end  
    end
end