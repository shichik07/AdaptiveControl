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
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% Workspace on Windows PC
% /C/Users/doex9445/Dateien/Julius/AdaptiveControl/AdaptiveControl



% % Linux folder locations
% eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
% ft                           = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/fieldtrip-20220514';
% %ft2                          =
% %'/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/fieldtrip'; %This
% %version doesn't work for transfering the data
% spm122                          = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/spm12';
% % set directories
% dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
% dirs.eegsave                 = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Data/EEG/Raw/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
% dirs.functions               = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Analysis/Analysis/AdaptiveControl/Analysis';

addpath(dirs.functions)
addpath(eegl);
addpath(ft);
addpath(spm122)
eeglab


% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

Participant_IDs          = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

% Make_new directory - needed to do it for a different hard drive - only apply it once
%
% for sub = 1: Part_N
%     dir_name = fullfile(dirs.eegsave,Participant_IDs{sub});
%     mkdir(dir_name)
% end


for sub = 1:Part_N
    % get file IDs
    folderID                    = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    folderID_t1                 = dir(fullfile(folderID, '**', '*_pl.mat'));
    folderID_t2                 = dir(fullfile(folderID, '**', '*_npl.mat'));
    folderID_t1                 = {folderID_t1(:).name}; %phase_locked files
    folderID_t2                 = {folderID_t2(:).name}; %non_phase locked files
    
    cd(fullfile(dirs.home, Participant_IDs{sub}))
    
    fprintf('Perform spm conversion for participant %s. \n',Participant_IDs{sub})
    for p_type = 1:2 % do this for the phase locked and non-phase locked data
        
        if p_type == 1
            files = folderID_t1;
            file_type = 'TF_phase';
        else
            files = folderID_t2;
            file_type = 'TF_non_phase';
        end
        
        for con = 1:length(files)
            
            load_loc = fullfile(dirs.eegsave, Participant_IDs{sub},  files{con});
            data = load(load_loc);
            if p_type == 1
                data = data.TF_phase;
            else
                data = data.TF_non_phase;
            end
            
            % convert data into fieldtrip like data structure
            data = get_pseudo_FT_struct(data);
            
            % define save location - we save on the fast disk again
            new_fn = files{con};
            new_fn = strcat(new_fn(1:end-4), '_spm.mat');
            save_loc = fullfile(dirs.home, Participant_IDs{sub}, new_fn);
            subn_l = length(Participant_IDs{sub}); % get length of subject name
            condition =  new_fn(subn_l+2:subn_l+10);
            fprintf('Perform spm conversion in participant %s for condition %s. \n',Participant_IDs{sub}, condition)
            
            % convert to spm format
            D = spm_eeg_ft2spm(data, new_fn);
            
            % Examples of providing additional information in a script
            % [] comes instead of an index vector and means that the command
            % applies to all channels/all trials.
            %--------------------------------------------------------------------------
            D = type(D, 'single');             % Sets the dataset type
            D = chantype(D, [], 'EEG');        % Sets the channel type
            D = conditions(D, [], condition);  % Sets the condition label
            
            % save
            %--------------------------------------------------------------------------
            
            save(D);
            clear data D
        end
    end
end
