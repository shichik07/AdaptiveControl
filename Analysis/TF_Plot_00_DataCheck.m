%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis -Plottinng %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Explore time frequency analysis results        %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After baseline correction we will have a look at the data we plot the
% list-wise condition congruent and incongruent and the interaction between
% both, in other words the difference. We do that at the FCz electrode both
% for the old and PD participants respectively

clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';
addpath(dirs.functions)

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
%addpath(ft);
%addpath(ft2); % this version does not work anymore to transform the data to spm


% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

%% Get the data


for sub = 1:Part_N
    % get file IDs
    folderID                    = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    folderID_t1                 = dir(fullfile(folderID, '**', '*_pl.mat'));
    folderID_t2                 = dir(fullfile(folderID, '**', '*_npl.mat'));
    folderID_t1                 = {folderID_t1(:).name}; %phase_locked files
    folderID_t2                 = {folderID_t2(:).name}; %non_phase locked files
    
    
    for p_type = 1:2 % do this for the phase locked and non-phase locked data
        % load both the phase locked and non_phase locked data to see if it has
        % been baseline corrected
        
        fprintf('Perform baseline correction for participant %s. \n',Participant_IDs{sub})
        save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t1{1});
        load(save_loc, 'TF_phase');
        save_loc2 = fullfile(dirs.home, Participant_IDs{sub},  folderID_t2{1});
        load(save_loc2, 'TF_non_phase');
        
        % now check if a baseline field already exists - if so break
        