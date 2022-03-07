%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %                               
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Epoching for Frequency Analysis                %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Short Script using EEGlab to epoch the automatically preprocessed data to
%segments of length 3 seconds pre-onset to 3- seconds post onset. 

%load eeg data

clear all; close all; clc
dbstop if error

eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
% set directories
dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden

%git directory name for me
%  /C/Users/doex9445/Dateien/Julius/AdaptiveControl/tVNS-Project

% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

%remove PD 16, index 68 - no dat available
Participant_IDs(68)          = []; 
Participant_IDs              = {Participant_IDs(:).name};  
Part_N                       = length(Participant_IDs); %number of participants

% event onset codes
onset =  {  'S 21'  'S 22'  'S 23'  'S 24'  'S 41'  'S 42'  'S 43' ...
    'S 44'  'S 61'  'S 62'  'S 63'  'S 64'  'S 81'  'S 82'  'S 83'  'S 84'  };
epoch_dur = [-2  2];


% INSERT PARTICIPANTS WE WISH NOT TO ANALYZE

%add EEGLAB path and start the program
addpath(eegl);
eeglab

%% Epoching

for sub = 1:Part_N
    % get file location and load data
    fileID                       = strcat(Participant_IDs{sub}, '_task-stroop_eeg_pruned_auto.set'); %get file ID
    folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    % load data
    EEG = pop_epoch( EEG,onset, epoch_dur, 'newname', 'Epochs for wavelet convolution', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    
    %save data
    EEG = pop_saveset(EEG, 'filename',  strcat(Participant_IDs{sub}, '_epoched_freq'),'filepath',[folderID]); 
    
    % TECHNICAL DEBT CALCULATE HOW MANY TRIALS ARE LOST
end