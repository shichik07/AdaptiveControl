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
dirs.home                    = 'G:\Julia Ficke\AdaptiveControl\CleanData_EEG'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'E:\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden

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
epoch_dur = [-1.3  2];


% INSERT PARTICIPANTS WE WISH NOT TO ANALYZE

%add EEGLAB path and start the program
addpath(eegl);
eeglab