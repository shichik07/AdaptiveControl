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

%Windows env file locations
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'E:\AdaptiveControl\Data\FrequencyData\'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.Julia                   = 'G:\Julia Ficke\AdaptiveControl\CleanData_EEG';
% %Linux env file locations
% eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
% % set directories
% dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
% dirs.eegsave                 = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden

%git directory name for me
%  /C/Users/doex9445/Dateien/Julius/AdaptiveControl/tVNS-Project

% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'derivatives'}));

Participant_IDs              = {Participant_IDs(:).name};  
Part_N                       = length(Participant_IDs); %number of participants

%add EEGLAB path and start the program
addpath(eegl);
eeglab

%Load a dataset of an old participant so we can use the full channel
%location file
sub = 4;
OldID                        = strcat(Participant_IDs{sub}, '_task-stroop_eeg_testset.set'); %get file ID
OldFile                      = fullfile(dirs.Julia,Participant_IDs{sub}, 'eeg');%get folder ID
EEG_OLD                      = pop_loadset('filename', OldID,'filepath',[OldFile]); % load file

% Because this old file doesn't have the FCz location, we add this manually
% here
EEG_OLD = pop_chanedit(EEG_OLD, 'insert',64,'changefield',{64,'labels','FCz'},'changefield',{64,'theta','0'},'changefield',{64,'radius','0.127777777777778'},'changefield',{64,'X','0.390731128489274'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.920504853452440'},'changefield',{64,'sph_theta','0'},'changefield',{64,'sph_phi','67'},'changefield',{64,'sph_radius','1'},'changefield',{64,'type','EEG'},'setref',{'1:128','FCz'});
EEG_OLD = pop_reref( EEG_OLD, [],'refloc',struct('labels',{'FCz'},'sph_radius',{1},'sph_theta',{0},'sph_phi',{67},'theta',{0},'radius',{0.12778},'X',{0.39073},'Y',{0},'Z',{0.9205},'type',{'EEG'},'ref',{'FCz'},'urchan',{[]},'datachan',{0}));


% event onset codes
onset =  {  'S 21'  'S 22'  'S 23'  'S 24'  'S 41'  'S 42'  'S 43' ...
    'S 44'  'S 61'  'S 62'  'S 63'  'S 64'  'S 81'  'S 82'  'S 83'  'S 84'  };
epoch_dur = [-1.2  2];


%% Epoching

for sub = 1:Part_N
    % get file location and load data
    fileID                       = strcat(Participant_IDs{sub}, '_task-stroop_eeg_pruned.set'); %get file ID
    folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    %Interpolate removed channels
    EEG = pop_interp(EEG, EEG_OLD.chanlocs, 'spherical');
    
    
    % load data
    EEG = pop_epoch( EEG,onset, epoch_dur, 'newname', 'Epochs for wavelet convolution', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    
    %save data
    EEG = pop_saveset(EEG, 'filename',  strcat(Participant_IDs{sub}, '_epoched_freq'),'filepath',[folderID]); 
    
    % TECHNICAL DEBT CALCULATE HOW MANY TRIALS ARE LOST
end