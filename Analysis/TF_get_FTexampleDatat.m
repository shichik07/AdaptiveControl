%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis -Plottinng %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Get fieldtrip example data                     %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error

%Linux folder locations
% eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
% % set directories
% dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
% dirs.eegsave                 = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
% dirs.functions               = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Analysis/Analysis/AdaptiveControl/Analysis';
% addpath(dirs.functions) 

%Windows folder locations
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'E:\AdaptiveControl\Data\FrequencyData\'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';
addpath(dirs.functions) 



% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

Participant_IDs          = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};  
Part_N                       = length(Participant_IDs); %number of participants

%add EEGLAB path and start the program
addpath(eegl);
addpath(ft);
addpath(spm122)
eeglab
cd('E:\AdaptiveControl\Data\FrequencyData\derivatives')
%% Export eeglab to fieldtrip

% Load one dataset to get parameters for wavelet anaylsis
sub                          = 3;
fileID                       = strcat(Participant_IDs{sub}, '_epoched_freq.set'); %get file ID
folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
fb = 'raw';
data = eeglab2fieldtrip(EEG, fb, 'none');


%% Perform Frequency Analysis

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials   ='yes'
TFRhann = ft_freqanalysis(cfg, data);


% next put the matrix into spm
fn = 'E:\AdaptiveControl\Data\FrequencyData\derivatives\abc.mat';
D = spm_eeg_ft2spm(TFRhann, fn)

dats = spm_eeg_load(fn);

display(D(channels))

% Examples of providing additional information in a script
% [] comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
D = type(D, 'single');             % Sets the dataset type
D = chantype(D, [], 'EEG');        % Sets the channel type 
D = conditions(D, [], 'Sound 1');  % Sets the condition label

% save
%--------------------------------------------------------------------------

save(D);
