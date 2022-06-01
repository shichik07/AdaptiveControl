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




EEGLABFILE = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Data/EEG/Raw/sub-CO_1_epoched_freq.set';

% hdr = ft_read_header( EEGLABFILE);
% data = ft_read_data( EEGLABFILE, 'header', hdr );
% events = ft_read_event( EEGLABFILE, 'header', hdr );
% 
fileN = fullfile(dirs.eegsave, 'testset.mat');
save(fileN, 'hdr', 'data', 'events')
load(fileN)

% to convert fieldtrip data to spm I have to get my data in shape. To
% accomplish this I use a code snippet by Karl Friston, Vladimir Litvak, 
% which i adapted to my purposes 











%% Actual Conversion
% define the output file name
%--------------------------------------------------------------------------
fname = 'spm_conv.mat';
save_loc = fullfile(dirs.eegsave, fname);

% Create the Fieldtrip raw struct
%--------------------------------------------------------------------------
ftdata = [];

for i = 1:size(data,3)
   ftdata.trial{i} = data(:, :, i);
   ftdata.time{i} = -2000:1000/hdr.Fs:2000-1000/hdr.Fs;
end

ftdata.fsample = hdr.Fs; %sampling rate
ftdata.label = hdr.label(:); %electrode labels
ftdata.hdr = hdr;

% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, save_loc);
spm_eeg_ft2spm(ftdata)
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

