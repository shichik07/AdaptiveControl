%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %                               
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Convert EEGlab data, 2 fieldtrip and export 2 spm %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error

eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
ft                           = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/fieldtrip-20220514';
%ft2                          =
%'/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/fieldtrip'; %This
%version doesn't work for transfering the data
spm122                          = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/spm12';
% set directories
dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Data/EEG/Raw/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Analysis/Analysis/AdaptiveControl/Analysis';
addpath(dirs.functions) 

addpath(eegl);
addpath(ft);
%addpath(ft2); % this version does not work anymore to transform the data to spm
addpath(spm122)
eeglab


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

