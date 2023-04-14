%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Calculate Trial Nrs           %                               
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
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
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

Participant_IDs          = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs          = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
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
response =  {  'S121'  'S122'  'S123'  'S124'  'S141'  'S142'  'S143' ...
    'S144'  'S161'  'S162'  'S163'  'S164'  'S181'  'S182'  'S183'  'S184'  };
epoch_dur = [-2  1.2];

% create table to save
sz = [Part_N 9];
varTypes = ["string","double","double", "double","double", "double","double", "double","double"];
varNames = ["Part_ID", "LWPC_incon_I", "LWPC_incon_C", "LWPC_con_I", "LWPC_con_C", "ISPC_incon_I", "ISPC_incon_C", "ISPC_con_I", "ISPC_con_C"];
Trial_info = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames)


%% Epoching

for sub = 1:Part_N
    % get file location and load data
    fileID                       = strcat(Participant_IDs{sub}, '_task-stroop_eeg_pruned.set'); %get file ID
    folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    savefolderID                     = fullfile(dirs.eegsave,Participant_IDs{sub});%get folder ID
    EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    % epoch data
    EEG = pop_epoch( EEG, response, epoch_dur, 'newname', 'Epochs for wavelet convolution', 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    
    try 
        trl_types = get_trlindices(EEG);
    
        % update table
        Trial_info(sub,1) = {Participant_IDs{sub}};
        Trial_info(sub,2) = {length(trl_types.LWPC_MI_I)};
        Trial_info(sub,3) = {length(trl_types.LWPC_MI_C)};
        Trial_info(sub,4) = {length(trl_types.LWPC_MC_I)};
        Trial_info(sub,5) = {length(trl_types.LWPC_MC_C)};
        Trial_info(sub,6) = {length(trl_types.ISPC_MI_I)};
        Trial_info(sub,7) = {length(trl_types.ISPC_MI_C)};
        Trial_info(sub,8) = {length(trl_types.ISPC_MC_I)};
        Trial_info(sub,9) = {length(trl_types.ISPC_MC_C)};
    catch
        warning(['Participant ', Participant_IDs{sub} , ' did not have a full dataset'])
    end
end

% We identified the datasets of the following participants on the basis
% that they had fewer than 40 trials in one category left: SubCO6 had fewer
% than 40 trials in the ISPC effect. In the PD group participant 21 has to
% be excluded entirely, participant PD_11 had fewer than 40 trials for the
% LWPC effect and participant 15 had fewer than 40 in the ISPC effect