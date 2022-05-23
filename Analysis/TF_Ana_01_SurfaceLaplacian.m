
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %                               
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Surface Laplacian Filter                       %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error

eegl                       = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
CSD_toolb                    = 'C:\Program Files\MATLAB\CSDtoolbox'; %toolbox for surface laplacian calculations

% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'E:\AdaptiveControl\Data\FrequencyData\'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden

% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs              = {Participant_IDs(:).name};  
Part_N                       = length(Participant_IDs); %number of participants

addpath(eegl)
addpath(genpath(CSD_toolb))
cd(dirs.home)

eeglab

%% Laplacian for EEG Data (https://sccn.ucsd.edu/pipermail/eeglablist/2020/015739.html)
for sub = 1:Part_N
    
    % get file location and load data
    fileID                       = strcat(Participant_IDs{sub}, '_epoched_freq.set'); %get file ID
    folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    %save names for converted channel locations
    from_Chan                    = fullfile(folderID, 'chanlocs.ced');
    to_Chan                      = fullfile(folderID, 'chanlocs.csd');
 
    %convert EEGlab location file from .ced or .locs to .csd
    EEG                          = pop_chanedit(EEG, 'save', from_Chan);
    ConvertLocations(from_Chan, to_Chan);
    
    % design table with channels labels
    E                            = cell(length(EEG.chanlocs),1);
    E(:,1)                       = {EEG.chanlocs(:).labels};
    
    %Extraction EEG montage to an array structure M using the CSD toolbox function
    M                            = ExtractMontage(to_Chan,E);
    
    %Generate Transformation Matrices G and H
    [G,H]                        = GetGH(M);
    
    % use single data precision
    data                         = single(repmat(NaN,size(EEG.data)));
    
    fprintf('Performing Surface-Laplacian for participant %s. \n',Participant_IDs{sub})
    % applying a CSD transform to a 3-D data matrix used in EEGlab
    for ne = 1:length(EEG.epoch)               % loop through all epochs
        myEEG                    = single(EEG.data(:,:,ne));      % reduce data precision to reduce memory demand
        MyResults                = CSD(myEEG,G,H);            % compute CSD for <channels-by-samples> 2-D epoch
        data(:,:,ne)             = MyResults;              % assign data output
    end
    looping_CSD_final            = double(data);          % final CSD data
    EEG.data                     =looping_CSD_final;
    
    % save data
    EEG = pop_saveset(EEG, 'filename',  strcat(Participant_IDs{sub}, '_epoched_freq_SLplc'),'filepath',[folderID]); 
    
    
end