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
dirs.home                    = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
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
Part_old                     = Participant_IDs(contains({Participant_IDs(:).name},{'sub-CO'}));
Part_PD                      = Participant_IDs(contains({Participant_IDs(:).name},{'sub-PD'}));
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

%% Get the data

%load one dataset to get dimensions
test                        = 17;
folderID                    = fullfile(dirs.home,Participant_IDs{test});%get folder ID
folderID_t1                 = dir(fullfile(folderID, '**', '*_npl.mat'));
folderID_t1                 = {folderID_t1(:).name}; %phase_locked files
save_loc                    = fullfile(dirs.home, Participant_IDs{test},  folderID_t1{8});
load(save_loc, 'TF_non_phase');
ds                          = size(TF_non_phase.power);

LW_MI.pl.congruent          = zeros(ds(2), ds(3), Part_N-1);
LW_MI.pl.incongruent        = zeros(ds(2), ds(3), Part_N-1);
LW_MI.npl.congruent         = zeros(ds(2), ds(3), Part_N-1);
LW_MI.npl.incongruent       = zeros(ds(2), ds(3), Part_N-1);

LW.Group = NaN(Part_N,1);

for sub = 2:Part_N
    % get file IDs
    folderID                    = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    folderID_t1                 = dir(fullfile(folderID, '**', '*_pl.mat'));
    folderID_t1                 = folderID_t1(contains({folderID_t1(:).name},{'LWPC_MI'}));
    folderID_t2                 = dir(fullfile(folderID, '**', '*_npl.mat'));
    folderID_t2                 = folderID_t2(contains({folderID_t2(:).name},{'LWPC_MI'}));
    folderID_t1                 = {folderID_t1(:).name}; %phase_locked files
    folderID_t2                 = {folderID_t2(:).name}; %non_phase locked files
    
    if contains(folderID_t1{1}, 'PD')
        LW.Group(sub-1) = 1;
    else
        LW.Group(sub-1) = 0;
    end
    
     fprintf('Load data of participant %s. \n',Participant_IDs{sub})
    
    
    for p_type = 1:2 % do this for the phase locked and non-phase locked data
        % load both the phase locked and non_phase locked data to see if it has
        % been baseline corrected
        for con = 1:2
            if p_type == 1
                fprintf('Load condition %s, of participant %s phase locked. \n',folderID_t1{con}(11:19), Participant_IDs{sub})
                save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t1{con});
                load(save_loc, 'TF_phase');
                % find FCz Channel
                channel_l = {TF_phase.chanlocs(:).labels};
                ch_t = find(strcmp(channel_l, 'FCz'));
                if  contains(folderID_t1{con}, 'LWPC_MI_C')
                    LW_MI.pl.congruent(:,:,sub) = squeeze(mean(TF_phase.power(ch_t,:,:,:), 4));
                else
                    LW_MI.pl.incongruent(:,:,sub) = squeeze(mean(TF_phase.power(ch_t,:,:,:), 4));
                    
                end
                clear  TF_phase
            elseif p_type ==2
                fprintf('Load condition %s, of participant %s non-phase locked. \n',folderID_t1{con}(11:19), Participant_IDs{sub})
                save_loc = fullfile(dirs.home, Participant_IDs{sub},  folderID_t2{con});
                load(save_loc, 'TF_non_phase');
                % find FCz Channel
                channel_l = {TF_non_phase.chanlocs(:).labels};
                ch_t = find(strcmp(channel_l, 'FCz'));
                if  contains(folderID_t1{con}, 'LWPC_MI_C')
                    LW_MI.npl.congruent(:,:,sub) = squeeze(mean(TF_non_phase.power(ch_t,:,:,:), 4));
                else
                    LW_MI.npl.incongruent(:,:,sub) = squeeze(mean(TF_non_phase.power(ch_t,:,:,:), 4));
                    
                end
                clear  TF_non_phase
            end
            
        end
    end
end
        % now check if a baseline field already exists - if so break
        
LW_MI.group = LW.Group;
LW_MI.time  = TF_phase.Time;
LW_MI.freq  = TF_phase.Frequencies;

save_name ='FCz_MI_data.mat';
Save_file = 'E:\AdaptiveControl\Data\FrequencyData\derivatives';
save_loc = fullfile(Save_file,save_name);
save(save_loc,'LW_MI','-v7.3');

%% Now let us plot first for the Parkinson Patients 

% Frist remove the first participant data, somehow that did not work as
% intended 
LW_MI.pl.congruent(:,:,1) = [];
LW_MI.pl.incongruent(:,:,1) = [];
LW_MI.npl.congruent(:,:,1) = [];
LW_MI.npl.incongruent(:,:,1) = [];

% also remove the last group index
LW_MI.group(end) = [];

% PD_group indices 
PD_ind = LW_MI.group == 1;

% First we plot congruent and incongruent for the phase locked data
max_scale = 5;%max(LW_MI.pl.congruent,[], 'all')
min_scale = 5;%min(LW_MI.pl.congruent,[], 'all')
freq_low = min(LW_MI.freq);
freq_up  = max(LW_MI.freq);
LW_MI.group == 1
data_MI_C = mean(LW_MI.pl.congruent(:, :, PD_ind), 3);

figure(1)
contourf(LW_MI.time,LW_MI.freq,data_MI_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

% next let us plot incongruent 
data_MI_I = mean(LW_MI.pl.incongruent(:, :, PD_ind), 3);

figure(2)
contourf(LW_MI.time,LW_MI.freq,data_MI_I,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

% now let us plot the Stroop effect
Stroop_MI_pl = data_MI_I - data_MI_C;

figure(3)
contourf(LW_MI.time,LW_MI.freq,Stroop_MI_pl,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

%% Now the same for the non_phase-locked data

data_MI_C_npl = mean(LW_MI.npl.congruent(:, :, PD_ind), 3);

figure(4)
contourf(LW_MI.time,LW_MI.freq,data_MI_C_npl,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

% next let us plot incongruent 
data_MI_I_npl = mean(LW_MI.npl.incongruent(:, :, PD_ind), 3);

figure(5)
contourf(LW_MI.time,LW_MI.freq,data_MI_I_npl,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

% now let us plot the Stroop effect
Stroop_MI_npl = data_MI_I_npl - data_MI_C_npl;

figure(9)
contourf(LW_MI.time,LW_MI.freq,Stroop_MI_npl,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar

        