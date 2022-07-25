%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Statistical Analysis using Fieldtrip           %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Explore beta distribution
clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
LIMO                         = 'C:\Program Files\MATLAB\LIMO';
% set directories
dirs.home                    = 'E:\AdaptiveControl\Data\FrequencyData\'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

addpath(LIMO)
addpath(eegl)
addpath(ft)
addpath(dirs.functions)



%% Load example data file for some useful information
examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_LWPC_MC_I_frequency_data_npl_bslC"));
data                        = load(examp, 'TF_non_phase');
Infos.chanlocs                    = data.TF_non_phase.chanlocs;
Infos.freqs                       = data.TF_non_phase.Frequencies;
Infos.time                        = data.TF_non_phase.Time;


% load test data set
load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, '_GLM_LWPC_pl.mat'));
load(load_loc, 'Betas');

Betas = reshape_3d_limo_fixed(Betas, Infos.freqs, Infos.time);

Beta_s = size(Betas);



%initialze data matrices
LW.Betas_intercept              = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Congruency             = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Block                  = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Interaction            = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Group                        = NaN(Part_N,1);


% Load the weight matrices and get the data
for sub = 2:Part_N
    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, '_GLM_ISPC_pl.mat'));
    load(load_loc, 'Betas');
    
    Betas                           = reshape_3d_limo_fixed(Betas, Infos.freqs, Infos.time);
    LW.Betas_intercept(:,:,:,sub)   = Betas(:,:,:, 4);
    LW.Betas_Congruency(:,:,:,sub)  = Betas(:,:,:, 1);
    LW.Betas_Block(:,:,:,sub)       = Betas(:,:,:, 2);
    LW.Betas_Interaction(:,:,:,sub) = Betas(:,:,:, 3);
    
    % get index variable
    if Participant_IDs{sub}(5:6) == 'CO'
        LW.Group(sub) = 0;
    elseif Participant_IDs{sub}(5:6) == 'PD'
        LW.Group(sub) = 1;
    end
    clear Betas
end
        
% calculate marginal means
MI_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + ( 0.25)*LW.Betas_Interaction;
MI_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + (-0.5)*LW.Betas_Block + (-0.25)*LW.Betas_Interaction;
MC_I = LW.Betas_intercept + (-0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + (-0.25)*LW.Betas_Interaction;
MC_C = LW.Betas_intercept + ( 0.5)*LW.Betas_Congruency + ( 0.5)*LW.Betas_Block + ( 0.25)*LW.Betas_Interaction;

% Caluclate Congruency Effects per block
MI_Con = MI_I - MI_C;
MC_Con = MC_I - MC_C;

% Average over 4-8Hz to make it one frequency
frequencies2ana         = [4, 8]; % we are interested 4 to 8 Hz 
frequencies2ana         = Infos.freqs >= frequencies2ana(1) & Infos.freqs <= frequencies2ana(2);
Infos.freqs(end+1)      = 1; % we add another frequency as a standin for our averaged theta data

MI_Con(:, end+1, :, :)  = squeeze(mean(MI_Con(:, frequencies2ana, :, :),2));
MC_Con(:, end+1, :, :)  = squeeze(mean(MC_Con(:, frequencies2ana, :, :),2));
% MC_C(:, end+1, :, :)  = squeeze(mean(MC_C(:, frequencies2ana, :, :),2));
% MC_I(:, end+1, :, :)  = squeeze(mean(MC_I(:, frequencies2ana, :, :),2));



% Get the data into fieldtrip format
Group = LW.Group == 0;
LWPC_MICon_CO = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
LWPC_MCCon_CO = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);
Group = LW.Group == 1;
LWPC_MICon_PD = get_pseudo_FT_struct_cluster(MI_Con(:,:,:,Group), Infos);
LWPC_MCCon_PD = get_pseudo_FT_struct_cluster(MC_Con(:,:,:,Group), Infos);

% %Test
% Group = LW.Group == 1;
% LWPC_MI_PD = get_pseudo_FT_struct_cluster(MC_I(:,:,:,Group), Infos);
% LWPC_MC_PD = get_pseudo_FT_struct_cluster(MC_C(:,:,:,Group), Infos);






%% Fieldtrip Cluster-based Permutation Test statistics

cfg = [];
cfg.channel          = {'EEG'};
cfg.latency          = [0 1];
cfg.frequency        = 1; % one the standin for our theta averaged data
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, LWPC_MICon_CO);


subj = size(LWPC_MICon_PD.powspctrm, 1);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, LWPC_MCCon_CO, LWPC_MICon_CO)


% get layout
cfg           = [];
cfg.elec      = LWPC_MICon_CO.elec;
[EasyCap128, cfg] = ft_prepare_layout(cfg)

% plot results of the analysis
cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.toi    = [0:0.1:1];
cfg.layout = EasyCap128;
ft_clusterplot(cfg, stat);

