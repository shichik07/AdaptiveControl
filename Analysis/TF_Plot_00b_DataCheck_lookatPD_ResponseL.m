clear all; close all; clc
dbstop if error

% Windows folder locations
spm122                       = 'C:\Program Files\MATLAB\EEGSoftware\spm12';
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware\eeglab2022.0';
ft                           = 'C:\Program Files\MATLAB\Fieldtrip';
LIMO                         = 'C:\Program Files\MATLAB\LIMO';
% set directories
dirs.home                    = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_9'})); % data looks really messed up
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

addpath(LIMO)
addpath(eegl)
addpath(dirs.functions)
eeglab

% Silly me did not save the channel locations and times, so we have to load
% a dataset that contains it

examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_LWPC_MC_I_frequency_data_npl_RL_bslC_RL.mat"));
data                        = load(examp);
% can_locs                    = data.TF_phase.chanlocs;
% freqs                       = data.TF_phase.Frequencies;
% Time                        = data.TF_phase.Time;


can_locs                    = data.TF_non_phase.chanlocs;
freqs                       = data.TF_non_phase.Frequencies;
Time                        = data.TF_non_phase.Time;


% load test data set
load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, '_GLM_LWPC_npl_RL.mat'));
load(load_loc, 'Betas');

Betas = reshape_3d_limo_fixed(Betas, freqs, Time);

Beta_s = size(Betas);



%initialze data matrices
LW.Betas_intercept              = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Congruency             = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Block                  = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Betas_Interaction            = NaN(Beta_s(1),Beta_s(2),Beta_s(3), Part_N);
LW.Group                        = NaN(Part_N,1);


% Load the weight matrices and get the data
for sub = 2:Part_N
    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, '_GLM_LWPC_npl_RL.mat'));
    load(load_loc, 'Betas');
    
    Betas                           = reshape_3d_limo_fixed(Betas, freqs, Time);
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

% Interaction 
M_Inter = MC_Con - MI_Con;

% get baseline
Baseline = LW.Betas_intercept;

% find FCz Channel and Group
channel_l = {can_locs(:).labels};
ch_t = find(strcmp(channel_l, 'Fz'));
Group = LW.Group == 0;
freqs2plt    = [3, 7];

max_scale = 2;%max(LW_MI.pl.congruent,[], 'all')
min_scale = -6;%min(LW_MI.pl.congruent,[], 'all')
freq_low = min(freqs);
freq_up  = max(freqs);
%% find participants with messy data by checking the baseline
ch_t = find(strcmp(channel_l, 'FCz'));
Group = LW.Group == 1;

%Define Data
data_total = squeeze(MC_I(ch_t,:,:,Group));
for part_n = 1:sum(Group)
    data_C = squeeze(data_total(:, :, part_n));
    figure(part_n)
    contourf(Time,freqs,data_C,40,'linecolor','none')
    set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
    title(['Power - PD participant ', num2str(part_n) ,' Congruence MI (FCz) ISPC'])
    colormap(jet);
    colorbar
end

%after looking through the data, participant 29 in the PD group looks
%really suspicious