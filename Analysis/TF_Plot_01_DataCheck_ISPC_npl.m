%% Plot the beta values for the ISPC effect for the non-phase locked activity
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
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_9'})); 
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'sub-PD_16'})); %remove PD 16, no data available and was excluded
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

examp                       = fullfile(dirs.home, Participant_IDs{4}, strcat(Participant_IDs{4} ,"_ISPC_MC_I_frequency_data_npl_bslC"));
data                        = load(examp);
can_locs                    = data.TF_non_phase.chanlocs;
freqs                       = data.TF_non_phase.Frequencies;
Time                        = data.TF_non_phase.Time;


% load test data set
load_loc = fullfile(dirs.home, Participant_IDs{2}, 'LIMO_Output', strcat(Participant_IDs{2}, '_GLM_LWPC_pl.mat'));
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
    load_loc = fullfile(dirs.home, Participant_IDs{sub}, 'LIMO_Output', strcat(Participant_IDs{sub}, '_GLM_LWPC_npl.mat'));
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


% find FCz Channel and Group
channel_l = {can_locs(:).labels};
ch_t = find(strcmp(channel_l, 'Cz'));
Group = LW.Group == 0;
freqs2plt    = [4, 8];
band_info    = 'Theta';
effect_info    = 'LWPC(npl)';

max_scale = 0.8;%max(LW_MI.pl.congruent,[], 'all')
min_scale = -0.8;%min(LW_MI.pl.congruent,[], 'all')
max_int = 0.5;%max(LW_MI.pl.congruent,[], 'all')
min_int = -0.5;%min(LW_MI.pl.congruent,[], 'all')
freq_low = min(freqs);
freq_up  = max(freqs);
%% Healthy participants Cz Electrode
ch_t = find(strcmp(channel_l, 'Cz'));
Group = LW.Group == 0;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(1)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MI (FCz) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(2)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MC (FCz) LWPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(3)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_int, max_int],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Interaction (Cz) LWPC')
colormap(jet)
colorbar

%% Parkinson participants Cz Electrode
ch_t = find(strcmp(channel_l, 'Cz'));
Group = LW.Group == 1;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(4)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MI (Cz) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(6)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MC (Cz) LWPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(6)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Interaction (Cz) LWPC')
colormap(jet)
colorbar

%% Healthy participants Pz Electrode
ch_t = find(strcmp(channel_l, 'Cz'));
Group = LW.Group == 0;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(7)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MI (Pz) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(8)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MC (Pz) ISPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(9)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Interaction (Pz) ISPC')
colormap(jet)
colorbar

%% Parkinson participants Pz Electrode
ch_t = find(strcmp(channel_l, 'Cz'));
Group = LW.Group == 1;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(10)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MI (Pz) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(11)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MC (Pz) ISPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(12)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_int, max_int],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Interaction (Pz) ISPC')
colormap(jet)
colorbar

%% Healthy participants FCZ Electrode
ch_t = find(strcmp(channel_l, 'FCz'));
Group = LW.Group == 0;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(7)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MI (FCZ) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(8)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Congruence MC (FCZ) ISPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(9)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Healthy Interaction (FCZ) ISPC')
colormap(jet)
colorbar

%% Parkinson participants FCZ Electrode
ch_t = find(strcmp(channel_l, 'FCz'));
Group = LW.Group == 1;

%Define Data
data_C = squeeze(mean(MI_Con(ch_t,:,:,Group),4));
figure(10)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MI (FCZ) ISPC')
colormap(jet)
colorbar


data_C = squeeze(mean(MC_Con(ch_t,:,:,Group),4));
figure(11)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Congruence MC (FCZ) ISPC')
colormap(jet)
colorbar



data_C = squeeze(mean(M_Inter(ch_t,:,:,Group),4));
figure(12)
contourf(Time,freqs,data_C,40,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power - Parkinson Interaction (FCZ) ISPC')
colormap(jet)
colorbar


%% topoplots MI_Con Theta CO LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = MI_Con;
Group = LW.Group == 0;


figure(1)
%tiledlayout(3,5)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.8 0.8]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end

colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect MI ', effect_info, ' CO Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_CO_',effect_info ,'_', band_info,'_MI.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)



%% topoplots MC_Con Theta CO LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = MC_Con;
Group = LW.Group == 0;


figure
%tiledlayout(3,5)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.8 0.8]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end

colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect MC ', effect_info, ' CO Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_CO_',effect_info ,'_', band_info,'_MC.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)

%% topoplots interaction Theta CO LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = M_Inter;
Group = LW.Group == 0;


figure(1)
%tiledlayout(6,6)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.5 0.5]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end
colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect Interaction ', effect_info, ' CO Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_CO_',effect_info ,'_', band_info,'_Interaction.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)


%% topoplots MC_Con Theta PD LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = MC_Con;
Group = LW.Group == 1;


figure
%tiledlayout(3,5)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.8 0.8]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end

colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect MC ', effect_info, ' PD Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_PD_',effect_info ,'_', band_info,'_MC.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)

%% topoplots MI_Con Theta PD LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = MI_Con;
Group = LW.Group == 1;


figure
%tiledlayout(3,5)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.8 0.8]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end

colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect MI ', effect_info, ' PD Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_PD_',effect_info ,'_', band_info,'_MI.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)

%% topoplots interaction Theta PD LWPC
frequencies2plot    = freqs2plt; % we plot 2 to 8 Hz data in intervals of 50ms
frequencies2plot    = freqs >= frequencies2plot(1) & freqs <= frequencies2plot(2);
times2plot          = [-0.150:0.075:0.900]; 
ch_t                = find(strcmp(channel_l, 'FCz'));
data2plot           = M_Inter;
Group = LW.Group == 1;


figure
%tiledlayout(3,5)
for i=1:length(times2plot)-1
    subplot(3,5,i)
    %nexttile; 
    
    %get index of time
    t_ind = Time >= times2plot(i) & Time <= times2plot(i+1);
    
    
    % extract EEG data and replace FC4 data with noise
    freqdata2plot = double(squeeze(mean(data2plot(:, frequencies2plot, t_ind, Group),4))); % average over participants
    freqdata2plot = double(squeeze(mean(freqdata2plot,3))); %average over time
    freqdata2plot = double(squeeze(mean(freqdata2plot,2))); %average over freqeuncies
    
    topoplot(freqdata2plot,can_locs,'maplimits',[-0.5 0.5]);
    title([ num2str(times2plot(i)*1000) ' to ' num2str(times2plot(i+1)*1000) ' ms' ])
end

colormap(flipud(brewermap([], "Rdbu")));
h = colorbar;
set(h,'Position',[0.93 0.05 0.03 0.9])
h.FontSize = 10;
h.FontWeight = 'bold';
sgt = sgtitle([band_info, ' Congruency Effect Interaction ', effect_info, ' PD Group']);
sgt.FontSize = 15;
sgt.FontWeight = 'bold'

f = gcf;

save_n = ['Topo_PD_',effect_info ,'_', band_info,'_Interaction.png'];
save_fl = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\Figures\StimulusLocked Data'
save_loc = fullfile(save_fl,save_n);
exportgraphics(f, save_loc, 'Resolution', 300)

