%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %                               
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% FreqAnalysis                                   %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assumptions about the FT
% Data are stationary - there is no change in frequencies over time

% Convolution theorem - Convolution in the time domain is the same as
% multiplication in the frequency domain. You can either flip the kernel
% and slide it along the time axis or - you can compute the fft and the
% kernel - multiply them and then take the inverse fft and calculate the
% signal back. WHy? Frequency convolution is fast where it is slow in the
% time domain. When you perform frequency by freqeuncy multiplication you
% are scaling the freqeuncy spectrum. This is also the reason why
% convolution can be understood as a freqeuncy domain filter.

% First sine wave is a sine wave of zero frequency also termed DC offset -
% shows offset of the overall signal. But which frequencies did we actually
% compute? The number of frequencies that we compute is exactly one half of
% the number of datapoints in the time series plus the zero frequency
% (Nyquist Theorem). So in our case 10/5 - the highest frequency is 5Hz.
% Power is the length of the vector squared.

%load eeg data

% Workspace on Windows PC
% /C/Users/doex9445/Dateien/Julius/AdaptiveControl/AdaptiveControl 


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

% INSERT PARTICIPANTS WE WISH NOT TO ANALYZE


% Trial indices in question

%add EEGLAB path and start the program
addpath(eegl);
eeglab

%% Wavelets

% Load one dataset to get parameters for wavelet anaylsis
sub                          = 3;
fileID                       = strcat(Participant_IDs{sub}, '_epoched_freq.set'); %get file ID
folderID                     = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file

% resample to reduce the size of our dataset
%EEG = pop_resample( EEG, 125);

% cd('/home/jules/Dropbox/PhD_Thesis/Studying/NeuralTimesSeries')
% load sampleEEGdata
% 

% Define params wavelet
length_wavelet               = 4; % wavelet lengths will be 4 seconds
wave_pnts                    = length_wavelet* EEG.srate;
cycle_num                    = 6; % number of wavelet cycles
freq_up                      = 30; % upper frequency limit in Hz
freq_low                     = 2; %lower freqeuncy limit in Hz
freq_num                     = 20; %number of freqeuncies to be estimated
freq_range                   = logspace(log10(freq_low), log10(freq_up), freq_num); %range of frequencies
limit                        = 262144; % limit of what matrix operations my PC is capable of performing - important for zero padding
bin_nr                       = 10; %number of bins in which we analyze the data

% time vector for wavlet (-2:2s)
time                         = -wave_pnts/EEG.srate/2: 1/EEG.srate : wave_pnts/EEG.srate/2;
%vector of standard deviations per frequency
%s                            = cycle_num./(2*pi*freq_range); % single cycle wavelet
s    = logspace(log10(3),log10(10),freq_num )./(2*pi*freq_range); % logarithmically spaced wavelets - more precision for higher frequency bands, higher temporal precision for lower frequency bands


% create wavelets
wavelets                     = zeros(freq_num, wave_pnts + 1);
for fr = 1:freq_num
    scaling_factor = sqrt(1/(s(fr)*sqrt(pi)));
    c_sine = exp(1i*2*pi*freq_range(fr).*time);
    gaus = exp(-(time.^2)./(2*(s(fr)^2)));
    wavelets(fr,:) = scaling_factor * c_sine.* gaus;
end


% Baseline indices prior to fixation presentation in samples
fix_baseline = [-0.2 -0.1]*EEG.srate;
Baseline_time = dsearchn(EEG.times',[-300 -100]'); 
keep_time =  dsearchn(EEG.times',[-200 1000]'); % time range we want to keep in the final dataset
New_trial_time = -200/1000: 1/(EEG.srate/2) : 1000/1000; % divide sampling rate because we downsample data
whole_trl_bsl =  dsearchn(EEG.times',[-200 1000]'); % we use the whole epoch before performing the baseline correction



% briefly inspect wavelets
 figure (5)
 plot(time, real(wavelets(1,:)))
%     

%% Perform Wavelet analysis on participant data

for sub = 2%:Part_N
    %load data set
    fileID              = strcat(Participant_IDs{sub}, '_epoched_freq.set'); %get file ID
    folderID            = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    EEG                 = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    % resample to reduce the size of our dataset
    %EEG = pop_resample( EEG, 125);
    
    
    
    
    Items = get_trlindices(EEG); %get indices of trials by condition
    Fnames = fieldnames(Items); % get associated condition names
    
    % create struct where we save our averaged frequency condition data
    %     vals = zeros(EEG.nbchan, freq_num , length(keep_time(1):keep_time(2)), EEG.trials);
    %     vals = repmat({vals},1,length(Fnames));
    %     args=[Fnames';vals];
    %     vals2 = zeros(EEG.nbchan, freq_num , length(keep_time(1):keep_time(2)));
    %     vals2 = repmat({vals2},1,length(Fnames));
    %     args3 = [Fnames';vals2];
    %     TF_phase.power = struct(args{:}); % create a struct with matrix for each category for power values
    %     TF_phase.itpc  = struct(args3{:}); % create a struct with matrix for each category for inter trial phase clustering values
    % %     TF_non_phase.power = struct(args{:}); % create a struct with matrix for each category for power values
    % %     TF_non_phase.itpc  = struct(args3{:}); % create a struct with matrix for each category for inter trial phase clustering values
    %     args2 = [Fnames';num2cell(structfun(@numel,Items))']; %get numer of items in case of weighting
    %     TF_phase.Item_nr = struct(args2{:});
    % %     TF_non_phase.Item_nr = struct(args2{:});
    %     TF_phase.chanlocs = EEG.chanlocs; % keep channel information
    %     TF_phase.Frequencies = freq_range; % keep frequency information that are modelled
    %     TF_phase.Time = New_trial_time; % keep time indices that we are tying to model
    % %     TF_non_phase.chanlocs = EEG.chanlocs; % keep channel information
    % %     TF_non_phase.Frequencies = freq_range; % keep frequency information that are modelled
    % %     TF_non_phase.Time = New_trial_time; % keep time indices that we are tying to model
    %
    %
    %     clear vals args args2 % clear these to free up workspace
    %
    fprintf('Performing wavelet convolution on participant %s. \n',Participant_IDs{sub})
    
    for con = 1:length(Fnames)
        item_nr = structfun(@numel,Items); % get nr of items we have to save
        item_nr = item_nr(con);
        vals = zeros(EEG.nbchan, freq_num , length(New_trial_time),item_nr); % prepare datamatrix
        vals2 = zeros(EEG.nbchan, freq_num , length(New_trial_time));
        
        
        trl_indices = getfield(Items, Fnames{con});
        fprintf('Analyzing trial subset %s of participant %s. \n',Fnames{con}, Participant_IDs{sub})
        % get bin data and parameter
        con_data = EEG.data(:,:,trl_indices);
        con_trl = size(con_data,3); % nr of trials for this category
        % get condition specific ERP
        con_ERP = mean(con_data,3);
        
        % define convolution parameters (used from example code by Mike C. Cohen)
        n_wavelet            = length(time);
        n_data               = EEG.pnts*con_trl;
        n_convolution        = n_wavelet+n_data-1;
        n_conv_pow2          = pow2(nextpow2(n_convolution));
        half_of_wavelet_size = (n_wavelet-1)/2;
        
        if n_conv_pow2 > limit
            % double check that PC can actually exectue concolution on data
            % this size. Otherwise break program instead of PC
            break
        end
        
        % perform convolution on first data bin
        for p_type = 1:2 %analyze both phaselocked and non_phaselocked data
            if p_type == 2
                fprintf('Analyzing trial non-phase-locked power in subset %s of participant %s. \n',Fnames{con}, Participant_IDs{sub})
                % create data structure to save results
                TF_non_phase.power = vals; %save power
                TF_non_phase.itpc = vals2; %save itpc
                TF_non_phase.Item_nr = item_nr; %save item_nr
                TF_non_phase.chanlocs = EEG.chanlocs; % keep channel information
                TF_non_phase.Frequencies = freq_range; % keep frequency information that are modelled
                TF_non_phase.Time = New_trial_time; % keep time indices that we are tying to model
            else
                fprintf('Analyzing trial phase-locked power in subset %s of participant %s. \n',Fnames{con}, Participant_IDs{sub})
                % create data structure to save results
                TF_phase.power = vals; %save power
                TF_phase.itpc = vals2; %save itpc
                TF_phase.Item_nr = item_nr; %save item_nr
                TF_phase.chanlocs = EEG.chanlocs; % keep channel information
                TF_phase.Frequencies = freq_range; % keep frequency information that are modelled
                TF_phase.Time = New_trial_time; % keep time indices that we are tying to model
            end
            for chan = 1: EEG.nbchan
                for freq = 1:freq_num
                    if p_type == 2
                        dat = squeeze(con_data(chan,:,:)) - repmat(con_ERP(chan,:)', 1, size(con_data,3)); %subtract condition sepcific ERP from each trial
                        dat = reshape(dat, [1, EEG.pnts*con_trl]); %reshape trial into vector
                    else
                        % reshape trials into one vector
                        dat = reshape(squeeze(con_data(chan,:,:)), [1, EEG.pnts*con_trl]);
                    end
                    % get fft of data
                    fft_dat = fft(dat, n_conv_pow2);
                    
                    fft_wavelet = fft(wavelets(freq,:),n_conv_pow2);
                    
                    % convolve data and wavelet
                    decomp = fft_dat.*fft_wavelet;
                    
                    % transform back
                    decomp = ifft(decomp);
                    decomp = decomp(1:n_convolution);
                    decomp = decomp(half_of_wavelet_size + 1:end - half_of_wavelet_size);
                    % get original trial shape back
                    decomp = reshape(decomp, [EEG.pnts,con_trl]);
                    
                    % extract itpc
                    %itpc(chan, freq, :) = abs(mean(exp(1i*angle(decomp)),2));
                    itpc = abs(mean(exp(1i*angle(decomp)),2));
                    
                    % take mean of all trials, compute magnitude and square to
                    % get power
                    decomp = abs(decomp).^2;
                    
                    
                    % find indices per trial that need to be adjusted
                    fix_time = zeros(con_trl, 2);
                    base_power = zeros(con_trl,1);
                    baseidx = [];
                    
                    for trl = 1:con_trl
                        trl_idx = trl_indices(trl);
                        event_lat = cell2mat(EEG.epoch(trl_idx).eventlatency);
                        baseidx = event_lat(dsearchn(event_lat',0)-1); % get the correct time stamp
                        baseidx = dsearchn(EEG.times',baseidx);% translate to indices
                        fix_time(trl,:) = baseidx + Baseline_time; %safe indices and subtract samples for baseline indices
                        %base_power(trl) = mean(decomp(fix_time(trl,1): fix_time(trl,2), trl)); %baseline blank screen period
                        base_power(trl) = mean(decomp(whole_trl_bsl(1): whole_trl_bsl(2), trl)); % baseline over trial period
                        %base_power(trl) = mean(decomp(Baseline_time(1): Baseline_time(2), trl)); %baseline during fixation
                        %base_power(trl)  = bp(chan);
                    end
                    
                    %                     temppower = zeros(EEG.pnts,con_trl);
                    %                     % perform baseline correction and convert to decible scale
                    %                     for trl = 1:con_trl
                    %                         temppower(:,trl) = 10*log10(decomp(:,trl)./base_power(trl));
                    %                     end
                    
                    %                     if p_type == 2
                    %                         %save power in power-matrix
                    %                         TF.power_non_phase.(Fnames{con})(chan,freq,:) = mean(temppower(keep_time(1):keep_time(2),:),2);
                    %                         %save itpc in itpc-matrix
                    %                         TF.itpc_non_phase.(Fnames{con})(chan,freq,:) = itpc(keep_time(1):keep_time(2)); %only keep period we are interested in
                    %                     else
                    %                         %save power in power-matrix
                    %                         TF.power_phase.(Fnames{con})(chan,freq,:) = mean(temppower(keep_time(1):keep_time(2),:),2);
                    %                         %save itpc in itpc-matrix
                    %                         TF.itpc_phase.(Fnames{con})(chan,freq,:) = itpc(keep_time(1):keep_time(2)); %only keep period we are interested in
                    %                     end
                    
                    % after the analysis we downsample to 125 Hz to have a smaller sample size but without loosing information
                    pow2keep = downsample(decomp(keep_time(1):keep_time(2),:),2);
                    itpc2keep = downsample(itpc(keep_time(1):keep_time(2)), 2);
                    
                    if p_type == 2
                        TF_non_phase.power(chan,freq,:, trl_indices) =  pow2keep;
                        %save itpc in itpc-matrix
                        TF_non_phase.itpc(chan,freq,:) = itpc2keep; %only keep period we are interested in
                        
                    else
                        %save not baseline corrected power in power-matrix
                        TF_phase.power(chan,freq,:, trl_indices) = pow2keep;
                        %save itpc in itpc-matrix
                        TF_phase.itpc(chan,freq,:) = itpc2keep; %only keep period we are interested in
                    end
                    
                end
                clear decomp temppower fft_dat fft_wavelet dat
            end
            fprintf('Saving frequency data of participant %s. \n',Participant_IDs{sub})
            
            % Save data
            if p_type == 2
                save_name = strcat(Participant_IDs{sub}, '_', Fnames{con} ,'_frequency_data_npl.mat');
                save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, save_name);
                save(save_loc,'TF_non_phase','-v7.3');
                clear TF_non_phase % remove data to save up space
            else
                save_name = strcat(Participant_IDs{sub}, '_', Fnames{con} ,'_frequency_data_pl.mat');
                save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, save_name);
                save(save_loc,'TF_phase','-v7.3');
                clear TF_phase % remove data to save up space
            end
            fprintf('Frequency data of participant %s has been saved. \n',Participant_IDs{sub})
        end
    end
    
    % save power and itpc data, baseline correction will be performed
    % afterwards
    %     fprintf('Saving frequency data of participant %s. \n',Participant_IDs{sub})
    %
    %     save_name = strcat(Participant_IDs{sub}, '_frequency_data_phaselocked.mat');
    %     save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, save_name);
    %     save(save_loc,'TF_phase','-v7.3');
    %
    %     fprintf('Frequency data of participant %s has been saved. \n',Participant_IDs{sub})
    %
    %remove intermediate variables to make some space
    clear EEG
    
end
%% Plot 

% get average power over trials
old_TF = TF;
LWPC_I_I_nphase = TF.power_phase.LWPC_MI_I - TF.power_non_phase.LWPC_MI_I;
LWPC_I_C_nphase = TF.power_phase.LWPC_MI_C - TF.power_non_phase.LWPC_MI_C;
av_eegpower = TF.power_phase.LWPC_MI_C -TF.power_phase.LWPC_MI_I;
chan2use = 'Pz';
chanidx = strcmpi(chan2use,{EEG.chanlocs.labels});
power_sub = squeeze(av_eegpower(chanidx, :, :));
%itpc_sub = squeeze(itpc(chanidx, :, :));

max_scale = max(power_sub,[], 'all') +0.5
min_scale = min(power_sub,[], 'all') -0.5
figure(6)
contourf(New_trial_time,freq_range,power_sub,60,'linecolor','none')
set(gca,'clim',[min_scale, max_scale],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)
colorbar


av_eegpower_old = old_TF.power.LWPC_MI_I;
power_sub_old = squeeze(av_eegpower_old(chanidx, :, :));

figure(4)
contourf(New_trial_time,freq_range,power_sub_old,20,'linecolor','none')
set(gca,'clim',[-5, 2],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Power')
colormap(jet)


%'xlim',[-200 800],
figure(3)
subplot(121)
contourf(New_trial_time,freq_range,power_sub,40,'linecolor','none')
set(gca,'clim',[-0.5 0.5],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
title('Logarithmic frequency scaling')


% subplot(122)
% contourf(EEG.times,freq_range,itpc_sub,40,'linecolor','none')
% set(gca,'clim',[0 .6],'xlim',[-200 1000])
% xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
% title('ITPC')


%% Tasks to be done
% remove the ERP from the data - to get non-phase locked pertubations

% Use the whole data baseline if possible - seems to be a computational
% problem as the dataset is too large to be processed in this manner though

% With regard to the baseline correction we will continue as follows.
% baseline correction so far was done incorrectly. We do not use the
% baseline of each trial but we HAVE to take the baseline average. Hence we
% continue as follows. We will keep all frequency transformed trials. Of
% these frequency transformed trials we can calculate the average baseline
% over all trials - covering the baseline activity from pre-stim to whole
% trial freqeuncy activy. This will be done in a second step. Hence we will
% have single trials which we can subsequently analyze using a linear
% regression. Moreover, keep the averaged ERSP image. So we will have
% single trials on which we can perform our regression analysis and in
% addition, we will have the correctly avereaged (prior to log transfrom)
% ERSP which we can use for visualization purposes.

% perhaps focus on a sub analysis in the theta band to look at scalp
% distribution for all electrodes midfrontal (Fz, FC1, FCz, FC2, Cz) or
% (AFz, F1, Fz, F2, FC1, FCz, FC2.), for the latter see https://www.nature.com/articles/s41598-021-95631-1

% Select a proper baseline - compare baseline analyses with a baseline over
% the whole trial and a baseline period during fixation. Use a baseline
% that goes over the whole task dataset - not a condition specific baseline
