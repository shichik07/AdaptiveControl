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

clear all; close all; clc
dbstop if error

eegl                         = '/home/jules/Dropbox/PhD_Thesis/EEG_Labor/EEG_Software/eeglab2021.1';
% set directories
dirs.home                    = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = '/media/jules/DriveJules/AdaptiveControl/Data/FrequencyData/'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = '/home/jules/Dropbox/PhD_Thesis/Adaptive_Control/Analysis/Analysis/AdaptiveControl/Analysis';
addpath(dirs.functions) 



%git directory name for me
%  /C/Users/doex9445/Dateien/Julius/AdaptiveControl/tVNS-Project

% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

%remove PD 16, index 68 - no data available and was excluded from the
%analysis
Participant_IDs(68)          = []; 
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
% cd('/home/jules/Dropbox/PhD_Thesis/Studying/NeuralTimesSeries')
% load sampleEEGdata
% 

% Define params wavelet
length_wavelet               = 4; % wavelet lengths will be 4 seconds
wave_pnts                    = length_wavelet* EEG.srate;
cycle_num                    = 6; % number of wavelet cycles
freq_up                      = 40; % upper frequency limit in Hz
freq_low                     = 6; %lower freqeuncy limit in Hz
freq_num                     = 30; %number of freqeuncies to be estimated
freq_range                   = logspace(log10(freq_low), log10(freq_up), freq_num); %range of frequencies
limit                        = 262144; % limit of what matrix operations my PC is capable of performing - important for zero padding
bin_nr                       = 10; %number of bins in which we analyze the data

% time vector for wavlet (-2:2s)
time                         = -wave_pnts/EEG.srate/2: 1/EEG.srate : wave_pnts/EEG.srate/2;
%vector of standard deviations per frequency
s                            = cycle_num./(2*pi*freq_range); % single cycle wavelet
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
Baseline_time = dsearchn(EEG.times',[-300 -100]'); 
keep_time =  dsearchn(EEG.times',[-200 1000]'); % time range we want to keep in the final dataset
New_trial_time = -200/1000: 1/EEG.srate : 1000/1000;



% briefly inspect wavelets
 figure (5)
 plot(time, real(wavelets(1,:)))
%     

%% Perform Wavelet analysis on participant data

for sub = 1:Part_N
    %load data set
    fileID              = strcat(Participant_IDs{sub}, '_epoched_freq.set'); %get file ID
    folderID            = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    EEG                 = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    Items = get_trlindices(EEG); %get indices of trials by condition
    Fnames = fieldnames(Items); % get associated condition names
    
    % create struct where we save our averaged frequency condition data
    vals = zeros(EEG.nbchan, freq_num , length(keep_time(1):keep_time(2)));
    vals = repmat({vals},1,length(Fnames));
    args=[Fnames';vals];
    TF.power = struct(args{:}); % create a struct with matrix for each category for power values
    TF.itpc  = struct(args{:}); % create a struct with matrix for each category for inter trial phase clustering values
    args2 = [Fnames';num2cell(structfun(@numel,Items))']; %get numer of items in case of weighting
    TF.Item_nr = struct(args2{:}); 
    TF.chanlocs = EEG.chanlocs; % keep channel information
    
    fprintf('Performing wavelet convolution on participant %s. \n',Participant_IDs{sub})
   
    for con = 1:size(Fnames,1)
        trl_indices = getfield(Items, Fnames{con});
        
        fprintf('Analyzing trial subset %s of participant %s. \n',Fnames{con}, Participant_IDs{sub})
        % get bin data and parameter
        con_data = EEG.data(:,:,trl_indices);
        con_trl = size(con_data,3); % nr of trials for this category
        
        % reshape trials into one vector
        con_data = reshape(con_data, [EEG.nbchan, EEG.pnts*con_trl]);
        
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
        for chan = 1: EEG.nbchan
            for freq = 1:freq_num
                % reshape trials into one vector
                dat = reshape(squeeze(con_data(chan,:,:)), [1, EEG.pnts*con_trl]);
        
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
                TF.itpc.(Fnames{con})(chan,freq,:) = itpc(keep_time(1):keep_time(2)); %only keep period we are interested in
                
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
                    base_power(trl) = mean(decomp(keep_time(1): keep_time(2), trl)); % baseline over trial period
                    %base_power(trl) = mean(decomp(Baseline_time(1): Baseline_time(2), trl)); %baseline during fixation
                end
                
                % perform baseline correction and convert to decible scale
                for trl = 1:con_trl
                    temppower(:,trl) = 10*log10(decomp(:,trl)./base_power(trl));
                end
                
                %save in matrix
                TF.power.(Fnames{con})(chan,freq,:) = mean(temppower(keep_time(1):keep_time(2),:),2);
            end
            clear decomp temppower fft_dat fft_wavelet dat
        end
    end
    % save power and itpc data, baseline correction will be performed
    % afterwards
    fprintf('Saving frequency data of participant %s. \n',Participant_IDs{sub})
    
    save_name = strcat(Participant_IDs{sub}, '_frequency_data_phaselocked.mat');
    save_loc = fullfile(dirs.eegsave, Participant_IDs{sub}, save_name);
    save(save_loc,'TF','-v7.3');
    
    fprintf('Frequency data of participant %s has been saved. \n',Participant_IDs{sub})
    
    %remove intermediate variables to make some space
    clear EEG TF
    
end
%% Plot 

% get average power over trials
old_TF = TF;
av_eegpower = TF.power.LWPC_MI_I;
chan2use = 'Pz';
chanidx = strcmpi(chan2use,{EEG.chanlocs.labels});
power_sub = squeeze(av_eegpower(chanidx, :, :));
%itpc_sub = squeeze(itpc(chanidx, :, :));

max_scale = max(power_sub,[], 'all') +0.5
min_scale = min(power_sub,[], 'all') -0.5
figure(8)
contourf(New_trial_time,freq_range,power_sub,20,'linecolor','none')
set(gca,'clim',[-5, 2],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
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

% restrict Analysis to midfrontal channels for theta analysis for now.
% There is no reason to perform analysis on all channels

% epoch datasets with sufficient margins to have good time frequency
% resolution

% remove the ERP from the data - to get non-phase locked pertubations

% always only use a chunk of trials small enough not to overload the
% workspace (for fft we zeropad to the power of two, let us keep this
% smaller for more efficient analyses)

% perhaps focus on a sub analysis in the theta band to look at scalp
% distribution for all electrodes midfrontal (Fz, FC1, FCz, FC2, Cz) or
% (AFz, F1, Fz, F2, FC1, FCz, FC2.), for the latter see https://www.nature.com/articles/s41598-021-95631-1

% restore online reference - which means the analysis has to be done
% again... that might be a task for later


A = zeros(9,2,3);
B = ones (9,2);
C = ones (9,2)+1;

A(:,:,2) = B;
A(:,:,3) = C;
B = reshape(A,[9,6]);

B = reshape(B,[9,2,3]);

A = [9,6,3;4,2,8;1,0,0]; 
B = [3,2,1];

bs = base_power';
temppower = bsxfun(@rdivide,decomp,bs(:)');
temppower = decomp./base_power';


a = temppower(:,1);
b = temppower(:,2);
c = decomp(:,1)./base_power(1);
d = decomp(:,2)./base_power(2);

mean(a == c)
mean(b == d)

a = temppower(keep_time(1):keep_time(2),:);