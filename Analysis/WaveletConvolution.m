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
% signal back. WHy? Frequency convotuion is fast where it is slow in the
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
dirs.home                    = 'G:\Julia Ficke\AdaptiveControl\CleanData_EEG'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'E:\AdaptiveControl\Data\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden

%git directory name for me
%  /C/Users/doex9445/Dateien/Julius/AdaptiveControl/tVNS-Project

% Participant IDs we want to analyze            
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'

%remove PD 16, index 68 - no dat available
Participant_IDs(68)          = []; 
Participant_IDs              = {Participant_IDs(:).name};  
Part_N                       = length(Participant_IDs); %number of participants

% event onset codes
onset =  {  'S 21'  'S 22'  'S 23'  'S 24'  'S 41'  'S 42'  'S 43' ...
    'S 44'  'S 61'  'S 62'  'S 63'  'S 64'  'S 81'  'S 82'  'S 83'  'S 84'  };
epoch_dur = [-1.3  2];


% INSERT PARTICIPANTS WE WISH NOT TO ANALYZE

%add EEGLAB path and start the program
addpath(eegl);
eeglab

%% Wavelets

% Load one dataset to get parameters for wavelet anaylsis
sub                          = 3;
fileID                       = strcat(Participant_IDs{sub}, '_task-stroop_eeg_pruned_auto.set'); %get file ID
folderID                     = fullfile(dirs.home,Participant_IDs{sub}, 'eeg' );%get folder ID
EEG                          = pop_loadset('filename', fileID,'filepath',[folderID]); % load file


% Define params wavelet
length_wavelet               = 4; % wavelet lengths will be 4 seconds
wave_pnts                    = length_wavelet* EEG.srate;
cycle_num                    = 4; % number of wavelet cycles
freq_up                      = 60; % upper frequency limit in Hz
freq_low                     = 4; %lower freqeuncy limit in Hz
freq_num                     = 30; %number of freqeuncies to be estimated
freq_range                   = logspace(log10(freq_low), log10(freq_up), freq_num); %range of frequencies


% time vector for wavlet (-2:2s)
time                         = -wave_pnts/EEG.srate/2: 1/EEG.srate : wave_pnts/EEG.srate/2;

%vector of standard deviations per frequency
s                            = cycle_num./(2*pi*freq_range);

% create wavelets
wavelets                     = zeros(freq_num, wave_pnts + 1);
for fr = 1:freq_num
    c_sine = exp(1i*2*pi*freq_range(fr).*time);
    gaus = exp(-(time.^2)./(2*(s(fr)^2)));
    wavelets(fr,:) = c_sine.* gaus;
end


% briefly inspect wavelets
% figure
% plot(time, real(wavelets(2,:)))
%     

%% Perform Wavelet analysis on participant data

for sub = 2:Part_N
    % load data set
    fileID              = strcat(Participant_IDs{sub}, '_task-stroop_eeg_pruned_auto.set'); %get file ID
    folderID            = fullfile(dirs.home,Participant_IDs{sub}, 'eeg' );%get folder ID
    EEG                 = pop_loadset('filename', fileID,'filepath',[folderID]); % load file
    
    n_data_sub1               = EEG.pnts;
    
    % define convolution parameters (used from example code by Mike C. Cohen)
    n_wavelet            = length(time);
    n_data               = EEG.pnts;
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;


    
%     % get fft of data - later i might just concatenate all channels to one
%     % array
%     fft_dat             = zeros(EEG.nbchan, n_conv_pow2);
%     for chan            = 1:EEG.nbchan
%         fft_dat(chan,:) = fft(EEG.data(chan,:), n_conv_pow2);
%     end
    
    % initialize
    eegpower            = zeros(EEG.nbchan, freq_num , EEG.pnts);
    itpc                = zeros(EEG.nbchan, freq_num , EEG.pnts);
    
    fprintf('Performing wavelet convolution on participant %s. \n',Participant_IDs{sub})
    
    
    for chan = 1: EEG.nbchan
        for freq = 1:freq_num
            % get fft transform of wavelet
            fft_wavelet = fft(wavelets(freq,:),n_conv_pow2);
            
            % get fft of data
            fft_dat = fft(EEG.data(chan,:), n_conv_pow2);
            
            % convolve data and wavelet
            decomp = fft_dat.*fft_wavelet;
            
            % transform back
            decomp = ifft(decomp);
            decomp = decomp(1:n_convolution);
            decomp = decomp(half_of_wavelet_size + 1:end - half_of_wavelet_size);
            % take mean of all trials, compute magnitude and square
            %decomp = reshape(decomp, [EEG.pnts,EEG.trials]);
            
            
            % extract ITPC
            itpc(chan, freq, :) = abs(mean(exp(1i*angle(decomp)),2));
            eegpower(chan, freq,:) = mean(abs(decomp),2).^2;
        end
    end
    % save power and itpc data, baseline correction will be performed
    % afterwards
    
    fprintf('Saving frequency data of participant %s. \n',Participant_IDs{sub})
    
    save_name = strcat(Participant_IDs{sub}, '_frequency_data.mat');
    save_loc = fullfile(dirs.eegsave, save_name);
    save(save_loc,'itpc','eegpower','-v7.3');
    
    fprintf('Frequency data of participant %s has been saved. \n',Participant_IDs{sub})
    
    
    
    %remove intermediate variables to make some space
    clear EEG itpc eegpower decomp fft_dat fft_wavelet
    
end
%    Baseline correction will be executed elsewhere
% %temppower = 10*log10(temppower./mean(temppower(bsidx(1):bsidx(2))));
            
% 
% % Baseline correction between -500ms and -200 ms
%     bsidx = dsearchn(EEG.times',[-500 -200]');
%     baseline = mean(eegpower(:,:,[bsidx(1):bsidx(2)]),3);
%     %
%     % calculate power in decibel relative to baseline period
%     for chan = 1:EEG.nbchan
%         base = baseline(chan,:)'; %electrode specific baseline
%         eegdata = squeeze(eegpower(chan,:,:)); % electrode specific data
%         eegpower(chan,:,:) = 10*log10(eegdata./base); % eeg power with respect to baseline activity
%     end
%     
% % Plot 
% chan2use = 'pz';
% chanidx = strcmpi(chan2use,{EEG.chanlocs.labels});
% power_sub = squeeze(eegpower(chanidx, :, :));
% itpc_sub = squeeze(itpc(chanidx, :, :));
% 
% figure(4)
% subplot(121)
% contourf(EEG.times,freq_range,power_sub,40,'linecolor','none')
% set(gca,'clim',[-3 3],'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(freq_low),log10(freq_up),6),'yticklabel',round(logspace(log10(freq_low),log10(freq_up),6)*10)/10)
% title('Power')
% 
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