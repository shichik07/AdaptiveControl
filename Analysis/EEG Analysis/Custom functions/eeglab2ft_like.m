function X = eeglab2ft_like(eeg ,ana_type)
% function to get my preprocessed EEG data into spm with fieldtrip-"like"
% data. Takes the EEG data struct from eeglab and turns it into a struct
% that can be processed by the spm_eeg_ft2spm() function. 

ftdata = [];

% if we want to perform a timelocked analysis
if ana_type == 'tl' 
    for i = 1:size(eeg.data,3)
        ftdata.trial{i} = eeg.data(:, :, i);
        ftdata.time{i} = eeg.xmin*1000:eeg.pnts/eeg.srate:eeg.xmax*1000;
    end
end

% if we want to perform a powerspectrum analysis
elseif ana_type == 'pw'
    for i = 1:size(eeg.data,3)
        ftdata.powspctrm{i} = eeg.data(:, :, :, :);
        ftdata.time{i} = eeg.xmin*1000:eeg.pnts/eeg.srate:eeg.xmax*1000;
    end


ftdata.fsample = eeg.srate; %sampling rate
ftdata.label = {eeg.chanlocs.labels}'; %electrode labels

elif
