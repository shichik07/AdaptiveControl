function ftdata = get_pseudo_FT_struct_cluster(data, Infos)
% Transform our preprocessed data into a structure that can be read by the
% spm_eeg_ft2spm.m function. Takes the preprocessed data that we used and
% outputs it into a FT like data structure that can be read by spm.
ftdata = [];

ftdata.label = {Infos.chanlocs(:).labels}; % 1 x elec ; type: cell
ftdata.dimord = 'subj_chan_freq_time';
ftdata.freq = Infos.freqs; % 1 x freq 
ftdata.time = Infos.time; % 1 x time 

% before we can include the power matrix, we have to permute it to get the
% trial x channel x frequency x time structure
ftdata.powspctrm =   permute(data,[4 1 2 3]);; % subj x chan x freq x time

% create a header section because the spm function reads that taken from
% the read_eeglabheader() function by 2008 Arnaud Delorme.

ftdata.hdr.Fs          =  1/(Infos.time(2)- Infos.time(1)); %sampling rate
ftdata.hdr.nChans      =  length(ftdata.label); %number of channels
ftdata.hdr.nSamples    =  length(ftdata.time); %number of samples
ftdata.hdr.nSamplesPre = -min(ftdata.time)* ftdata.hdr.Fs; %prestim sample number
% 
ftdata.hdr.label        = ftdata.label;

ind = 1;
for i = 1:length(Infos.chanlocs)
    ftdata.hdr.elec.label{ind, 1} = Infos.chanlocs(i).labels;
    % this channel has a position
    ftdata.hdr.elec.elecpos(ind,1) = Infos.chanlocs(i).X;
    ftdata.hdr.elec.elecpos(ind,2) = Infos.chanlocs(i).Y;
    ftdata.hdr.elec.elecpos(ind,3) = Infos.chanlocs(i).Z;
    ind = ind+1;
    
end
ftdata.elec = ftdata.hdr.elec;

end