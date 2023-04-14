%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% and Julia Ficke (julia.ficke@uni-oldenburg.de) %
% Preprocessing corrected                        %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc
dbstop if error


% set directories
dirs.home                    = 'G:\Julia Ficke\AdaptiveControl\BIDS Format'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'H:\KricheldorffJ\AdaptiveControl\PreprocessedNew\FrequencyData'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden
dirs.functions               = 'C:\Users\doex9445\Dateien\Julius\AdaptiveControl\AdaptiveControl\Analysis';

% Participant IDs we want to analyze
Participant_IDs              = dir(dirs.home);
Participant_IDs              = Participant_IDs([Participant_IDs(:).isdir]); % remove all files (isdir property is 0)
Participant_IDs              = Participant_IDs(~ismember({Participant_IDs(:).name},{'.','..'}));% remove '.' and '..'
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

cd('G:\Julia Ficke\AdaptiveControl\CleanData_EEG')
eegl                         = 'C:\Program Files\MATLAB\EEGSoftware';
% set directories
dirs.home                    = 'G:\Julia Ficke\AdaptiveControl\BIDS Format'; %hier habe ich das Gruppenlaufwerk gespeichert - du müsstest hier deinen Speicherort für die Daten eintragen
dirs.eegsave                 = 'G:\Julia Ficke\AdaptiveControl\CleanData_EEG'; % hier Ordner zum Speichern der Ergebenisse - im Gruppenlaufwerk unter PipelineValidate zu finden


% Trigger Namen
Block.onset  = {'S  2', 'S  3', 'S  4', 'S  5', 'S  6', 'S  7', 'S  8'};
Block.offset = {'S101', 'S102', 'S103', 'S104', 'S105', 'S106', 'S107'};

% add fieldtrip and set defaults
addpath(eegl);
eeglab

for sub = 1:Part_N
    ALLCOM = []; ALLEEG = []; EEG = []; STUDY = [];
    
    %construct folder name to load the dataset
    folderID = fullfile(dirs.home,'BIDS',Participant_IDs{sub}, 'eeg');
    fileID = dir(fullfile(folderID, '*.vhdr')); %find header name
    
    
    % Load EEG Data
    EEG = pop_loadbv(folderID, fileID.name);
    
    % Edit Channel Labels
    EEG=pop_chanedit(EEG, 'settype',{'1:127','EEG'},'settype',{'128:129','EOG'});
    
    
    %% Codierung der Trigger/Events  (Welche Benennungen werden gebraucht)
    % Es werden die einzelnen Spalten für die Informationen zu den Events erstellt
    
    %Länge der neuen Spalten wird bestimmt
    number_event = length(EEG.event)
    
    %renaming event Onset ( if Response is incorrect, the onset Event should be
    %renamed)
    for ii=1:(number_event-1)
        
        switch EEG.event(ii+1).type
            case 'S131'
                EEG.event(ii).type = 'S 31';
            case 'S132'
                EEG.event(ii).type = 'S 32';
            case 'S133'
                EEG.event(ii).type = 'S 33';
            case 'S134'
                EEG.event(ii).type = 'S 34';
            case 'S151'
                EEG.event(ii).type = 'S 51';
            case 'S152'
                EEG.event(ii).type = 'S 52';
            case 'S153'
                EEG.event(ii).type = 'S 53';
            case 'S154'
                EEG.event(ii).type = 'S 54';
            case 'S171'
                EEG.event(ii).type = 'S 71';
            case 'S172'
                EEG.event(ii).type = 'S 72';
            case 'S173'
                EEG.event(ii).type = 'S 73';
            case 'S174'
                EEG.event(ii).type = 'S 74';
            case 'S191'
                EEG.event(ii).type = 'S 91';
            case 'S192'
                EEG.event(ii).type = 'S 92';
            case 'S193'
                EEG.event(ii).type = 'S 93';
            case 'S194'
                EEG.event(ii).type = 'S 94';
                
        end
        
    end
    
    part  = cell(number_event,1);  % onset, response
    analysetype = cell(number_event,1);  % MI, MC, main_incon, main_con
    congruency  = cell(number_event,1);  % congruent, incongruent
    trial     = cell(number_event,1);  % inducer, diagnostic
    answer = cell(number_event,1);  % correct, incorrect
    
    % Es wird die Benennung durch SXXX codiert
    for i=1:number_event
        if strcmp(EEG.event(i).code, 'Stimulus')
            digit1 = str2double(EEG.event(i).type(2));
            digit2 = str2double(EEG.event(i).type(3));
            digit3 = str2double(EEG.event(i).type(4));
            
            
            %definition für part (onset, response, fixation dot)
            
            if isnan(digit1)
                part{i} = 'Onset';
                answer{i} = 'Onset';
            end
            
            switch digit1
                case 1
                    part{i} = 'Response';
                case 2
                    part{i} = 'Fixation';
            end
            
            
            if isnan(digit2) || isnan(digit2) || isnan(digit2)
                part{i}     = 'Start Block';
                analysetype{i} = 'Start Block';
                congruency{i}     = 'Start Block';
                trial{i} = 'Start Block';
                answer{i} = 'Start Block';
                continue
            end
            
            if digit2==0 && digit1 == 1
                part{i}     = 'End Block';
                analysetype{i} = 'End Block';
                congruency{i}     = 'End Block';
                trial{i} = 'End Block';
                answer{i} = 'End Block';
                continue
            end
            
            %defintion of answer (correct, false)
            
            switch digit2
                case {2, 4, 6, 8}
                    answer{i} = 'correct';
                case {3, 5, 7, 9}
                    answer{i} = 'incorrect';
            end
            
            
            %defintion of congruency
            
            if strcmp(part{i}, 'Fixation')
                % the interpretation of digit2 is not given for Fixation
                congruency{i} = 'Fixation'; % incongruent or congruent
                trial{i} = 'Fixation'; % inducer or diagnostic
                analysetype{i} = 'Fixation'; % MI, MC, main_incon, main_con
                answer{i} = 'Fixation';
            else
                switch digit3
                    case {1, 3}
                        congruency{i} = 'incongruent';
                    case {2, 4}
                        congruency{i} = 'congruent';
                end
                
                switch digit3
                    case {1, 2}
                        trial {i} = 'inducer';
                    case {3, 4}
                        trial{i} = 'diagnostic';
                end
            end
            
            %defintion of analyse type (MI, MC, main_incon, main_con)
            switch digit2
                case {2, 3}
                    analysetype{i} = 'MI';
                case {4, 5}
                    analysetype{i} = 'MC';
                case {6, 7}
                    analysetype{i} = 'main_incon';
                case {8, 9}
                    analysetype{i} = 'main_con';
            end
            
            
        end
    end % for
    
    
    % die neu erstellten Spalten werden in die Tabelle ALL EEG eingefügt
    for i=1:number_event
        EEG.event(i).analysetype = analysetype (i)
        EEG.event(i).part = part (i)
        EEG.event(i).congruency = congruency (i)
        EEG.event(i).trial = trial (i)
        EEG.event(i).answer = answer (i)
    end
    
    
    %% Resample to 250 Hz
    EEG = pop_resample( EEG, 250);
    
    
    %remove EOG Elektroden
    EEG = pop_select( EEG, 'nochannel',{'HEOG','VEOG'});
    
    %add channel FCz
    EEG = pop_chanedit(EEG, 'insert',64,'changefield',{64,'labels','FCz'},'changefield',{64,'theta','0'},'changefield',{64,'radius','0.127777777777778'},'changefield',{64,'X','0.390731128489274'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.920504853452440'},'changefield',{64,'sph_theta','0'},'changefield',{64,'sph_phi','67'},'changefield',{64,'sph_radius','1'},'changefield',{64,'type','EEG'},'setref',{'1:128','FCz'});
    
    
    
    %Save Data
    pdest1 = fullfile(dirs.eegsave, Participant_IDs{sub}, 'eeg')
    mkdir (pdest1)
    EEG = pop_saveset(EEG,'filename', strcat(fileID.name(1:end-5), '_testset.set'),'filepath',[pdest1]);
    
    
    
    
    % Reject data periods not containing experimental data
    
    for ii = 1:length(Block.onset)
        try
            Blockstart = find(ismember({EEG.event.type},Block.onset{ii})) +1;
            Blockstop = find(ismember({EEG.event.type},Block.offset{ii}));
            EEG = eeg_eegrej(EEG, [EEG.event(Blockstop).latency+(1.5*EEG.srate) EEG.event(Blockstart).latency-(1.5*EEG.srate)]); % +-500 samples um 2 Sekunden Daten extra zu behalten
        end
    end
    
    % We also reject data first instance to last epoch and last to end
    % recording6
    try
        LastBlock = find(ismember({EEG.event.type},'S108'));
        EEG = eeg_eegrej(EEG, [EEG.event(LastBlock).latency+(1.5*EEG.srate) EEG.pnts]);
    end
    try
        FirstBlock = find(ismember({EEG.event.type},'S  1')) +1;
        EEG = eeg_eegrej(EEG, [EEG.event(1).latency EEG.event(FirstBlock).latency-(1.5*EEG.srate)]); % again +-500 samples um 2 Sekunden Daten extra zu behalten
    end
    % Highpass-filter data with 0.1 Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'chantype',{'EEG'});
    
    %detrend data
    EEG.data(:,:) = detrend(EEG.data(:,:)')'
    
    
    % Reject electrodes automatically and store index of rejected electrodes
    % under EEG.reject.indelec so we can use the indices later for
    % interpolation. For example EEG = eeg_interp(EEG,EEG.reject.indelec)
    
    EEGold = EEG;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',5,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
    
    % if there are no rejected channels, skip this part
    try
        
        EEG.reject.indelec = EEG.etc.clean_channel_mask;
        % name rejected channels
        EEG.reject.remChan = {EEGold.chanlocs(~EEG.reject.indelec).labels};
        fprintf('Channel %s was removed. \n',EEG.reject.remChan{:})
    end
    % write the FCz Channel back into the datafile
    EOG = [find(ismember({EEG.chanlocs.type},{'EOG'}))]
    
    %rereference data with average
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'FCz'},'sph_radius',{1},'sph_theta',{0},'sph_phi',{67},'theta',{0},'radius',{0.12778},'X',{0.39073},'Y',{0},'Z',{0.9205},'type',{'EEG'},'ref',{'FCz'},'urchan',{[]},'datachan',{0}));
    
    % Remove line noise with cleanline
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',2.5,'winstep',1);
    
    % Use checkset to see if we messed something up
    EEG = eeg_checkset(EEG);
    
    % Save data at this point at a location of your choice
    EEG = pop_saveset(EEG, 'filename', strcat(fileID.name(1:end-5), '_imported'),'filepath',[pdest1]);
    
end
%% Automatic Noise rejection - Noise entfernen


for sub = 1:length(Part_N)
    ALLCOM = []; ALLEEG = []; EEG = []; STUDY = [];
    
    %construct folder name to load the dataset
    folderID = fullfile(dirs.home,'BIDS',Participant_IDs{sub}, 'eeg');
    fileID = dir(fullfile(folderID, '*.vhdr')); %find header name
    pdest1 = fullfile(dirs.eegsave, Participant_IDs{sub}, 'eeg')
    
    EEG = pop_loadset('filename', strcat(fileID.name(1:end-5), '_imported.set'),'filepath',[pdest1]);
    
    %Perform ASR
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off' ,'BurstCriterion',80,'WindowCriterion','off','BurstRejection','on','Distance','Euclidian');
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    EEG = pop_saveset(EEG, 'filename', strcat(fileID.name(1:end-5), '_noise_removed_automatic'),'filepath',[pdest1]);
end

eeglab redraw
%% Perform ICA

for sub = 1:length(Part_N)
    ALLCOM = []; ALLEEG = []; EEG = []; STUDY = [];
    
    %construct folder name to load the dataset
    folderID = fullfile(dirs.home,'BIDS',Participant_IDs{sub}, 'eeg');
    fileID = dir(fullfile(folderID, '*.vhdr')); %find header name
    pdest1 = fullfile(dirs.eegsave, Participant_IDs{sub}, 'eeg');
    
    % Load Data
    EEG = pop_loadset('filename', strcat(fileID.name(1:end-5), '_noise_removed_automatic.set'),'filepath',[pdest1]);
    
    % Perform ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', {'EEG'},'extended',1,'interrupt','on'); %all chans except eye chans
    
    
    % Save
    EEG = pop_saveset(EEG, 'filename', strcat(fileID.name(1:end-5), '_decomposed_auto'),'filepath',[pdest1]);
    
    
end

%% Manuell Components entfernen und Daten Set als '_pruned.set' abspeichern
% Hier werden einmal Epochen erstellt um besser über die Componenten zu
% entscheiden. Es werden sich die ersten 35 Componenten bei IC label
% angeschaut und die Störfaktoren rausgeschrieben


%% Interpolate Channels - after ICA rejection and epoch data

%Quelle: https://sccn.ucsd.edu/pipermail/eeglablist/2016/011616.html#:~:text=To%20interpolate%20channels%20you%20would%20load%20up%20one,list%20is.%20Go%20to%20Tools%20%3E%20Interpolate%20Channels.


for sub = 1:length(Part_N)
    ALLCOM = []; ALLEEG = []; EEG = []; STUDY = [];
    
    %construct folder name to load the dataset
    pdest1 = fullfile(dirs.eegsave, Participant_IDs{sub}, 'eeg')
    folderID = fullfile(dirs.home,'BIDS',Participant_IDs{sub},'eeg');
    fileID = dir(fullfile(folderID, '*.vhdr')); %find header name
    
    
    
    EEG2 = pop_loadset('filename', strcat(fileID.name(1:end-5), '_testset.set'),'filepath',[pdest1]);
    EEG = pop_loadset('filename', strcat(fileID.name(1:end-5), '_pruned.set'),'filepath',[pdest1]);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    %Interpolate removed channels
    EEG = pop_interp(EEG, EEG2.chanlocs, 'spherical');
    
    %Extract epochs (Onset coded)
    
    EEG = pop_epoch( EEG, { 'S 21'  'S 22' 'S 23'  'S 24' 'S 41'  'S 42' 'S 43'  'S 44' 'S 61'  'S 62' 'S 63'  'S 64' 'S 81'  'S 82' 'S 83'  'S 84'  }, [-0.2           1.0], 'newname', ' resampled pruned with ICA epochs', 'epochinfo', 'yes');
    
    % Save
    EEG = pop_saveset(EEG, 'filename', strcat(fileID.name(1:end-5), '_epoched_freq'),'filepath',[pdest1]);
end
    
  