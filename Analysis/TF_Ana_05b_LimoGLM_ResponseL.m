%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% Adaptive Control Frequency Analysis            %
% Julius Kricheldorff(julius.kricheldorff@uol.de)%
% Statistical Analysis using Limo                %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'derivatives'}));
Participant_IDs              = Participant_IDs(~contains({Participant_IDs(:).name},{'sub-CY'})); %remove Data of young participants
Participant_IDs              = {Participant_IDs(:).name};
Part_N                       = length(Participant_IDs); %number of participants

addpath(genpath(LIMO))
addpath(eegl)
addpath(dirs.functions)


for sub = 2:Part_N
    
    % switch to sub directory
    cd(fullfile(dirs.home, Participant_IDs{sub}))
    
    % get file IDs
    folderID                    = fullfile(dirs.home,Participant_IDs{sub});%get folder ID
    folderID_t1                 = dir(fullfile(folderID, '**', '*_pl_RL_bslC_RL.mat'));
    folderID_t2                 = dir(fullfile(folderID, '**', '*_npl_RL_bslC_RL.mat'));
    LWPC_t1                     = folderID_t1(contains({folderID_t1(:).name},{'LWPC'}));
    ISPC_t1                     = folderID_t1(contains({folderID_t1(:).name},{'ISPC'}));
    LWPC_t1                     = {LWPC_t1(:).name}; %phase_locked files LWPC
    ISPC_t1                     = {ISPC_t1(:).name}; %phase_locked files ISPC
    LWPC_t2                     = folderID_t2(contains({folderID_t2(:).name},{'LWPC'}));
    ISPC_t2                     = folderID_t2(contains({folderID_t2(:).name},{'ISPC'}));
    LWPC_t2                     = {LWPC_t2(:).name}; %non_phase_locked files LWPC
    ISPC_t2                     = {ISPC_t2(:).name}; %non_phase_locked files ISPC
    
    
    % first we load all the data for one condition (for now we will do this
    % only for the listwide effect and the phase-locked data)
    
    fprintf('Loading data of participant %s. \n',Participant_IDs{sub})
    
    for effect = 1:4
        if effect == 3
            files = LWPC_t2;
            file_type = 'TF_non_phase';
            trial_counts = zeros(1, length(files));
            data = [];
            for con = 1:length(files)
                load_loc = fullfile(dirs.home, Participant_IDs{sub},  files{con});
                clear TF_non_phase % remove data to save up space
                load(load_loc, 'TF_non_phase');
                data = cat(4, data, TF_non_phase.power);
                trial_counts(con) = size(TF_non_phase.power,4);
                save_n = strcat(Participant_IDs{sub}, '_GLM_LWPC_npl_RL.mat');
            end
        elseif effect == 4
            files = LWPC_t1;
            file_type = 'TF_phase';
            trial_counts = zeros(1, length(files));
            data = [];
            for con = 1:length(files)
                load_loc = fullfile(dirs.home, Participant_IDs{sub},  files{con});
                clear TF_phase % remove data to save up space
                load(load_loc, 'TF_phase');
                data = cat(4, data, TF_phase.power);
                trial_counts(con) = size(TF_phase.power,4);
                save_n =  strcat(Participant_IDs{sub}, '_GLM_LWPC_pl_RL.mat');
            end
        elseif effect == 1
            files = ISPC_t2;
            file_type = 'TF_non_phase';
            trial_counts = zeros(1, length(files));
            data = [];
            for con = 1:length(files)
                load_loc = fullfile(dirs.home, Participant_IDs{sub},  files{con});
                clear TF_non_phase % remove data to save up space
                load(load_loc, 'TF_non_phase');
                data = cat(4, data, TF_non_phase.power);
                trial_counts(con) = size(TF_non_phase.power,4);
                save_n =  strcat(Participant_IDs{sub}, '_GLM_ISPC_npl_RL.mat');
            end
        elseif effect == 2
            files = ISPC_t1;
            file_type = 'TF_phase';
            trial_counts = zeros(1, length(files));
            data = [];
            for con = 1:length(files)
                load_loc = fullfile(dirs.home, Participant_IDs{sub},  files{con});
                clear TF_phase % remove data to save up space
                load(load_loc, 'TF_phase');
                data = cat(4, data, TF_phase.power);
                trial_counts(con) = size(TF_phase.power,4);
                save_n = strcat(Participant_IDs{sub}, '_GLM_ISPC_pl_RL.mat');
            end
        end
        
        % make limo directory
        limo_dir = 'LIMO_Output';
        if ~exist(limo_dir)
            mkdir(limo_dir);
        end
        
        % get vecotr for categories
        condition_order = {};
        for con = 1:length(files)
            if contains(files{con}, 'MC_C')
                condition_order{con} = 'MC_C';
            elseif contains(files{con}, 'MC_I')
                condition_order{con} = 'MC_I';
            elseif contains(files{con}, 'MI_C')
                condition_order{con} = 'MI_C';
            elseif contains(files{con}, 'MI_I')
                condition_order{con} = 'MI_I';
            end
        end
        
        %         Y                  = squeeze(data(:,1,:,:)); %EEG data with format electrodes x frames x trials/subjects
        %         Cat                = conditions; %vector describing the different conditions (if no conditions Cat = 0)
        %         Cont               = 0; % matrix describing the different covariates (if no covariates Cont = 0)
        %         directory          = fullfile(dirs.home, Participant_IDs{sub}, limo_dir);%path of folder where the outputs will be saved (see below)
        %         zscoring           = 1 ; % we only use categorical regressors [0/1] - %if 1 (default) continuous regressors are zscored, which means that the betas coefficients will have units micro-volts per std of the predictor variable.
        %         full_factorial     = 1 ;%[0/1] - %if 1 create interaction terms between the factors described in Cat
        %         flag               = 1; %figure on/off [1/0]
        %         %chanlocs           = % we can also put in chanlocs?
        %         [X, nb_conditions, nb_interactions, nb_continuous] = limo_design_matrix(Y, Cat, ...
        %             Cont,directory, zscoring,full_factorial,flag);
        %
        X                  = contrast_mat_AC(trial_counts, condition_order);
        nb_conditions      = 3; %change the number of conditions that have to be estimated
        nb_interactions    = 0;
        nb_continuous      = 0;
        
        % we go for each subject and each channel seperately - tf and time are
        % concatenated so the data Y has the shape trials x cat(time, freq). X
        % will be a two dimensional design-matrix
        
        % reshape dimension of our data
        forced_dim      = [size(data,1), size(data,2)*size(data,3), size(data,4)]; % collapse frequencies and time dimensions
        Y_data               = limo_tf_4d_reshape(data,forced_dim);
        
        X               = X; %2 dimensional design matrix
        nb_conditions   = nb_conditions; %a vector indicating the number of conditions per factor
        nb_interactions = nb_interactions; %a vector indicating number of columns per interactions
        nb_continuous   = nb_continuous; %number of covariates
        methods          = 'OLS'; %weighted least squares
        n_freqs         = size(data,2); % number of frequency bins
        n_times         = size(data,3); %number of time bins
        analysis_type   = 'Time-Frequency'; % type of analysis
        
        Betas = zeros(size(Y_data,1), size(Y_data, 2),size(X,2));
        %Yhat  = zeros(size(Y_data));
        %Res   = zeros(size(Y_data));
        R2    = zeros(size(Y_data,1),size(Y_data, 2),3);
        
        % now perform glm electrode by electrode
        for elec = 1:size(Y_data,1)
            fprintf('Performing GLM on electrode %s of participant %s for effect %s. \n',num2str(elec), Participant_IDs{sub}, num2str(effect ))
            Y = squeeze(Y_data(elec, :, :));
            model = limo_glm(Y', X, nb_conditions, nb_interactions, ...
                nb_continuous, methods, analysis_type, n_freqs, n_times);
            Betas(elec,:,:,1)       = model.betas';
            %fitted_data             = X*model.betas;
            %Yhat(elec, :, :)        = fitted_data';
            %Res(elec, :, :)         = Y - fitted_data';
            R2(elec,:,1)       = model.R2_univariate;
            R2(elec,:,2)       = model.F;
            R2(elec,:,3)       = model.p;
        end
        
        % so limo_tf_4d has concatenated the frequencies for each time point,
        % i.e. 20frequencies for timepoint one, followed by 20 frequencies by
        % timepoint two and so forth. We have the same structure in the beta
        % weights and all other, which we now want to transform back into a
        % matrix that contains frequencies.
        
        
        % save data
        save_loc = fullfile(dirs.home, Participant_IDs{sub},limo_dir, save_n);
        save(save_loc,'Betas', 'R2', '-v7.3') %'Yhat', 'Res', 'R2', '-v7.3')
        clear data Y_data forced_dim condition_order trial_counts
        
    end
end
