function Items = get_trlindices(EEG)
% inputs EEGlab structre and returns struct with trial indices
Items.ISPC_MC_C =[];
Items.ISPC_MC_I =[];
Items.ISPC_MI_C =[];
Items.ISPC_MI_I =[];
Items.LWPC_MC_C =[];
Items.LWPC_MC_I =[];
Items.LWPC_MI_C =[];
Items.LWPC_MI_I =[];

%Correct_check = {}; 



for i=1:EEG.trials
    % get 0 trigger
    onset = cellfun(@(x)x == 0, EEG.epoch(i).eventlatency);
    
    if EEG.epoch(i).eventanswer{onset} == "correct" % proceed only with correctly answered trials
        Correct_check{end+1} = EEG.epoch(i).eventanswer{onset};
        if EEG.epoch(i).eventanalysetype{onset} == "MI"
            if EEG.epoch(i).eventcongruency{onset} == "incongruent"
                Items.LWPC_MI_I(end+1) = i;
            else
                Items.LWPC_MI_C(end+1) = i;
            end
        elseif EEG.epoch(i).eventanalysetype{onset} == "MC"
            if EEG.epoch(i).eventcongruency{onset} == "incongruent"
                Items.LWPC_MC_I(end+1) = i;
            else
                Items.LWPC_MC_C(end+1) = i;
            end
        elseif EEG.epoch(i).eventanalysetype{onset} == "main_incon"
            if EEG.epoch(i).eventcongruency{onset} == "incongruent"
                Items.ISPC_MI_I(end+1) = i;
            else
                Items.ISPC_MI_C(end+1) = i;
            end
        elseif EEG.epoch(i).eventanalysetype{onset} == "main_con"
            if EEG.epoch(i).eventcongruency{onset} == "incongruent"
                Items.ISPC_MC_I(end+1) = i;
            else
                Items.ISPC_MC_C(end+1) = i;
            end
        end
    end
end

%unique(cellfun(@unique, Correct_check)) % I believe at this poitn we only
%have correct trials - but just in case


end
  