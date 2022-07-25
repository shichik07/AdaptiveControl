function X = contrast_mat_AC(conditions, condition_order)
% transforms our conditions into a contrast matrix with contrasts for
% congruency, block_type(list-wise or item-specific) and the interaction between both).

X           = zeros(sum(conditions), 4); %initialize matrix
X(:,4)      = 1; %intercept

for con = 1:length(condition_order)
    % first we get the index of interest
    if con == 1
        index = 1:conditions(con);
    else
        start = sum(conditions(1:con-1));
        index = start + 1 : start + conditions(con);
    end
    
    % next we determine the index per condition - the first columns for the
    % congruency effect, the second columns for the main effect for the
    % main effect for block and the third column for the interaction
    % effect.
    
    
    if condition_order{con} == 'MI_I'
        X(index, 1) = -0.5; % main congruency
        X(index, 2) = -0.5; % main block 
        X(index, 3) =  0.5; % interaction
    elseif condition_order{con} == 'MI_C'
        X(index, 1) =  0.5;
        X(index, 2) = -0.5;
        X(index, 3) = -0.5;
    elseif condition_order{con} == 'MC_I'
        X(index, 1) = -0.5;
        X(index, 2) =  0.5;
        X(index, 3) = -0.5;
    elseif condition_order{con} == 'MC_C'
        X(index, 1) =  0.5;
        X(index, 2) =  0.5;
        X(index, 3) =  0.5;
    end
end
end