function reshaped = reshape_3d_limo_fixed(reshape_in, freqs, times)
% The Limo function requires a Limo object to reshape 3d to 4d data, but
% not the other way around. Hence we just copy the code here so that it
% works without the Limo object.
[n_elec, n_freq_times, N] = size(reshape_in);
n_freqs                   = size(freqs,2);
n_times                   = size(times, 2);
if n_freq_times ~= n_freqs*n_times
    error('dimensions disagreement to reshape freq*time')
else
    reshaped = nan(n_elec,n_freqs, n_times, N);
end

% For each trial, take the stack of (elec x (freqs timepoint))
% and split it into frames of (elec x freqs x times), then populate the 4D Y with this

for tr = 1:N
    eft_3d = NaN(n_elec,n_freqs,n_times);
    for tm = 1:n_times
        this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf
        eft_3d(:,:,tm) =reshape_in(:,this_freq_start_index:this_freq_start_index+n_freqs-1,tr);
    end
    reshaped(:,:,:,tr) = eft_3d;
end


end