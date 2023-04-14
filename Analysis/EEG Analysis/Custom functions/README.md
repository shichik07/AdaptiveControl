Custom functions for the time frequency analysis of the EEG data. Attention - some functions are not my intellectual work but functions from the fieldtrip and LIMO toolbox, merely adapted so that we could use their functionality for our data

contrast_mat_AC.m - function to construct a contrast matrix for our data used for the regression analysis
eeglab2ft_like.m - function to transform eeglab data to fieldtrip data. Not sure if actually used or left-over
find_cluster.m - function to extract electrodes and time indices of the significant cluster identified via the fieldtrip cluster permutation output
get_pseudo_FT_struct.m - function to transform raw EEG data into a struct that fieldtrip can work with to perform cluster permutation testing
get_pseudo_FT_struct_cluster.m - I assume similar as above? One may be redundant
get_trial_numbers.m - function to extract trial numbers - not relevant for the current scripts
get_trlindices.m - relevant for the wavelet convolution 
get_pseudo_FT_struct_cluster.m - LIMO function rewritten, so that I can use the Limo functions in question without my data having to be in the LIMO format
theta_average_to_R.m - function to write eeg summary data so that I can import it in R
theta_mat_to_R.m - same, but in a matrix format 