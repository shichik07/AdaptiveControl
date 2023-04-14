Time Frequency Analysis of the EEG data

TF_00_1_Preprocessing_Correction.m - Preprocessing script up until ICA correction 

TF_Ana_00_EpochFreq.m - epoching and interpolation of prior rejected channels stimulus locked

TF_Ana_00b_EpochFreq_ResponseLocked.m - epoching and interpolation of prior rejected channels responses locked

TF_Ana_01_SurfaceLaplacian.m - Surface Laplacian Transformation of stimulus locked data

TF_Ana_01b_SurfaceLaplacian_ResponseL.m - Surface Laplacian Transformation of reponse locked data

TF_Ana_02_WaveletConvolution.m - Frequency Transformation via Morlet wavelet convolution of the stimulus locked data

TF_Ana_02b_WaveletConvolution_ResponseL.m - Frequency Transformation via Morlet wavelet convolution of the response locked data

TF_Ana_03_baseline_correction.m - Whole trial decibel conversion for the stimulus locked data

TF_Ana_03b_baseline_correctionResponseL.m - Whole trial decibel conversion for the response locked data

TF_Ana_04_LimoGLM.m - Linear model analysis of the stimulus locked data using functions from the LIMO toolbox

TF_Ana_04b_LimoGLM_ResponseL.m - Linear model analysis of the response locked data using functions from the LIMO toolbox

TF_Ana_05_FT_ClusterPermJoined.m - Cluster Permutation statistics of for both response locked and stimulus locked data to test proactive and reactive control effects in the theta band