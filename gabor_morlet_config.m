function [plot_params, freq, gabor]=gabor_morlet_config(freq_step, load_wavelets_freq, load_wavelets_gabor)

% plot_params contain the parameters for the spectrogram
plot_params=zeros(5,1);

% plot_params(1) is the frequency step 
% For maximal frequency resolution, use sampling rate of signal
plot_params(1)=freq_step;

% plot_params(2) is the lowest frequency to consider
plot_params(2)=1;

% plot_params(3) is the highest frequency to consider
plot_params(3)=100;

% plot_params(4) is the number of steps to use for the spectrogram
plot_params(4)=200;

% plot_params(5) is the bandwidth to use for the spectrogram
plot_params(5)=1/4;

% If load_wavelets_freq and load_wavelets_gabor are []
% need to create gabor-morlet wavelets,
% else load wavelets from these files
load_wavelets_freq
load_wavelets_gabor
if isempty(load_wavelets_freq) || isempty(load_wavelets_gabor)
    [freq gabor]=create_gabormorlet(plot_params(1),plot_params(2),plot_params(3),plot_params(4),plot_params(5));
    % Save freq and gabor in .mat files
    save 'gbwavelets_freq.mat' freq
    save 'gbwavelets_gabor.mat' gabor
else
    freq0=load(load_wavelets_freq);
    gabor0=load(load_wavelets_gabor);
    freq=freq0.freq;
    gabor=gabor0.gabor;
end
