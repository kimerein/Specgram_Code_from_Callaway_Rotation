function [freqs, ti]=matlab_spectrogram_config

% freqs returns the frequencies at which to evaluate the spectrogram
freqs=zeros(1,3);

% freqs(1) is the lowest frequency for the spectrogram
freqs(1)=1;

% freqs(2) is the frequency step to use
freqs(2)=0.1;

% freqs(3) is the highest frequency for the spectrogram
freqs(3)=100;

% ti is the title for the spectrogram
ti=['Spectrogram - Frequency Resolution=' num2str(freqs(2), '%e')];