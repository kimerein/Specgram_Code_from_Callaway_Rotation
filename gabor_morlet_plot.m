function LFP_specgram=gabor_morlet_plot(sig, gb_params, freq, gabor, do_plot_sig, sig_to_plot, sig_title, X_div, Y_div, yes_logY, Y_ticks, sampling_rate)

% Sig is an array of data segments to average

% gb_params give the gabor-morlet parameters used to calculate the
% spectrogram

% freq and gabor are variables that specify the gabor-morlet filters

% do_plot_sig equals 1 if want to plot signal beneath spectrogram, 
% but only works for non-averages

% sig_to_plot is the signal to plot beneath the spectrogram

% X_div gives the distance between tick marks for the time axis

% Y_div gives the distance between tick marks for the logarithmic frequency
% axis

% sampling_rate is the LFP signal sampling rate

% Fs is the frequency step to use for spectrogram
% Low and high frequencies for spectrogram are set when gabor-morlet
% filters are created (see gabor_morlet_config.m)
% Nsteps are the number of steps used to create the gabor-morlet wavelets
% Bandwidth is the bandwidth used for the spectrogram
% Bandwidth trades off resolution in the time and frequency domains
Fs=gb_params(1);
Nsteps=gb_params(4);
Bandwidth=gb_params(5);

sum_sigs=[];
for i=1:length(sig)
    sprintf('Calculating spectrogram for sample %d\n',i)
    signal=sig{i}';
    Nsamples = length(signal);
    tf  = gmfilterfast(signal,gabor);
    if i==1
        sum_sigs=tf;
    else
        % Cut all subsequent tf's to size of first tf
        if ~isequal(size(tf),size(sum_sigs))
            tf=tf(:,1:size(sum_sigs,2));
        end
        sum_sigs=sum_sigs+tf;
    end
end
av_sigs=sum_sigs/length(sig);
tf=av_sigs;

mtf = repmat(mean(tf,2),1,size(tf,2));
xtf = (tf-mtf) ./ (tf+mtf);

if do_plot_sig~=0 && length(sig)==1
    figure;
    subplot(2,1,1);
    plot_time_freq(tf, 0, Nsamples/sampling_rate, X_div, freq(1), freq(end), Y_div, yes_logY, Y_ticks);
    ti=sprintf('Time-Frequency Analysis of Data - Frequency Step=%e, Num. Steps for Calculation=%d, Bandwidth=%.4f',Fs,Nsteps,Bandwidth); 
    title(ti);
    subplot(2,1,2);
    plot(0:(Nsamples/Fs)/(length(sig_to_plot)-1):Nsamples/Fs, sig_to_plot);
	axis([0,Nsamples/Fs,min(sig_to_plot),max(sig_to_plot)]);
	xlabel('Time (sec)');
	ylabel(sig_title);
else
    figure;
    plot_time_freq(tf, 0, Nsamples/sampling_rate, X_div, freq(1), freq(end), Y_div, yes_logY, Y_ticks);
    ti=sprintf('Time-Frequency Analysis of Data - Frequency Step=%e, Num. Steps for Calculation=%d, Bandwidth=%.4f',Fs,Nsteps,Bandwidth); 
    title(ti);
end

if do_plot_sig~=0 && length(sig)==1
    figure;
    subplot(2,1,1);
    plot_time_freq(xtf, 0, Nsamples/sampling_rate,  X_div, freq(1), freq(end), Y_div, yes_logY, Y_ticks);
    ti=sprintf('Time-Frequency Analysis of Data (Normalized) - Frequency Step=%e, Num. Steps for Calculation=%d, Bandwidth=%.4f',Fs,Nsteps,Bandwidth); 
    title(ti);
    subplot(2,1,2);
    plot(0:(Nsamples/Fs)/(length(sig_to_plot)-1):Nsamples/Fs, sig_to_plot);
	axis([0,Nsamples/Fs,min(sig_to_plot),max(sig_to_plot)]);
	xlabel('Time (sec)');
	ylabel(sig_title);
else
    figure;
    plot_time_freq(xtf, 0, Nsamples/sampling_rate,  X_div, freq(1), freq(end), Y_div, yes_logY, Y_ticks);
    ti=sprintf('Time-Frequency Analysis of Data (Normalized) - Frequency Step=%e, Num. Steps for Calculation=%d, Bandwidth=%.4f',Fs,Nsteps,Bandwidth); 
    title(ti);
end
LFP_specgram=tf;