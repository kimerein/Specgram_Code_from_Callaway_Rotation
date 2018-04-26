function regions=compare_freq_bands()
% Finds UP states according to the LFP power ratio

% Get LFP data
LFP_data=load('C:\MATLAB7\work\Code for Callaway Lab Rotation\92291_ch2_0to718s.mat');
low_pass_data=LFP_data.lue_light_on_eyes_on_every_30_seconds_up_and_down_states_Ch2.values;

lowPass_samplingRate=1/LFP_data.lue_light_on_eyes_on_every_30_seconds_up_and_down_states_Ch2.interval;
clear LFP_data

% Time to analyze is from t1 to t2
t1=20; % in seconds
t2=30; % in seconds
sig=cell(1,1);
sig{1}=low_pass_data(floor(t1*lowPass_samplingRate):floor(t2*lowPass_samplingRate));
 
% Plot using Gabor-Morlet wavelets (code modified from
% dxjones)
show_signal=1;
[gb_params, freq, gabor]=gabor_morlet_config(lowPass_samplingRate,'gbwavelets_freq.mat','gbwavelets_gabor.mat');
if show_signal && length(sig)==1
    LFP_specgram=gabor_morlet_plot(sig,gb_params,freq,gabor,1,sig{1},'LFP',1,10,1,[1 10 20 30 60 100],lowPass_samplingRate);
else
    LFP_specgram=gabor_morlet_plot(sig,gb_params,freq,gabor,0,[],'LFP',1,10,1,[1 10 20 30 60 100],lowPass_samplingRate);
end

% Show UP states according to following definition
% divide gives the distinction between "low" and "high" frequencies
band1_lower=30; % in Hz
band1_upper=50; % in Hz
band2_lower=50; % in Hz
band2_upper=80; % in Hz

% For each time point, plot integral (area under power spectrum)
% of frequencies greater than divide over the integral of frequencies 
% less than divide
% Use UP_thresh as threshold to identify UP states from this power ratio
% analysis
UP_thresh=5;

% Find row index in LFP_specgram that corresponds to divide
divide_ind1_lower=0;
divide_ind1_upper=0;
divide_ind2_lower=0;
divide_ind2_upper=0;
for i=1:length(freq)
    if freq(i)>=band1_lower
        divide_ind1_lower=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band1_upper
        divide_ind1_upper=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band2_lower
        divide_ind2_lower=i;
        break
    end
end
for i=1:length(freq)
    if freq(i)>=band2_upper
        divide_ind2_upper=i;
        break
    end
end
divide_ind1_lower
divide_ind1_upper
divide_ind2_lower
divide_ind2_upper

high_integral=zeros(size(LFP_specgram,2),1);
size_high_int=divide_ind1_upper-divide_ind1_lower+1;
size_low_int=divide_ind2_upper-divide_ind2_lower+1;
low_integral=zeros(size(LFP_specgram,2),1);
for i=1:size(LFP_specgram,2)
    high_area=sum(LFP_specgram(divide_ind1_lower:divide_ind1_upper,i));
    low_area=sum(LFP_specgram(divide_ind2_lower:divide_ind2_upper,i));
    high_integral(i)=high_area;
    low_integral(i)=low_area;
end

% Normalize high and low integral by number of frequencies included
power_ratio=(high_integral/size_high_int)./(low_integral/size_low_int);

figure;
intervals=(t2-t1)/(length(high_integral)-1);
times=t1:intervals:t2;
plot(times,power_ratio);
ti=sprintf('Power Ratio (Power In %d to %d Hz Band over Power in %d to %d Hz Band) from %.5f Secs to %.5f\nUP Threshold in Red',band1_lower,band1_upper,band2_lower,band2_lower,t1,t2);
title(ti);
xlim([t1 t2]);
xlabel('Time (s)');
ylabel('Power Ratio');
% Draw UP threshold
line([times(1) times(end)], [UP_thresh UP_thresh], 'Color', 'r');

% Return UPs
in_UP=0;
UP_starts=[];
UP_ends=[];
for i=1:size(LFP_specgram,2)    
    if in_UP==0
        if power_ratio(i)>=UP_thresh
            in_UP=1;
            UP_starts=[UP_starts; times(i)];
        end
    else
        if power_ratio(i)<UP_thresh
            in_UP=0;
            UP_ends=[UP_ends; times(i)];
        end
    end
end
if in_UP==1
    UP_ends=[UP_ends; times(end)];
end
regions=[UP_starts UP_ends];
        