% plot_spectrograms.m

% Plots spectrograms for data processed with analyze_data.m
% Data parsed into regimes: UP and DOWN states; before, during and after
% stimulation
% Data stored in .mat files from initial analysis

%% Load processed data
stim_triggered_data0=load('92291_all_stim_triggered_data.mat');
stim_triggered_data=stim_triggered_data0.stimTriggeredData;
spontaneous_UPs0=load('92291_all_spontaneous_UPs.mat');
spontaneous_UPs=spontaneous_UPs0.spontaneousUPs;
nonSpontaneous_UPs0=load('92291_all_nonSpontaneous_UPs.mat');
nonSpontaneous_UPs=nonSpontaneous_UPs0.nonSpontaneousUPs;
all_UPs0=load('92291_all_all_UPs.mat');
all_UPs=all_UPs0.allUPs;

%% Load LFP data
% Load .mat data imported from Spike2
LFP_data=load('C:\MATLAB7\work\Code for Callaway Lab Rotation\92291_LFP_all.mat');
%low_pass_data=LFP_data.lue_light_on_eyes_on_every_30_seconds_up_and_down_states_Ch2.values;
low_pass_data=LFP_data.LFP_all;
%lowPass_samplingRate=1/LFP_data.lue_light_on_eyes_on_every_30_seconds_up_and_down_states_Ch2.interval;
lowPass_samplingRate=1/3.4000e-004;
clear LFP_data

%% Data analysis - Plot spectrograms
% Option=1. Plot using Matlab's spectrogram function 
% OR Option=2. Plot using Gabor-Morlet wavelets (code modified from dxjones)
option=2;

% To show signal below spectogram, choose show_signal=1
% else choose show_signal=0
% Will only show signal for single trace
show_signal=0;

% Analyze a specific data segment, or average a particular regime over multiple trials
% Regimes:
% regime=1  specific time regime
% regime=2  time regime triggered to stimulus
% regime=3  spontaneous or evoked UP state not coinciding with stimulation
% regime=4  non-spontaneous UP state coinciding with stimulation
% regime=5  all UP states
% regime=6  DOWN state coinciding with stimulation
% regime=7  DOWN state not coinciding with stimulation
% regime=8  all DOWN states
% regime=9  UP state(s) before stimulus
% regime=10 DOWN state(s) before stimulus
% regime=11 UP state(s) after stimulus
% regime=12 DOWN state(s) after stimulus
regime=2;

switch regime
    case 1
        % Time to analyze is from t1 to t2
        t1=20; % in seconds
        t2=30; % in seconds
        sig=cell(1,1);
        sig{1}=low_pass_data(floor(t1*lowPass_samplingRate):floor(t2*lowPass_samplingRate));
    case 2
        % Analyze a time regime triggered by stimulus
        % Window from t_before seconds prior to the onset of stimulation
        % to t_after seconds after the end of stimulation
        t_before=0.2; % in seconds, can be negative
        t_after=1; % in seconds, can be negative
        % To average over all trials (all stim. events), choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='some';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:length(stim_triggered_data);
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        elseif strcmp(which_to_average, 'some')
            % Choose which trials to average
            % Trials are numbered according to stim. event, e.g.,
            % first trial is triggered to the first stim. event in trace
            average_these_trials=55:64; 
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if stim_triggered_data(i).stimTrainOnset-t_before<0 && 1+floor((stim_triggered_data(i).stimTrainEnd+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif stim_triggered_data(i).stimTrainOnset-t_before<0
                sig{j}=low_pass_data(1:1+floor((stim_triggered_data(i).stimTrainEnd+t_after)*lowPass_samplingRate));
            elseif 1+floor((stim_triggered_data(i).stimTrainEnd+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((stim_triggered_data(i).stimTrainOnset-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((stim_triggered_data(i).stimTrainOnset-t_before)*lowPass_samplingRate):1+floor((stim_triggered_data(i).stimTrainEnd+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end
    case 3
        % Analyze spontaneous or evoked UP state(s) not coinciding with
        % stimulation
        % Window from t_before seconds prior to the onset of UP
        % to t_after seconds after the end of UP
        t_before=0.2; % in seconds
        t_after=0.5; % in seconds
        % To average over all UPs not coinciding with stimulation, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='some';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(spontaneous_UPs,1);
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        elseif strcmp(which_to_average, 'some')
            % Choose which UPs to average
            % UPs are numbered from beginning of trace to end
            average_these_trials=1801:2000; 
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if spontaneous_UPs(i,1)-t_before<0 && 1+floor((spontaneous_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif spontaneous_UPs(i,1)-t_before<0
                sig{j}=low_pass_data(1:1+floor((spontaneous_UPs(i,2)+t_after)*lowPass_samplingRate));
            elseif 1+floor((spontaneous_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((spontaneous_UPs(i,1)-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((spontaneous_UPs(i,1)-t_before)*lowPass_samplingRate):1+floor((spontaneous_UPs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end    
    case 4
        % Analyze UP state(s) coinciding with stimulation
        % Window from t_before seconds prior to the onset of UP
        % to t_after seconds after the end of UP
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all UPs coinciding with stimulation, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(nonSpontaneous_UPs,1);
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        elseif strcmp(which_to_average, 'some')
            % Choose which UPs to average
            % UPs are numbered from beginning of trace to end
            average_these_trials=1:5; 
            disp(sprintf('Total Number of Trials to Average: %d\n', length(average_these_trials)));
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if nonSpontaneous_UPs(i,1)-t_before<0 && 1+floor((nonSpontaneous_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif nonSpontaneous_UPs(i,1)-t_before<0
                sig{j}=low_pass_data(1:1+floor((nonSpontaneous_UPs(i,2)+t_after)*lowPass_samplingRate));
            elseif 1+floor((nonSpontaneous_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((nonSpontaneous_UPs(i,1)-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((nonSpontaneous_UPs(i,1)-t_before)*lowPass_samplingRate):1+floor((nonSpontaneous_UPs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end    
    case 5
        % Average over some or all UP states
        % Window from t_before seconds prior to the onset of UP
        % to t_after seconds after the end of UP
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all UPs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(all_UPs,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which UPs to average
            % UPs are numbered from beginning of trace to end
            average_these_trials=[2 3 4 5]; 
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if all_UPs(i,1)-t_before<0 && 1+floor((all_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif all_UPs(i,1)-t_before<0
                sig{j}=low_pass_data(1:1+floor((all_UPs(i,2)+t_after)*lowPass_samplingRate));
            elseif 1+floor((all_UPs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((all_UPs(i,1)-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((all_UPs(i,1)-t_before)*lowPass_samplingRate):1+floor((all_UPs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end    
    case 6
        % Average over DOWN state(s) coinciding with stimulation
        % Find DOWN state(s) coinciding with stimulation
        all_stim_DOWNs=[];
        for i=1:length(stim_triggered_data)
            between_UPs=all_UPs(intersect(find(all_UPs(:,1)>=stim_triggered_data(i).stimTrainOnset),find(all_UPs(:,1)<stim_triggered_data(i).stimTrainEnd)),:);
            DOWN_ends=[between_UPs(:,1); stim_triggered_data(i).stimTrainEnd];
            DOWN_starts=[stim_triggered_data(i).stimTrainOnset; between_UPs(:,2)];
            DOWNs=[DOWN_starts DOWN_ends];
            all_stim_DOWNs=[all_stim_DOWNs; DOWNs];
        end
        % Window from t_before seconds prior to the onset of DOWN
        % to t_after seconds after the end of DOWN
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all DOWNs coinciding with stimulation, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(DOWNs,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which DOWNs to average
            % DOWNs are numbered from beginning of trace to end
            average_these_trials=[1 2]; 
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if all_stim_DOWNs(i,1)-t_before<0 && 1+floor((all_stim_DOWNs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif all_stim_DOWNs(i,1)-t_before<0
                sig{j}=low_pass_data(1:1+floor((all_stim_DOWNs(i,2)+t_after)*lowPass_samplingRate));
            elseif 1+floor((all_stim_DOWNs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((all_stim_DOWNs(i,1)-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((all_stim_DOWNs(i,1)-t_before)*lowPass_samplingRate):1+floor((all_stim_DOWNs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end
    case 7
        % Average over DOWN state(s) not coinciding with stimulation
        % Find DOWN state(s) not coinciding with stimulation
        all_DOWN_starts=[0; all_UPs(:,2)];
        all_DOWN_ends=[all_UPs(:,1); length(low_pass_data)/lowPass_samplingRate];
        stimTrainOnsets=zeros(length(stim_triggered_data),1);
        stimTrainEnds=zeros(length(stim_triggered_data),1);
        for i=1:length(stim_triggered_data)
            stimTrainOnsets(i)=stim_triggered_data(i).stimTrainOnset;
            stimTrainEnds(i)=stim_triggered_data(i).stimTrainEnd;
        end
        [spont_DOWNs, nonspont_DOWNs]=sort_DOWN_states(stimTrainOnsets, stimTrainEnds, all_DOWN_starts, all_DOWN_ends, 0:(1/lowPass_samplingRate):length(low_pass_data)/lowPass_samplingRate, 0);
        % Window from t_before seconds prior to the onset of DOWN
        % to t_after seconds after the end of DOWN
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all DOWNs not coinciding stimulation, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='some';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(spont_DOWNs,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which DOWNs to average
            % DOWNs are numbered from beginning of trace to end
            average_these_trials=1:200; 
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if spont_DOWNs(i,1)-t_before<0 && floor((spont_DOWNs(i,2)+t_after)*lowPass_samplingRate)+1>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif spont_DOWNs(i,1)-t_before<0
                sig{j}=low_pass_data(1:floor((spont_DOWNs(i,2)+t_after)*lowPass_samplingRate)+1);
            elseif floor((spont_DOWNs(i,2)+t_after)*lowPass_samplingRate)+1>length(low_pass_data)
                sig{j}=low_pass_data(floor((spont_DOWNs(i,1)-t_before)*lowPass_samplingRate)+1:end);
            else
                sig{j}=low_pass_data(floor((spont_DOWNs(i,1)-t_before)*lowPass_samplingRate)+1:1+floor((spont_DOWNs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end
    case 8
        all_DOWN_starts=[0; all_UPs(:,2)];
        all_DOWN_ends=[all_UPs(:,1) length(low_pass_data)/lowPass_samplingRate];
        all_DOWNs=[all_DOWN_starts all_DOWN_ends];
        % Window from t_before seconds prior to the onset of DOWN
        % to t_after seconds after the end of DOWN
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all DOWNs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(all_DOWNs,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which DOWNs to average
            % DOWNs are numbered from beginning of trace to end
            average_these_trials=[1 2]; 
        end
        sig=cell(1,length(average_these_trials));
        j=1;
        for i=average_these_trials
            if all_DOWNs(i,1)-t_before<0 && 1+floor((all_DOWNs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1:end);
            elseif all_DOWNs(i,1)-t_before<0
                sig{j}=low_pass_data(1:1+floor((all_DOWNs+t_after)*lowPass_samplingRate));
            elseif 1+floor((all_DOWNs(i,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                sig{j}=low_pass_data(1+floor((all_DOWNs(i,1)-t_before)*lowPass_samplingRate):end);
            else
                sig{j}=low_pass_data(1+floor((all_DOWNs(i,1)-t_before)*lowPass_samplingRate):1+floor((all_DOWNs(i,2)+t_after)*lowPass_samplingRate));
            end
            j=j+1;
        end
    case 9
        % Window from t_before seconds prior to the onset of UP
        % to t_after seconds after the end of UP
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all UPs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(stim_triggered_average,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which stimulus trials to consider
            average_these_trials=[1 2]; 
        end
        % Consider n UPs before each stimulus
        n=1;
        sig=cell(1,length(average_these_trials)*n);
        j=1;
        for i=average_these_trials
            for k=1:n
                if stim_triggered_data(i).UPs_before_stim(k,1)-t_before<0 && 1+floor((stim_triggered_data(i).UPs_before_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1:end);
                elseif stim_triggered_data(i).UPs_before_stim(k,1)-t_before<0
                    sig{j}=low_pass_data(1:1+floor((stim_triggered_data(i).UPs_before_stim(k,2)+t_after)*lowPass_samplingRate));
                elseif 1+floor((stim_triggered_data(i).UPs_before_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).UPs_before_stim(k,1)-t_before)*lowPass_samplingRate):end);
                else
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).UPs_before_stim(k,1)-t_before)*lowPass_samplingRate):1+floor((stim_triggered_data(i).UPs_before_stim(k,2)+t_after)*lowPass_samplingRate));
                end
                j=j+1;
            end
        end
    case 10
        % Window from t_before seconds prior to the onset of DOWN
        % to t_after seconds after the end of DOWN
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all DOWNs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(stim_triggered_data,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which stimulus trials to consider
            average_these_trials=[1 2]; 
        end
        % Consider n UPs before each stimulus
        n=1;
        sig=cell(1,length(average_these_trials)*n);
        j=1;
        for i=average_these_trials
            for k=1:n
                if stim_triggered_data(i).DOWNs_before_stim(k,1)-t_before<0 && 1+floor((stim_triggered_data(i).DOWNs_before_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1:end);
                elseif stim_triggered_data(i).DOWNs_before_stim(k,1)-t_before<0
                    sig{j}=low_pass_data(1:1+floor((stim_triggered_data(i).DOWNs_before_stim(k,2)+t_after)*lowPass_samplingRate));
                elseif 1+floor((stim_triggered_data(i).DOWNs_before_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).DOWNs_before_stim(k,1)-t_before)*lowPass_samplingRate):end);
                else
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).DOWNs_before_stim(k,1)-t_before)*lowPass_samplingRate):1+floor((stim_triggered_data(i).DOWNs_before_stim(k,2)+t_after)*lowPass_samplingRate));
                end
                j=j+1;
            end
        end
    case 11
        % Window from t_before seconds prior to the onset of UP
        % to t_after seconds after the end of UP
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all UPs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(stim_triggered_data,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which stimulus trials to consider
            average_these_trials=1:5; 
        end
        % Consider nth UPs after each stimulus
        numUPs_perStim=1;
        n=2;
        %n=[1 2];
        sig=cell(1,length(average_these_trials)*numUPs_perStim);
        j=1;
        for i=average_these_trials
            for k=n
                k
                if stim_triggered_data(i).UPs_after_stim(k,1)-t_before<0 && 1+floor((stim_triggered_data(i).UPs_after_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1:end);
                elseif stim_triggered_data(i).UPs_after_stim(k,1)-t_before<0
                    sig{j}=low_pass_data(1:1+floor((stim_triggered_data(i).UPs_after_stim(k,2)+t_after)*lowPass_samplingRate));
                elseif 1+floor((stim_triggered_data(i).UPs_after_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).UPs_after_stim(k,1)-t_before)*lowPass_samplingRate):end);
                else
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).UPs_after_stim(k,1)-t_before)*lowPass_samplingRate):1+floor((stim_triggered_data(i).UPs_after_stim(k,2)+t_after)*lowPass_samplingRate));
                end
                j=j+1;
            end
        end
    case 12
        % Window from t_before seconds prior to the onset of DOWN
        % to t_after seconds after the end of DOWN
        t_before=0; % in seconds
        t_after=0; % in seconds
        % To average over all DOWNs, choose
        % which_to_average='all'; else which_to_average='some'
        which_to_average='all';
        if strcmp(which_to_average, 'all')
            average_these_trials=1:size(stim_triggered_data,1);
        elseif strcmp(which_to_average, 'some')
            % Choose which stimulus trials to consider
            average_these_trials=[1 2]; 
        end
        % Consider n UPs before each stimulus
        n=1;
        sig=cell(1,length(average_these_trials)*n);
        j=1;
        for i=average_these_trials
            for k=1:n
                if stim_triggered_data(i).DOWNs_after_stim(k,1)-t_before<0 && 1+floor((stim_triggered_data(i).DOWNs_after_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1:end);
                elseif stim_triggered_data(i).DOWNs_after_stim(k,1)-t_before<0
                    sig{j}=low_pass_data(1:1+floor((stim_triggered_data(i).DOWNs_after_stim(k,2)+t_after)*lowPass_samplingRate));
                elseif 1+floor((stim_triggered_data(i).DOWNs_after_stim(k,2)+t_after)*lowPass_samplingRate)>length(low_pass_data)
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).DOWNs_after_stim(k,1)-t_before)*lowPass_samplingRate):end);
                else
                    sig{j}=low_pass_data(1+floor((stim_triggered_data(i).DOWNs_after_stim(k,1)-t_before)*lowPass_samplingRate):1+floor((stim_triggered_data(i).DOWNs_after_stim(k,2)+t_after)*lowPass_samplingRate));
                end
                j=j+1;
            end
        end
end

% Option=1. Plot using Matlab's spectrogram function 
% OR Option=2. Plot using Gabor-Morlet wavelets (code modified from
% dxjones)

% Cut all elements in sig to match the length of sig{1}
% Spectrograms must have same sizes to average
max_sig_size=0;
size(sig{1});
for i=1:length(sig)
   if length(sig{i})>max_sig_size
       max_sig_size=length(sig{i});
   end
end
%for i=2:length(sig)
%    this_sig=sig{i};
%    sig{i}=this_sig(1:length(sig{1}));
%end
for i=1:length(sig)
   this_sig=sig{i};
   if length(sig{i})<max_sig_size
        sig{i}=[this_sig; zeros(max_sig_size-length(this_sig),1)];
    end
end

%min_sig_size=length(sig{1});
min_sig_size=1.7*lowPass_samplingRate;
for i=1:length(sig)
  this_sig=sig{i};
  if length(sig{i})>min_sig_size
      sig{i}=this_sig(1:min_sig_size);
  end
end

if option==1
    % Get frequencies at which to evaluate spectrogram
    % And title for the plot
    [freqs, title]=matlab_spectrogram_config;
	F = freqs(1):freqs(2):freqs(3);
    sumP=[];
    for i=1:length(sig)
        sprintf('Calculating spectrogram for sample: %d\n',i)
        [Y,F,T,P]=spectrogram(sig{i},256,250,F,lowPass_samplingRate,'yaxis');
        if k==1
            sumP=P;
        else
            sumP=sumP+P;
        end
    end
    avP=sumP/length(sig);
	plot_one_spectrogram(avP, T, F, 1, title);
elseif option==2
    % X_div gives the time between tick marks for the spectrogram
    X_div=0.1; % in seconds
    % Y_div gives the frequency bands between tick marks for the
    % spectrogram
    Y_div=10; % in Hz
    if exist('gbwavelets_freq.mat') && exist('gbwavelets_gabor.mat')
        [gb_params, freq, gabor]=gabor_morlet_config(lowPass_samplingRate,'gbwavelets_freq.mat','gbwavelets_gabor.mat');
    else
        [gb_params, freq, gabor]=gabor_morlet_config(lowPass_samplingRate,[],[]);
    end
    if show_signal && length(sig)==1
        gabor_morlet_plot(sig,gb_params,freq,gabor,1,sig{1},'LFP',X_div,Y_div,1,[1 10 20 30 40 50 60 80 100],lowPass_samplingRate);
    else
        gabor_morlet_plot(sig,gb_params,freq,gabor,0,[],'LFP',X_div,Y_div,1,[1 10 20 30 40 50 60 80 100],lowPass_samplingRate);
    end
end