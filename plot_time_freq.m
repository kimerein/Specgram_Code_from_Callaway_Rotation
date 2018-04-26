function plot_time_freq(data, TimeLo, TimeHi, TimeStep, FreqLo, FreqHi, FreqStep, yes_logY, Y_ticks)

[Nfreq Ntime] = size(data);

x = linspace(TimeLo,TimeHi,Ntime);
Xtick = TimeLo : TimeStep : TimeHi;

if yes_logY==1
    % Logarithmic frequency axis
    logy = linspace(log(FreqLo),log(FreqHi),Nfreq);
    Ytick = Y_ticks;
    imagesc(x,logy,data);
    set(gca, ...
        'XTick', Xtick, ...
        'XTickLabel', arrayfun(@num2str, Xtick, 'UniformOutput', false), ...
        'Ydir','normal', ...
        'YTick', log(Ytick), ...
        'YTickLabel', arrayfun(@num2str, Ytick, 'UniformOutput', false));
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
else
    % Linear frequency axis
    y = linspace(FreqLo,FreqHi,Nfreq);
    Ytick = FreqLo : FreqStep : FreqHi;
    imagesc(x,y,data);
    set(gca, ...
        'XTick', Xtick, ...
        'XTickLabel', arrayfun(@num2str, Xtick, 'UniformOutput', false), ...
        'Ydir','normal', ...
        'YTick', Ytick, ...
        'YTickLabel', arrayfun(@num2str, Ytick, 'UniformOutput', false));
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
end