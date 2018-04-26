function plot_one_spectrogram(power_spectral_density, times, freqs, freqs_on_y_axis, plotName)

figure;
h=surf(times, freqs, 10*log10(abs(power_spectral_density)), 'EdgeColor', 'none');
set(get(h,'Parent'),'YScale','log');
axis xy; axis tight; colormap(jet); view(0,90);
if freqs_on_y_axis~=0
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(plotName);
else
    xlabel('Frequency (Hz)');
    ylabel('Time (s)');
    title(plotName);
end
