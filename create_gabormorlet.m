% create_gabormorelet.m
%
% This function generates a series of Gabor/Morlet filters
%
% Fs = sampling frequency in Hz (e.g. 1000)
% Flo = lowest frequency Gabor in Hz (e.g. 5)
% Fhi = highest frequency Gabor in Hz (e.g. 60)
% Nsteps = number of frequency steps from Flo to Fhi
% Bandwidth = Gabor bandwidth in octaves (e.g. 1/4)
% M = parameter that controls tradeoff in resolution
% between time and frequency (e.g. M = 12)

function [freq gabor] = create_gabormorlet(Fs,Flo,Fhi,Nsteps,Bandwidth)

% convert Gabor Bandwidth to "M" parameter used by Tallon-Baudry
K = sqrt(2*log(2));
M = K * (2^Bandwidth + 1) / (2^Bandwidth - 1);
B = log2((1+K/M)/(1-K/M));

% Nsteps_recommended = 2 * log2(Fhi/Flo) / Bandwidth

freq = exp(linspace(log(Flo),log(Fhi),Nsteps));
phase = pi/4;
gabormorlet = [];
for i = 1:Nsteps
    f = freq(i);
    sigma_freq = f/M;
    sigma_time = 1/(2*pi*sigma_freq);
    w = ceil(3.5*sigma_time*Fs);
    t = (-w:w)*(1/Fs);
    env = exp(-t.^2/(2*sigma_time^2));
    wave = cos(2*pi*f*t+phase);
    g = env .* wave;
    h = fliplr(g);
    amp = sqrt(sum(g.^2+h.^2));
    g = (1/amp) * g;
    h = (1/amp) * h;
    figure(1), plot(t,g,'ro-', t,h,'go-'), axis tight
    gabor{i}.g = g;
    gabor{i}.h = h;
end

