function [gcis, degg]= getGCIs(egg, fs)
% Detects Glottal Closure Instants (GCIs) from an EGG signal using using
% its derivative. Performs wavelet denoising using MATLAB wdenoise function with the default parameters,
% followed by a 30-tap FIR bandpass differentiator (20–7000 Hz).
% Inputs:
%   egg     - [1xN]  Electroglottography signal (EGG) 
%   fs      - Sampling frequency 
% Outputs:
%   gcis     - [1×M] Indices of detected Glottal Closure Instants (GCIs)
%   degg     - [1×N] Differentiated EGG signal (dEGG)   

    % Denoise the signal without altering its shape
    eggDEN = wdenoise(egg);

    % Obtain the dEGG using a bandpass differentiator to avoid amplifying
    % high-frequency noise in the resulting signal
    Nf = 30;
    Fpass = 20;
    Fstop = 7000;
    d = designfilt('differentiatorfir','FilterOrder',Nf, ...
        'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
        'SampleRate',fs);
    dt = 1/fs;
    degg  = filter(d, eggDEN)/dt;
    
    % Adjust delay to correct GCI timing.
    % shift has no impact since edges are unvoiced.
    delay = mean(grpdelay(d));
    degg = circshift(degg, -delay); 

    % Normalize the dEGG signal to [-1,1]
    degg = 2*((degg - min(degg))./(max(degg)-min(degg))) - 1;
    
    % Define an initial threshold on the height based on the signal's maximum
    threshold   = 0.25*max(degg);
    [~,locs]  = findpeaks(degg,'MinPeakHeight',threshold);
    
    % Estimate the period and do a second pass to discard false positive
    % defining a second threshold on the distance
    T  = round(median(diff(locs)));
    [~,gcis] = findpeaks(degg, ...
                           'MinPeakHeight',threshold, ...
                           'MinPeakDistance',floor(0.8*T));
end