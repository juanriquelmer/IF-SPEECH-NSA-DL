function plotSpectrogram(x, fs)
    t = (0:length(x)-1) / fs;
    window_length = round(0.025 * fs);   % 25 ms
    overlap = round(0.015 * fs);         % 15 ms (10 ms stride)
    nfft = 2^nextpow2(window_length);
    window = hamming(window_length, 'periodic');
    
    figure;
    subplot(2,1,1);
    plot(t, x);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original Signal');
    grid on;

    subplot(2,1,2);
    spectrogram(x, window, overlap, nfft, fs, 'yaxis');
    title('Spectrogram of EGG');
    ylabel('Frequency (kHz)');
    xlabel('Time (s)');
    colormap(gca, 'hot');
    ylim([0 3]);
end