function plotFFT(x, fs)
    t = (0:length(x)-1) / fs;
    n = length(x);
    X = abs(fft(x));
    X_dB = 20*log10(X + eps);
    f = (0:n-1)*(fs/n);
    
    half_n = floor(n/2);
    f_pos = f(1:half_n) / 1000;
    X_dB_pos = X_dB(1:half_n);

    figure;
    subplot(2,1,1);
    plot(t, x);
    title('Time Domain');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2,1,2);
    plot(f_pos, X_dB_pos);
    title('Frequency Domain');
    xlabel('Frequency (kHz)');
    ylabel('Magnitude (dB)');
end
