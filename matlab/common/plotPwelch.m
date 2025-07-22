function plotPwelch(x, xLabelStr, fs, startFraction, endFraction)
    N = length(x);
    startIdx = floor(startFraction * N) + 1;
    endIdx = floor(endFraction * N);

    t = (0:N-1) / fs;
    selectedTime = t(startIdx:endIdx);

    figure;
    
    subplot(2,1,1)
    plot(t, x, 'b')
    hold on
    plot(selectedTime, x(startIdx:endIdx), 'r', 'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('Amplitude')
    title(xLabelStr)
    grid on

    subplot(2,1,2)
    pwelch(x(startIdx:endIdx), [], [], [], fs)
    title('Power Spectral Density')
end
