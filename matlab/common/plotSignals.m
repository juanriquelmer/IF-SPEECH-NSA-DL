function plotSignals(x, xLabelStr, y, yLabelStr, fs, showFull, tInit, tEnd)
    t = (0:length(x)-1) / fs;

    figure;
    ax1 = subplot(2,1,1);
    plot(t, x, 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(xLabelStr);
    grid on;
    axis tight;
    if ~showFull
        xlim([tInit, tEnd]);
    end

    ax2 = subplot(2,1,2);
    plot(t, y, 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(yLabelStr);
    grid on;
    axis tight;
    if ~showFull
        xlim([tInit, tEnd]);
    end

    linkaxes([ax1, ax2], 'x');
end
