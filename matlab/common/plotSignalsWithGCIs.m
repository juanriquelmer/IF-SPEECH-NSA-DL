function plotSignalsWithGCIs(x, xLabelStr, y, yLabelStr, fs, gcis, showFull, tInit, tEnd)
    t = (0:length(x)-1) / fs;

    figure;
    ax1 = subplot(2,1,1);
    plot(t, x, 'LineWidth', 1.5);
    hold on;
    plot(t(gcis), x(gcis), 'rv', 'MarkerFaceColor', 'r');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title([xLabelStr, ' with Marked GCIs']);
    grid on;
    axis tight;
    if ~showFull
        xlim([tInit, tEnd]);
    end

    ax2 = subplot(2,1,2);
    plot(t, y, 'LineWidth', 1.5);
    hold on;
    plot(t(gcis), y(gcis), 'ko', 'MarkerFaceColor', 'k');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title([yLabelStr, ' with Marked GCIs']);
    grid on;
    axis tight;
    if ~showFull
        xlim([tInit, tEnd]);
    end

    linkaxes([ax1, ax2], 'x');
end
