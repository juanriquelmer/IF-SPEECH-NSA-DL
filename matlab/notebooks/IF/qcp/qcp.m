function [timeMarks] = computeTimeMarks( ...
    speech, fs, ...
    F0min, F0max, F0MarginFactor, F0_method, GCImethod, ...
    selectedPercent, ...
    WindType, WinLength, WinShift, phi, ...
    fileF0 ...
    )
% computeTimeMarks Compute time marks for speech analysis 

% Inputs:
%   speech          - The speech signal as a vector of samples
%   fs              - Sampling frequency in Hz
%
%   F0min, F0max    - Lower/upper pitch bounds in Hz
%   F0MarginFactor  - Number of semitones to expand around fileF0
%   F0_method       - Pitch method (e.g., 'NCF', 'PEF', 'CEP', 'LHS', 'SRH')
%   GCImethod       - GCI detection method (e.g., 'SEDREAMS' or 'GLOAT')
%
%   selectedPercent - Two-element vector indicating percent range of signal
%                     to consider, e.g. [5 95]
%
%   WindType        - Window mode: 'Constant Frame Rate' or 'Pitch Synchronous'
%   WinLength       - Frame length in seconds (if WindType = 'Constant Frame Rate')
%   WinShift        - Frame shift in seconds (if WindType = 'Constant Frame Rate')
%   phi             - Phase in % of pitch period for pitch marks 
%                     (only used in pitch-synchronous analysis)
%
%   fileF0          - The “ground truth” pitch stored in the file name 
%                     (used if you want to center your pitch range around this value)
%
% Output:
%   timeMarks       - struct with the following fields:
%       .analysisFramesSamples  - 2×N matrix, each column is [startIdx; endIdx]
%       .pitchMarks             - 2×(numPeriods) matrix, pitch-synchronous boundaries
%       .selectedSamplesIdxs    - The sample indices that lie in the selectedPercent range
%       .gcisTimes              - GCIs in seconds
%       .gcisSamples            - GCIs in samples
%       .isGCISelectedArray     - Boolean array indicating which GCIs fall inside selectedPercent
%       .numGCISelected         - How many GCIs are in the selectedPercent portion
%       .f0Min, .f0Max          - Actual F0 range used
%       .F0value                - Median pitch (Hz) of the entire signal
%       .F0Pulse                - Instantaneous pitch from consecutive GCIs
%       .semitonesDifferenceGT  - Semitone difference from median pitch to fileF0
%       .semitonesDifferencePulse - Per-pulse semitone difference from median pitch
%       .semitonesDifferencePulseGT - Per-pulse semitone difference from fileF0

%% 1) Handle F0 range
f0ValInFile = str2double(fileF0);

if ~isnan(F0MarginFactor) && ~isnan(f0ValInFile)
    %F0MarginFactor can be used to recenter the range around an expected pitch (extracted from the file info) 
    % by shifting it ± some number of semitones.
    f0MinFinal = f0ValInFile * 2 ^ (-F0MarginFactor / 12);
    f0MaxFinal = f0ValInFile * 2 ^ ( F0MarginFactor / 12);
    % F0 * 2^(± F0MarginFactor / 12) will scale F0 by some factor in semitones.
else
    % Fall back to direct values if something is missing
    f0MinFinal = F0min;
    f0MaxFinal = F0max;
end

timeMarks.f0Min = f0MinFinal;
timeMarks.f0Max = f0MaxFinal;


%% 2) Pitch computation (using MATLAB's built-in pitch function)

f0Contour = pitch( ...
    speech, fs, ...
    'Range', [timeMarks.f0Min timeMarks.f0Max], ...
    'Method', F0_method ...
);
% A single representative pitch for the entire signal.
timeMarks.F0value = median(f0Contour);


%% 3) GCIs computation

% ORIGINAL
% switch GCImethod
%     case 'SEDREAMS'
%         timeMarks.gcisTimes = gci_sedreams( ...
%             speech, fs, timeMarks.F0value, 1 ...
%         );
%         timeMarks.gcisSamples = round(timeMarks.gcisTimes * fs) + 1;
% 
%     case 'GLOAT' %(Drugman, T and Dutoit, T: "Glottal Closure andOpening Instant Detection from Speech Signals", 2009)
%         [timeMarks.gcisSamples, goi] = gci( ...
%             speech, timeMarks.F0value, fs ...
%         );
%         timeMarks.gcisTimes = timeMarks.gcisSamples / fs;
% 
%     otherwise
%         error(['Unknown GCImethod: ' GCImethod]);
% end

% From covarep/glottalsource/gci_sedreams
% sd_gci = gci_sedreams(x,fs,median(srh_f0),1);        % SEDREAMS
medianF0 = timeMarks.F0value(1,1);
 %speech, fs, timeMarks.F0value, 1 ...
timeMarks.gcisTimes = gci_sedreams( ...
            speech, fs, medianF0, 1 ...
         );
timeMarks.gcisSamples = round(timeMarks.gcisTimes * fs) + 1;

%% 4) Derive per-GCI F0 and semitone differences

% Once we have a series of GCIs in time, 
% we can derive the instantaneous pitch period by subtracting consecutive GCIs.
periodPulses = [0 timeMarks.gcisTimes];
timeMarks.F0Pulse = 1 ./ (periodPulses(2:end) - periodPulses(1:end-1));
% F0Pulse is effectively the pitch for each inter-GCI interval.

% If fileF0 is numeric, compare with your median F0
numericFileF0 = str2double(fileF0);
if isnan(numericFileF0)
    numericFileF0 = timeMarks.F0value; % fallback if parsing fails
end

% semitonesDifferenceGT: difference between the global (median) 
% pitch and the “ground truth” pitch from the file info.
timeMarks.semitonesDifferenceGT = 12 * log2(timeMarks.F0value / numericFileF0);

% semitonesDifferencePulse: difference between each local pitch (per GCI) and the median pitch.
timeMarks.semitonesDifferencePulse = 12 * log2(timeMarks.F0Pulse / timeMarks.F0value);

% semitonesDifferencePulseGT: difference between each local pitch
% and the ground-truth pitch from the file name.
timeMarks.semitonesDifferencePulseGT = 12 * log2(timeMarks.F0Pulse / numericFileF0);


%% 5) Pitch marks (pairs that form pitch-synchronous intervals)

% The variable pitchMarks is commonly used in pitch-synchronous analysis 
% to define the boundaries of small segments (pitch-synchronous frames).

% phi is some phase offset (in percent of the pitch period), 
% controlling exactly where within the GCI–GCI interval you place that mark.

% pitchMarks here is stored as pairs: each column is [start; end] of a pitch period region.

% The first boundary is at sample 1, the last boundary is at the very end of the speech signal. 
% The intermediate boundaries are derived from the GCIs plus that offset factor.

% The pitch mark is usually at a fraction 'phi' of the GCI–GCI interval
timeMarks.pitchMarks = [ ...
    1, ...
    round(timeMarks.gcisSamples(1:end-1) ...
          + (phi / 100) .* (timeMarks.gcisSamples(2:end) ...
                            - timeMarks.gcisSamples(1:end-1)) ), ...
    length(speech) ...
];
% Convert the 1D vector to pairs [start; end]
timeMarks.pitchMarks = [ ...
    timeMarks.pitchMarks(1:end-1); ...
    timeMarks.pitchMarks(2:end) ...
];


%% 6) Selected Samples

% selectedSamplesIdxs is the range of sample indices that lie within the chosen percentage cut.

if isempty(selectedPercent) || numel(selectedPercent) < 2
    selectedPercent = [0 100];
end
lowerIdx = round((selectedPercent(1) / 100) * length(speech)) + 1;
upperIdx = round((selectedPercent(2) / 100) * length(speech));
timeMarks.selectedSamplesIdxs = lowerIdx : upperIdx;

% Mark which GCIs are within that range
lowerSec = (lowerIdx - 1) / fs;
upperSec = (upperIdx - 1) / fs;

% isGCISelectedArray marks which GCIs lie within that time range.
timeMarks.isGCISelectedArray = ...
    (timeMarks.gcisTimes > lowerSec) & (timeMarks.gcisTimes < upperSec);

% numGCISelected is how many GCIs ended up in that subset of the signal.
timeMarks.numGCISelected = ...
    sum(timeMarks.isGCISelectedArray);


%% 7) Define analysis frames

% we either split the speech into fixed-length frames (e.g., 50 ms length, 25 ms shift) 
% or pitch-synchronous frames (where each frame spans a certain number of glottal pulses).

switch WindType
    case 'Constant Frame Rate'
        wind = round(WinLength * fs);
        step = round(WinShift * fs);

        % Using length(speech), not length(speech*fs)
        % speech is already in samples
        nSamples = length(speech);
        np = ceil((nSamples - wind) / step);

        starts = 1 : step : (step*np + 1);
        ends   = wind : step : (wind + step*np + 1);

        % Make sure last end = length(speech)
        ends(end) = nSamples;

        timeMarks.analysisFramesSamples = [starts; ends];

    case 'Pitch Synchronous'
        % Each frame includes 2 pulses
        % That means from pitchMark(i) to pitchMark(i+2)
        if size(timeMarks.pitchMarks, 2) < 3
            % If we have fewer than 3 pitch marks, we can't form a "2 pulses" frame 
            timeMarks.analysisFramesSamples = [];
        else
            timeMarks.analysisFramesSamples = [ ...
                timeMarks.pitchMarks(1, 1:end-2); ...
                timeMarks.pitchMarks(2, 3:end) ...
            ];
        end

    otherwise
        error(['Unknown WindType: ' WindType]);
end

end