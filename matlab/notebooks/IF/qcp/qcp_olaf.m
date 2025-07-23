function [G, Gd, H_VTA, H_G] = qcp_olaf(s, fs, VTorder, GSorder, lipRad, DQ, PQ, RQ, causality, STcompens, tmarks, gcis)
% Quasi Closed-Phase glottal inverse filtering with OLA processing
%
% Inputs:
%   s         - [1×N] Speech signal (row vector)
%   fs        - Sampling frequency (Hz)
%   VTorder   - LPC order for vocal tract (e.g., 48)
%   GSorder   - LPC order for glottal source (e.g., 3)
%   lipRad    - Leaky integration coefficient (e.g., 0.99)
%   DQ        - Duration Quotient (0.4–1.0)
%   PQ        - Position Quotient (0–0.2)
%   RQ        - Ramp Quotient (e.g., 0.14)
%   causality - 'causal' or 'noncausal' for FIR filters
%   STcompens - Flag for spectral tilt compensation (true/false)
%   tmarks    - [2×Nf] frame boundaries (samples)
%   gcis      - Glottal closure instants (sample indices)
%
% Outputs:
%   G     - [1×N]   Reconstructed glottal flow signal
%   Gd    - [1×N]   Glottal flow derivative signal
%   H_VTA - {1×Nf}  Cell array; each element is VT LPC coefficients [1×VTorder]
%   H_G   - {1×Nf}  Cell array; each element is GS LPC coefficients [1×GSorder]

% Ensure column vector becomes row
if iscolumn(s)
    s = s';
end

% Initialize outputs and filter memories
N = length(s);
G  = zeros(1, N);
Gd = zeros(1, N);
nf = size(tmarks, 2);
H_VTA = cell(1, nf);
H_G   = cell(1, nf);

% Frame processing time marks
tcenter = round((tmarks(1, :) + tmarks(2, :)) / 2);
tmc = round((tcenter(2:end) + tcenter(1:end-1)) / 2);
tmc = [tmarks(1, 1), tmc, tmarks(2, end)];

% Inverse filter memories
Z_vt   = zeros(1, VTorder);
Z_pre  = 0;
Z_int  = 0;
Z_st   = 0;
Z_der  = 0;

% Window function for LPC
winfunc = 'hanning';

wmin = 1e-5;

Lpf = VTorder + 1;

for n = 1:nf
    % Extract frame for weighting
    x = s(tmarks(1, n):tmarks(2, n));
    pos_gci_frame = find((gcis >= tmarks(1, n)) & (gcis <= tmarks(2, n)));
    gci_ins = gcis(pos_gci_frame) - tmarks(1, n) + 1;
    T0 = mean(diff(gci_ins));
    Nramp = round(RQ * T0);
    w = makeW(x, VTorder, DQ, PQ, wmin, Nramp, gci_ins, fs);

    % Prepare signal segment for inverse filtering
    if n == 1
        segment = [linspace(-s(tmc(n)), s(tmc(n)), Lpf), s(tmc(n):(tmc(n+1)-1))];
        idx = (Lpf+1):length(segment);
    else
        segment = s(tmc(n):(tmc(n+1)-1));
        idx = 1:length(segment);
    end

    % Pre-emphasis
    [s2, Z_pre] = filter([1 -1], 1, x, Z_pre);

    % Vocal tract LPC on windowed, weighted frame
    sw = win(s2, winfunc);
    [Hvt, ~] = wlp(sw, w, VTorder);

    % Remove real poles if requested
    if false
        [z_p, p_p, k_p] = tf2zp(1, Hvt);
        p_p = p_p(p_p < 0 | abs(imag(p_p)) > 1e-15);
        [~, Hvt] = zp2tf(z_p, p_p, k_p);
    end

    % Inverse filtering to get glottal source estimate
    [dg, Z_vt] = filter(Hvt, 1, segment, Z_vt);
    dg = dg(idx);
    %[sg, Z_int] = integrate(dg, lipRad, causality, Z_int);
    [sg, ~] = integrate(dg, lipRad, causality, Z_int);

    % Glottal source LPC
    Hg = lpc(win(sg, winfunc), GSorder);

    if STcompens
        % Spectral tilt compensation
        %[bvtt, avtt] = invfreqz(freqz(1, Hvt), 1024, 0, 1);
        [bvtt, avtt] = invfreqz(freqz(1, Hvt), linspace(0, pi, 512), 0, 1);
        [sg, Z_st] = filter(bvtt, avtt, sg, Z_st);
        [dg, Z_der] = filter([1 -lipRad], 1, sg, Z_der);
        Hg = lpc(win(sg, winfunc), GSorder);
    end

    H_VTA{n} = Hvt(:)';
    H_G{n}   = Hg(:)';
    G(tmc(n):(tmc(n+1)-1))  = sg;
    Gd(tmc(n):(tmc(n+1)-1)) = dg;
end
end

function y = win(x, winfunc)
    y = feval(winfunc, length(x))' .* x;
end

function [y, Z] = integrate(x, rho, causality, Z)
    if strcmp(causality, 'causal')
        [y, Z] = filter(1, [1 -rho], x, Z);
    else
        [y, Z] = filter(1, [1 -rho], flip(x), Z);
        y = -flip(y);
    end
end

function p = get_if_order(fs)
    % Get suitable inverse filter order for given sampling frequency
    p = round(fs / 1000) + 2;
end

function w = makeW(x, p, DQ, PQ, d, Nramp, gci_ins, fs)
    % Create AME weight function for frame x
    N = length(x);
    if Nramp > 0
        UPramp = linspace(d, 1, 2 + Nramp);
        UPramp = UPramp(2:end-1);
        DOWNramp = UPramp(end:-1:1);
    end

    if DQ + PQ > 1
        DQ = 1 - PQ;
    end

    w = d * ones(1, N + p);
    if isempty(gci_ins)
        T2 = 0;
        T1 = 1;
    else
        for i = 1:(length(gci_ins) - 1)
            T = gci_ins(i+1) - gci_ins(i);
            T1 = round(DQ * T);
            T2 = round(PQ * T);
            while T1 + T2 > T
                T1 = T1 - 1;
            end
            w(gci_ins(i) + T2 : gci_ins(i) + T2 + T1 - 1) = 1;
            if Nramp > 0
                w(gci_ins(i) + T2 : gci_ins(i) + T2 + Nramp - 1) = UPramp;
                if gci_ins(i) + T2 + T1 - Nramp > 0
                    w(gci_ins(i) + T2 + T1 - Nramp : gci_ins(i) + T2 + T1 - 1) = DOWNramp;
                end
            end
        end
        T = gci_ins(end);
        T1 = round(DQ * T);
        T2 = round(PQ * T);
        Nend = N - (T2 + gci_ins(end));
        if T2 + gci_ins(end) < N
            if T1 + T2 < Nend
                w(gci_ins(end) + T2 : gci_ins(end) + T2 + T1 - 1) = 1;
                if Nramp > 0
                    w(gci_ins(end) + T2 : gci_ins(end) + T2 + Nramp - 1) = UPramp;
                    w(gci_ins(end) + T2 + T1 - Nramp : gci_ins(end) + T2 + T1 - 1) = DOWNramp;
                end
            else
                T1 = Nend - T2;
                w(gci_ins(end) + T2 : gci_ins(end) + T2 + T1 - 1) = 1;
                if Nramp > 0
                    w(gci_ins(end) + T2 : gci_ins(end) + T2 + Nramp - 1) = UPramp;
                end
            end
        end
    end
end

function [a, E] = wlp(x, w, p)
%WLP Perform Weighted Linear Prediction (WLP)
%
%   [a, s_hat, residuals, E] = weightedLinearPrediction(x, w, p)
%
%   Inputs:
%       x - Input signal vector (Nx1)
%       w - Window weights vector (Nx1), same length as x
%       p - Prediction order (integer)
%
%   Outputs:
%       a         - Linear prediction coefficients (px1)
%       s_hat     - Predicted signal values ((N-p)x1)
%       residuals - Prediction residuals ((N-p)x1)
%       E         - Weighted sum of squared residuals (scalar)

    % Ensure column vectors
    s = x(:);
    w = w(:);
    k = p;
    Nw = length(s);

    %if length(w) ~= Nw
    %    error('Input signal and window must be the same length.');
    %end

    if k >= Nw
        error('Prediction order p must be less than the length of the signal.');
    end

    % Valid indices for prediction
    validIdx = (k+1):Nw;
    L = length(validIdx);

    % Build predictor matrix S
    S = zeros(L, p);
    for row = 1:L
        n = validIdx(row);
        S(row, :) = s(n-1:-1:n-k).';
    end

    % Target vector
    y = s(validIdx);

    % Diagonal weight matrix
    Wmat = diag(w(validIdx));

    % Weighted least squares solution
    A = S' * Wmat * S;
    b = S' * Wmat * y;
    a = A \ b;

    % Predicted values and residuals
    s_hat = S * a;
    residuals = y - s_hat;

    % Weighted error
    E = residuals' * (Wmat * residuals);
    a = [1; -a]; % Added
end
