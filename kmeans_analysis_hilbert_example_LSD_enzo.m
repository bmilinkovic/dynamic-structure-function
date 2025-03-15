%% Phase Coherence Analysis for LSD Study
% This script performs phase coherence analysis on LSD vs Placebo fMRI data
% using Hilbert transform and k-means clustering to identify brain states.
%
% Author: Borjan Milinkovic
% Last Updated: 2025
%
% Required Data Files:
% - LSD_TS_FC.mat: Contains time series data for LSD and Placebo conditions
% - Structural.mat: Contains structural connectivity data
%
% Parameters:
% - n_state: Number of brain states to identify (default: 5)
% - T_shift: Time shift for phase analysis (default: 9)
% - TR: Repetition time in seconds (default: 2)

%% Initialize workspace
clc; clear; close all;

% Define analysis parameters
N_STATES = 5;
T_SHIFT = 9;
TR = 2;
DO_FILTER = 0;  % Set to 1 to enable bandpass filtering

%% Load and prepare data
try
    load('LSD_TS_FC.mat', 'DataCorrel');
    load('Structural.mat', 'SC');
catch ME
    error('Required data files not found. Please ensure LSD_TS_FC.mat and Structural.mat are in the path.');
end

% Concatenate LSD and Placebo data
LSD_data = cellfun(@(x) x.LSD_TS, num2cell(DataCorrel), 'UniformOutput', false);
PCB_data = cellfun(@(x) x.PCB_TS, num2cell(DataCorrel), 'UniformOutput', false);
LSD_cat = horzcat(LSD_data{:});
PCB_cat = horzcat(PCB_data{:});
example_data = horzcat(LSD_cat, PCB_cat);

%% Signal Processing
% Setup bandpass filter
fnq = 1/(2*TR);                    % Nyquist frequency
flp = 0.01;                        % lowpass frequency
fhi = 0.1;                         % highpass frequency
Wn = [flp/fnq fhi/fnq];           % normalized frequency range
[bfilt2, afilt2] = butter(2, Wn);  % 2nd order Butterworth filter

% Z-score time series
TS = zscore(example_data, [], 2);
[N, L] = size(TS);
Phases = zeros(N, L);

% Compute Hilbert phase for each time series
for seed = 1:N
    % Detrend and demean
    x = demean(detrend(TS(seed,:)));
    
    % Mirror padding to reduce edge effects
    xp = [x(end:-1:1) x x(end:-1:1)];
    
    % Apply filter if enabled
    if DO_FILTER
        filter_xp = filtfilt(bfilt2, afilt2, xp);
    else
        filter_xp = xp;
    end
    
    % Extract central portion and compute Hilbert transform
    timeseriedata(seed,:) = filter_xp((length(x)+1):(2*length(x)));
    Xanalytic = hilbert(demean(timeseriedata(seed,:)));
    Phases(seed,:) = angle(Xanalytic);
end

%% Phase Pattern Analysis
T = (T_SHIFT+1):(size(Phases,2)-T_SHIFT);
Isubdiag = find(tril(ones(N),-1));
pattern = zeros(length(Isubdiag), length(T));

% Compute phase patterns
for t = 1:length(T)
    % Global synchronization
    kudata = sum(complex(cos(Phases(:,T(t))), sin(Phases(:,T(t))))/N);
    syncdata(t) = abs(kudata);
    
    % Pairwise phase relationships
    patt = zeros(N);
    for i = 1:N
        for j = 1:i-1
            patt(i,j) = cos(adif(Phases(i,T(t)), Phases(j,T(t))));
        end
    end
    pattern(:,t) = patt(Isubdiag);
end

%% Data Cleaning
all_pattern2D = pattern';

% Remove empty patterns
good_pattern = sum(abs(all_pattern2D), 2) > 0;
all_pattern2D = all_pattern2D(good_pattern, :);

% Remove outliers using cityblock distance
D = squareform(pdist(all_pattern2D, 'cityblock'));
D = zscore(mean(D));
good_pattern = D < 3;
all_pattern2D = all_pattern2D(good_pattern, :);

%% K-means Clustering
% Configure clustering parameters
opts = statset('Display', 'final', 'MaxIter', 100, 'UseParallel', 1);
[cidx_Pha, ctrs_Pha, sum_D_Pha] = kmeans(all_pattern2D, N_STATES, ...
    'Distance', 'cityblock', 'Replicates', 5, 'Options', opts);

%% Analysis of Brain States
% Compute correlation with structural connectivity
VC = im2double(SC(:));
CCA = zeros(1, N_STATES);
for i = 1:N_STATES
    QQ = squareform(ctrs_Pha(i,:));
    VA = QQ(:);
    MC = corrcoef(VA, VC);
    CCA(i) = MC(1,2);
end
[B, I] = sort(CCA);

% Calculate state probabilities
lsd_length = length(LSD_cat);
for bst = 1:N_STATES
    rate(bst) = sum(cidx_Pha==I(bst))/(L-2*T_SHIFT);
    ratea(bst) = sum(cidx_Pha(1:lsd_length)==I(bst))/(L-2*T_SHIFT);
    rateb(bst) = sum(cidx_Pha((lsd_length+1):end)==I(bst))/(L-2*T_SHIFT);
end

%% Visualization
figure('Position', [100 100 1200 800]);

% Plot brain states
for i = 1:N_STATES
    subplot(3, N_STATES+1, i)
    imagesc(squareform(ctrs_Pha(I(i),:)), [-1 1])
    colormap(jet)
    axis square
    title(sprintf('State %d (SFC: %.3f)', i, B(i)))
    colorbar
end

% Plot structural connectivity
subplot(3, N_STATES+1, N_STATES+1)
imagesc(SC, [-1 1])
axis square
title('Structural Connectivity')
colorbar

% Plot state probabilities
subplot(3, N_STATES+1, 2*N_STATES+2+1)
bar(ratea./sum(ratea))
ylim([0 0.5])
ylabel('Probability (LSD)')
xlabel('Brain State')

subplot(3, N_STATES+1, 2*N_STATES+2+2)
bar(rateb./sum(rateb))
ylim([0 0.5])
ylabel('Probability (Placebo)')
xlabel('Brain State')

% Plot correlation matrices for each subject
figure('Position', [100 100 1200 1500]);
n_subjects = length(DataCorrel);
for k = 1:n_subjects
    subplot(ceil(n_subjects/3), 3, k)
    imagesc(corrcoef(DataCorrel(k).LSD_TS'))
    colormap(jet)
    axis square
    title(sprintf('Subject %d - LSD', k))
    colorbar
end

%% Helper Functions
function y = demean(x)
    % Remove mean from signal
    y = x - mean(x);
end

function d = adif(x,y)
    % Calculate angular difference
    d = mod(x-y+pi, 2*pi) - pi;
end 