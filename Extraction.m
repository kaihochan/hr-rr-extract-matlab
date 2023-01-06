loadPVDF = readtable('Sample.xlsx','Range','A1:F7501');
PVDFdata = table2array(loadPVDF);
raw = PVDFdata(1:4500,2);

Fs = 1000;          % Sample frequency
T = 1/Fs;           % Sample preiod
L = length(raw);    % Length of data arry

HR_filtered = bandpass(raw, [1.00, 3.67], Fs);   % From 60 to 220 bpm
RR_filtered = bandpass(raw, [0.20, 0.33], Fs);   % From 12 to 20 breath/min

% Plot filtered signals
figure; tiledlayout(3,1);
nexttile;
plot(raw);
xlabel('Millisecond');
ylabel('V');
title('Data captured from 3x3 coupler phase 0');
nexttile;
plot(HR_filtered);
xlabel('Millisecond');
title('Heart rate filter applied');
nexttile;
plot(RR_filtered);
xlabel('Millisecond');
title('Resperation rate filter applied');

% Frequency domain f
f = Fs*(0:(L/2))/L;

% FFT for HR and RR filtered signal
HR_FFT = fft(HR_filtered);
RR_FFT = fft(RR_filtered);

% Compute the two-sided spectrum P_XR_2side
% Then compute the single-sided spectrum P_XR based on P_XR_2side and the even-valued signal length L.
P_HR_2side = abs(HR_FFT/L);
P_HR = P_HR_2side(1:L/2+1);
P_HR(2:end-1) = 2*P_HR(2:end-1);

P_RR_2side = abs(RR_FFT/L);
P_RR = P_RR_2side(1:L/2+1);
P_RR(2:end-1) = 2*P_RR(2:end-1);

[~, HR_f] = max(P_HR);
[~, RR_f] = max(P_RR);

HR = HR_f*Fs/L*60
RR = RR_f*Fs/L*60
