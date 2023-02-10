% load data from excel file
% col B is sensor data
% col C is reference data
loadPVDF = readtable('./Data/Sit_with_cushion_Test_1_John.xlsx','Range','B:C');
PVDFdata = table2array(loadPVDF);
senRaw = PVDFdata(1:end,1);
refRaw = PVDFdata(1:end,2);

% length of data array, 30 sec
L = length(senRaw);
% sample frequency, either 5kHz or 1kHz
Fs = L/30;
% sample period
T = 1/Fs;

% show raw data
figure; tiledlayout(2,1); 
nexttile;
plot(senRaw); title('Sensor Raw Data');
nexttile;
plot(refRaw); title('Reference Raw Data');

% FFT on reference signal
Y = fft(refRaw);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% show amplitude spectrum
figure; plot(f,P1) 
title("Single-Sided Amplitude Spectrum of Reference Raw Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% calculate HR based on reference signal
[~, reffL] = max(P1(2:200));
reffHR = reffL*Fs/L;
HR = 60/(1/(reffHR))