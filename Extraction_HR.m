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

% bandpass filter for the sensor signal
% acceptance range 0.83Hz to 2.5Hz
% sampling range 5kHz or 1kHz, depends on data array
senFlt = bandpass(senRaw, [0.83 2.5], Fs);

% show raw data
figure; tiledlayout(3,1); 
nexttile;
plot(senRaw); title('Sensor Raw Data');
nexttile;
plot(senFlt); title('Sensor Filtered Data');
nexttile;
plot(refRaw); title('Reference Raw Data');

% FFT on reference signal
refY = fft(refRaw);
refP2 = abs(refY/L);
refP1 = refP2(1:L/2+1);
refP1(2:end-1) = 2*refP1(2:end-1);
reff = Fs*(0:(L/2))/L;

% FFT on filtered sensor signal
senY = fft(senFlt);
senP2 = abs(senY/L);
senP1 = senP2(1:L/2+1);
senP1(2:end-1) = 2*senP1(2:end-1);
senf = Fs*(0:(L/2))/L;

% show amplitude spectrum
figure; tiledlayout(2,1); 
nexttile; plot(reff,refP1);
title("Single-Sided Amplitude Spectrum of Reference Raw Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")
nexttile; plot(senf,senP1);
title("Single-Sided Amplitude Spectrum of Filtered Sensor Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% calculate HR based on reference signal
[~, reffL] = max(refP1(2:200));
reffHR = reffL*Fs/L;
refHR = 60/(1/(reffHR))

% calculate HR based on filtered sensor signal
[~, senfL] = max(senP1(2:200));
senfHR = senfL*Fs/L;
senHR = 60/(1/(senfHR))
