% load data from m4a file
[refSound, refFs] = audioread('Data\Sit_without_anything_Test_4_John.m4a');

% load data from excel file
% col B is sensor data
% col C is reference data
loadPVDF = readtable('./Data/Sit_without_anything_Test_4_John.xlsx','Range','B:C');
PVDFdata = table2array(loadPVDF);
senRaw = PVDFdata(1:end,1);

% length of data array, 30 sec
L = length(senRaw);
% sample frequency, either 5kHz or 1kHz
Fs = L/30;
% sample period
T = 1/Fs;

% bandpass filter for the sensor signal
% acceptance range 0.16Hz to 0.66Hz
% sampling range 2.5kHz or 0.5kHz, depends on data array
senFlt = bandpass(senRaw, [0.16 0.66], Fs/2);

% show raw data
figure; tiledlayout(3,1); 
nexttile;
plot(senRaw); title('Sensor Raw Data');
nexttile;
plot(senFlt); title('Sensor Filtered Data');
nexttile;
plot(refSound); title('Reference Raw Data');

% FFT on filtered sensor signal
senY = fft(senFlt);
senP2 = abs(senY/L);
senP1 = senP2(1:L/2+1);
senP1(2:end-1) = 2*senP1(2:end-1);
senf = Fs*(0:(L/2))/L;

% show amplitude spectrum
figure; plot(senf,senP1);
title("Single-Sided Amplitude Spectrum of Filtered Sensor Signal")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% calculate RR based on filtered sensor signal
[~, senfL] = max(senP1(2:200));
senfRR = senfL*Fs/L;
senRR = 60/(1/(senfRR))

clear;

% this function is serve to find the best frequency range by bruteforce
% lower limit starts at 0.1 to 0.2, with step 0.01
% higher limit starts at 0.4 to 0.6, with step 0.01
function bruteforceFindFreq
    fileID = fopen('./Output/output.txt','w');
    
    % load data from m4a file
    [refSound, refFs] = audioread('Data\Sit_without_anything_Test_4_John.m4a');
    
    % load data from excel file
    % col B is sensor data
    % col C is reference data
    loadPVDF = readtable('./Data/Sit_without_anything_Test_4_John.xlsx','Range','B:C');
    PVDFdata = table2array(loadPVDF);
    senRaw = PVDFdata(1:end,1);
    
    % length of data array, 30 sec
    L = length(senRaw);
    % sample frequency, either 5kHz or 1kHz
    Fs = L/30;
    % sample period
    T = 1/Fs;
    
    for i = 0.1:0.01:0.3
        for j = 0.4:0.01:0.6
            % bandpass filter for the sensor signal
            % acceptance range 0.16Hz to 0.66Hz
            % sampling range 2.5kHz or 0.5kHz, depends on data array
            senFlt = bandpass(senRaw, [i j], Fs/2);
            
            % FFT on filtered sensor signal
            senY = fft(senFlt);
            senP2 = abs(senY/L);
            senP1 = senP2(1:L/2+1);
            senP1(2:end-1) = 2*senP1(2:end-1);
            senf = Fs*(0:(L/2))/L;
            
            % calculate RR based on filtered sensor signal
            [~, senfL] = max(senP1(2:200));
            senfRR = senfL*Fs/L;
            senRR = 60/(1/(senfRR));
            
            fprintf(fileID,'i=%0.2f\tj=%0.2f\tRR=%0.0f\n',i,j,senRR);
        end
    end

    fclose(fileID);
end
