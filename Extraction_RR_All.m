% this is the script that extract RR from all excel files with sound reference
% with excel files and sound files locate at ./Data/
myFiles = dir(fullfile('Data/','*.xlsx'));
disp('Reading all excel files in ./Data/');
disp('Respiration rate based on SENSOR SIGNAL');
for i = 1:length(myFiles)
    baseFileName = myFiles(i).name;
    if isfile(append('Data/', baseFileName(1:end-5), '.m4a'))
        fprintf('Sen=%0.0f\t %s\n', showRR(baseFileName, false), baseFileName(1:end-5));
    end
end

clear;

% return RR using FFT based on reference signal
% params:   fileName - the excel file name
%           flag - shows raw data plot & amplitude spectrum of ref signal
function senRR = showRR(filename, flag)
    % load data from m4a file
    [refSound, refFs] = audioread(append('./Data/',filename(1:end-5),'.m4a'));
    
    % load data from excel file
    % col B is sensor data
    % col C is reference data
    loadPVDF = readtable(append('./Data/',filename),'Range','B:C');
    PVDFdata = table2array(loadPVDF);
    senRaw = PVDFdata(1:end,1);
    
    % length of data array, 30 sec
    L = length(senRaw);
    % sample frequency, either 5kHz or 1kHz
    Fs = L/30;
    % sample period
    T = 1/Fs;
    
    % bandpass filter for the sensor signal
    % acceptance range 0.14Hz to 0.58Hz
    % sampling range 2.5kHz or 0.5kHz, depends on data array
    senFlt = bandpass(senRaw, [0.14 0.58], Fs/2);
    
    % FFT on filtered sensor signal
    senY = fft(senFlt);
    senP2 = abs(senY/L);
    senP1 = senP2(1:L/2+1);
    senP1(2:end-1) = 2*senP1(2:end-1);
    senf = Fs*(0:(L/2))/L;
    
    % show raw data
    figure; tiledlayout(3,1); 
    nexttile;
    plot(senRaw); title('Sensor Raw Data');
    nexttile;
    plot(senFlt); title('Sensor Filtered Data');
    nexttile;
    plot(refSound); title(append('Reference Raw Data ', filename(1:end-5)), 'Interpreter','none');
    
    % calculate RR based on filtered sensor signal
    [~, senfL] = max(senP1(2:200));
    senfRR = senfL*Fs/L;
    senRR = 60/(1/(senfRR));
    
    % show data plot when flag is true
    if flag
        % show amplitude spectrum
        figure; plot(senf,senP1);
        title("Single-Sided Amplitude Spectrum of Filtered Sensor Signal")
        xlabel("f (Hz)")
        ylabel("|P1(f)|")
    end
end
