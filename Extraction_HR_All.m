myFiles = dir(fullfile('Data/','*.xlsx'));
disp('Reading all excel files in ./Data/');
disp('Heartrate based on REFERENCE SIGNAL');
for i = 1:length(myFiles)
    baseFileName = myFiles(i).name;
    fprintf('%s = %0.0f \n', baseFileName(1:end-5), showRefHR(baseFileName, false));
end

% return HR using FFT based on reference signal
% params:   fileName - the excel file name
%           flag - shows raw data plot & amplitude spectrum of ref signal
function HR = showRefHR(fileName, flag)
    % load data from excel file
    % col B is sensor data
    % col C is reference data
    loadPVDF = readtable(append('./Data/',fileName),'Range','B:C');
    PVDFdata = table2array(loadPVDF);
    senRaw = PVDFdata(1:end,1);
    refRaw = PVDFdata(1:end,2);
    
    % length of data array, 30 sec
    L = length(senRaw);
    % sample frequency, either 5kHz or 1kHz
    Fs = L/30;
    % sample period
    T = 1/Fs;
    
    if flag
        % show raw data
        figure; tiledlayout(2,1); 
        nexttile;
        plot(senRaw); title('Sensor Raw Data');
        nexttile;
        plot(refRaw); title('Reference Raw Data');
    end
    
    % FFT on reference signal
    Y = fft(refRaw);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    if flag
        % show amplitude spectrum
        figure; plot(f,P1) 
        title("Single-Sided Amplitude Spectrum of Reference Raw Signal")
        xlabel("f (Hz)")
        ylabel("|P1(f)|")
    end

    % calculate HR based on reference signal
    [~, reffL] = max(P1(2:200));
    reffHR = reffL*Fs/L;
    HR = 60/(1/(reffHR));
end



