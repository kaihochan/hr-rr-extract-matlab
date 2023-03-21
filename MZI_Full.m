myFiles = dir(fullfile('./Data/MZI','*.xlsx'));
disp('Reading all excel files in ./Data/MZI');
disp('Heartrate and resperation rate based on REFERENCE SIGNAL');
for i = 1:length(myFiles)
    baseFileName = myFiles(i).name;
    % load data from excel file
    % col B is sensor data
    loadPVDF = readtable(append('./Data/MZI/',baseFileName),'Range','C:E');
    PVDFdata = table2array(loadPVDF);
    % calculate HR and RR based on readback data
    [HR,RR,~,~] = MZI_Extraction(PVDFdata,true);
    fprintf('HR= %0.0f\t RR=%0.0f\t %s\n',HR,RR,baseFileName(1:end-5));
end

clear;

% this function is to find the best frequency range of HR filter by bruteforce
% lower limit starts at 0.7 to 0.9, with step 0.01
% higher limit starts at 2.4 to 3.67, with step 0.01
function bruteforceFindHRFreq
    fileID = fopen('./Output/output.txt','w');

    % load data from excel file
    % col B is sensor data
    % col C is reference data
    loadPVDF = readtable('./Data/MZI/Sit_Without_Cushion_Avoid_Ref_Jason_03.xlsx','Range','B:E');
    PVDFdata = table2array(loadPVDF);
    sen1 = PVDFdata(1:end,2);
    sen2 = PVDFdata(1:end,3);
    sen3 = PVDFdata(1:end,4);
    
    % length of data array, 30 sec
    L = length(sen1);
    % sample frequency, either 5kHz or 1kHz
    Fs = L/30;
    % sample period
    T = 1/Fs;

    % demoduation on phased signal
    sen = (sen1+sen2+sen3)/3;
    a = sen1-sen;
    b = sen2-sen;
    c = sen3-sen;
    d = diff(a)*Fs;
    e = diff(b)*Fs;
    f = diff(c)*Fs;
    g = a(1:end-1).*(e-f);
    h = b(1:end-1).*(f-d);
    i = c(1:end-1).*(d-e);
    z = g+h+i;
    m = a.^2+b.^2+c.^2;
    u = z./(m(1:end-1));
    demodArray = -cumtrapz(u);
    
    for i = 0.7:0.01:0.9
        for j = 2.4:0.01:3.67
            % bandpass filter for the sensor signal
            % acceptance range i to j
            % sampling range 5kHz or 1kHz, depends on data array
            senFlt = bandpass(demodArray, [i j], Fs);
        
            % FFT on filtered sensor signal
            senY = fft(senFlt);
            senP2 = abs(senY/L);
            senP1 = senP2(1:L/2+1);
            senP1(2:end-1) = 2*senP1(2:end-1);
            senf = Fs*(0:(L/2))/L;
            
            % calculate HR based on filtered sensor signal
            [~, senfL] = max(senP1(2:200));
            senfHR = senfL*Fs/L;
            senHR = 60/(1/(senfHR));
            
            fprintf(fileID,'i=%0.2f\tj=%0.2f\tHR=%0.0f\n',i,j,senHR);
        end
    end

    fclose(fileID);
end

% this function is to find the best frequency range of RR filter by bruteforce
% lower limit starts at 0.1 to 0.2, with step 0.01
% higher limit starts at 0.4 to 0.6, with step 0.01
function bruteforceFindRRFreq
    fileID = fopen('./Output/output.txt','w');
    
    % load data from m4a file
    [refSound, refFs] = audioread('./Data/MZI/Sit_Cushion_Avoid_Ref_Jason_01.m4a');
    
    % load data from excel file
    % col B is sensor data
    % col C is reference data
    loadPVDF = readtable('./Data/MZI/Sit_Cushion_Avoid_Ref_Jason_01.xlsx','Range','B:E');
    PVDFdata = table2array(loadPVDF);
    sen1 = PVDFdata(1:end,2);
    sen2 = PVDFdata(1:end,3);
    sen3 = PVDFdata(1:end,4);
    
    % length of data array, 30 sec
    L = length(sen1);
    % sample frequency, either 5kHz or 1kHz
    Fs = L/30;
    % sample period
    T = 1/Fs;

    % demoduation on phased signal
    sen = (sen1+sen2+sen3)/3;
    a = sen1-sen;
    b = sen2-sen;
    c = sen3-sen;
    d = diff(a)*Fs;
    e = diff(b)*Fs;
    f = diff(c)*Fs;
    g = a(1:end-1).*(e-f);
    h = b(1:end-1).*(f-d);
    i = c(1:end-1).*(d-e);
    z = g+h+i;
    m = a.^2+b.^2+c.^2;
    u = z./(m(1:end-1));
    demodArray = -cumtrapz(u);
    
    for i = 0.1:0.01:0.3
        for j = 0.4:0.01:0.6
            % bandpass filter for the sensor signal
            % acceptance range i to j
            % sampling range 2.5kHz or 0.5kHz, depends on data array
            senFlt = bandpass(demodArray, [i j], Fs/2);
            
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