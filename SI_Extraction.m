function [HR,RR,hrFlt,rrFlt] = SI_Extraction(dataArray)
%SI Extract HR and RR from Sagnac interferometer

    % length of data array
    L = length(dataArray);
    % sample frequency, either 5kHz or 1kHz
    Fs = 5000;
    % sample period
    T = 1/Fs;
    
    % bandpass filter for HR
    % acceptance range 0.83Hz to 2.5Hz
    % sampling range 5kHz or 1kHz, depends on data array
    hrFlt = bandpass(dataArray,[0.83 2.5],Fs);
    
    % bandpass filter for RR
    % acceptance range 0.14Hz to 0.58Hz
    % sampling range 2.5kHz or 0.5kHz, depends on data array
    rrFlt = bandpass(dataArray,[0.14 0.58],Fs/2);
    
    % FFT on HR filtered sensor signal
    hrY = fft(hrFlt);
    hrP2 = abs(hrY/L);
    hrP1 = hrP2(1:L/2+1);
    hrP1(2:end-1) = 2*hrP1(2:end-1);
    hrf = Fs*(0:(L/2))/L;
    
    % FFT on RR filtered sensor signal
    rrY = fft(rrFlt);
    rrP2 = abs(rrY/L);
    rrP1 = rrP2(1:L/2+1);
    rrP1(2:end-1) = 2*rrP1(2:end-1);
    rrf = Fs*(0:(L/2))/L;
    
    % calculate HR based on filtered sensor signal
    [~, hrfL] = max(hrP1(2:200));
    fHR = hrfL*Fs/L;
    HR = 60/(1/(fHR));
    
    % calculate RR based on filtered sensor signal
    [~, rrfL] = max(rrP1(2:200));
    fRR = rrfL*Fs/L;
    RR = 60/(1/(fRR));
end
