function [HR,RR,hrFlt,rrFlt] = MZI_Extraction(PVDFdata,showFigure)
%MZI Extract HR and RR from Mach–Zehnder interferometer

    % load data from excel file
    % col 1 is reference signal
    % col 2-4 are phased sensor signal
    sen1 = PVDFdata(1:end,1);
    sen2 = PVDFdata(1:end,2);
    sen3 = PVDFdata(1:end,3);

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
    
    % bandpass filter for HR
    % acceptance range 0.83Hz to 2.5Hz
    % sampling range 5kHz or 1kHz, depends on data array
    hrFlt = bandpass(demodArray,[0.83 2.5],Fs);
    
    % bandpass filter for RR
    % acceptance range 0.14Hz to 0.58Hz
    % sampling range 2.5kHz or 0.5kHz, depends on data array
    rrFlt = bandpass(demodArray,[0.14 0.42],Fs/2);
    
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

    % show raw signal, filtered signal and filtered amplitude spectrum
    if showFigure
        figure; tiledlayout(5,1,'TileSpacing','tight','Padding','tight'); 

        nexttile; plot(0:T:30-2*T,demodArray); 
        title('Sensor Raw Signal');
        xlabel('s')
        ylabel('V')

        nexttile; plot(0:T:30-2*T,hrFlt); 
        title('HR Filtered Sensor Signal');
        xlabel('s')
        ylabel('V')

        nexttile; plot(0:T:30-2*T,rrFlt); 
        title('RR Filtered Sensor Signal');
        xlabel('s')
        ylabel('V')

        nexttile; plot(hrf,hrP1);
        title('Single-Sided Amplitude Spectrum of HR Filtered Sensor Signal')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')

        nexttile; plot(rrf,rrP1);
        title('Single-Sided Amplitude Spectrum of RR Filtered Sensor Signal')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
    end
end
