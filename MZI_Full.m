myFiles = dir(fullfile('./Data/MZI','*.xlsx'));
disp('Reading all excel files in ./Data/MZI');
disp('Heartrate and resperation rate based on REFERENCE SIGNAL');
for i = 1:length(myFiles)
    baseFileName = myFiles(i).name;
    % load data from excel file
    % col B is sensor data
    loadPVDF = readtable(append('./Data/MZI/',baseFileName),'Range','B:E');
    PVDFdata = table2array(loadPVDF);
    % calculate HR and RR based on readback data
    [HR,RR,refHR] = MZI_Extraction(PVDFdata,true);
    fprintf('senHR= %0.0f\t senRR= %0.0f\t refHR=%0.0f\t %s\n',HR,RR,refHR,baseFileName(1:end-5));
end

clear;
