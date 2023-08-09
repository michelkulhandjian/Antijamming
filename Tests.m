%close all; clear all;

offsetFolder  = 0;
SNR           = [0, 15];

txPathDir  = 'Tx';
rxPathDir  = 'Rx2';


for s = 1 : length(SNR)
    snr           = SNR(s);
    if ispc
        curRxPathDir = [rxPathDir, '\snr', int2str(snr)];
    else
        curRxPathDir = [rxPathDir, '/snr', int2str(snr)];
    end
    runType       = 'SIM_CHANNEL';
    Antijamming (txPathDir, curRxPathDir, [], snr, runType);
    runType       = 'RX';
    Antijamming (txPathDir, curRxPathDir, [], [], runType);
    
end



