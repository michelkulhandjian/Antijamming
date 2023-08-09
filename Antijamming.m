function Antijamming (txPathDir, rxPathDir, dataPathDir, snr, runType)

% Enable debug logs
rP.debugLog  = true;
% Initialize the System Class
pR = process; % initialize the process class

if ~exist('runType','var')
    runType   = "SIM_CHANNEL"; % Options: {"TX", "RX", "DATA_VISUALIZATION", "SIMULATE", "SIM_CHANNEL" }
end

switch runType

    case "TX"
        if rP.debugLog
            fprintf("Generating waveform Signals... \n");
        end
        start = clock;
        if ~exist('txPathDir','var')
            txPathDir = 'tx_signals4';
        end
        if rP.debugLog
            fprintf( 'txPathDir: %s\n', txPathDir);
        end
        rP.MFSK = true; rP.QAM = true; rP.LDS = true; rP.LDSO=true;
        rP.isOFDM = true; rP.isSynch = true;
        pR.generateTx(txPathDir, rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        if rP.debugLog
            fprintf( 'Finished Generating waveform signals %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
        end

    case "RX"

        if rP.debugLog
            fprintf(" Processing Received waveform Signals... \n");
        end
        start = clock;
        if ~exist('txPathDir','var')
            txPathDir = 'tx_signals4';
            %txPathDir = 'input_signals';
            txPathDir = 'Tx';
        end
        if rP.debugLog
            fprintf( 'txPathDir: %s\n', txPathDir);
        end
        if ~exist('rxPathDir','var')
            %rxPathDir = 'rx_signals';
            %rxPathDir = 'output_signals';
            rxPathDir = 'Rx_GR_OTA';
        end
        if rP.debugLog
            fprintf( 'rxPathDir: %s\n', rxPathDir);
        end

        rP.MFSK = true; rP.QAM = true; rP.LDS = true; rP.LDSO=true;
        rP.isOFDM = true; rP.isSynch = true;
        pR.processRx(txPathDir, rxPathDir, rP);

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        if rP.debugLog
            fprintf( 'Finished Processing Received waveform signals %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
        end

    case "SIM_CHANNEL"

        if rP.debugLog
            fprintf(" Processing Tx waveform Signals... \n");
        end
        start = clock;
        if ~exist('txPathDir','var')
            txPathDir = 'Tx';
        end
        if rP.debugLog
            fprintf( 'txPathDir: %s\n', txPathDir);
        end
        if ~exist('rxPathDir','var')
            rxPathDir = 'Rx2';
        end
        if rP.debugLog
            fprintf( 'rxPathDir: %s\n', rxPathDir);
        end

        if ~exist('snr','var')
            snr = 0;
        end
        if rP.debugLog
            fprintf( 'snr : %d\n', snr);
        end

        rP.MFSK = true; rP.QAM = true; rP.LDS = true; rP.LDSO=true;
        rP.isOFDM = true; rP.isSynch = true;
        pR.processSimChannel(txPathDir, rxPathDir, snr, rP);

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        if rP.debugLog
            fprintf( 'Finished Processing Received waveform signals %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
        end

    case "DATA_VISUALIZATION"
        if rP.debugLog
            fprintf(" Loading data... \n");
        end
        start = clock;
        if ~exist('dataPathDir','var')
            %dataPathDir = [pwd];
            %dataPathDir = [dataPathDir '\input_signals'];
            dataPathDir = 'output_signals';
            dataPathDir = 'tx_signals4';
        end
        if rP.debugLog
            fprintf( 'dataPathDir: %s\n', dataPathDir);
        end

        pR.visualize(dataPathDir, rP);

        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        if rP.debugLog
            fprintf( 'Finished Data visuzaliation %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
        end

    case "SIMULATE"
        if rP.debugLog
            fprintf("Start Simulation... \n");
        end
        start = clock;

        rP.MFSK = true; rP.QAM = true; rP.LDS = true; rP.LDSO=true;
        rP.isOFDM = true; rP.isSynch = true;
        pR.simulate(rP);
        time_s = (clock - start)*[0 0 24*60^2 60.^[2 1 0]]';
        elapsed = fix(mod(time_s, [0 24*60^2, 60^2, 60])./[24*60^2, 60^2, 60, 1]);
        if rP.debugLog
            fprintf( 'Finished Simulation %d days, %d hours, %d min, %d sec \n', elapsed(1), elapsed(2), elapsed(3), elapsed(4));
        end

    otherwise
        disp("Option not available");
end
end









