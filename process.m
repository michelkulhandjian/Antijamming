classdef process < handle
    properties
        M
        N
        K
        d_v
        S
        LDSO
        LDSO_f
        fskMod
        fskDemod
        freqSep
        MFSKSample
        isOFDM
        isSynch
        Es_N0_dB_f  %fixed Es/N0 for the pulsed jamming
        % QAM
        AveEnQAM
        N0QAMatEsN0
        % MFSK
        AveEnMFSK
        N0MFSKatEsN0
        % LDS
        AveEnLDS
        N0LDSatEsN0
        % LDSO
        AveEnLDSO
        N0LDSOatEsN0
        dSizeAdj
        MAX_TX_BUFF_SIZE
        N_ZPAD_PRE
        TX_SCALE
        NFFT
        hNFFT
        preamble
        lts_t
    end
    methods
        %-------------------------------------------------
        function obj = process()
            obj.M          = 8;
            obj.N          = 6;
            obj.K          = 8;
            obj.d_v        = 2;
            obj.freqSep    = 100;
            obj.MFSKSample = 8;
            obj.Es_N0_dB_f = 10;  %fixed Es/N0 for the pulsed jamming
            obj.isOFDM     = true;
            obj.isSynch    = true;
            obj.dSizeAdj   = 1;
            obj.S          = [];  % LDS code set
            obj.LDSO       = [];  % LDSO code set
            obj.LDSO_f     = [];
            obj.fskMod     = [];  % MFSK modulator
            obj.fskDemod   = [];  % MFSK demodulator
            % QAM
            obj.AveEnQAM     = 1;
            obj.N0QAMatEsN0  = 1;
            % MFSK
            obj.AveEnMFSK    = 1;
            obj.N0MFSKatEsN0 = 1;
            % LDS
            obj.AveEnLDS     = 1;
            obj.N0LDSatEsN0  = 1;
            % LDSO
            obj.AveEnLDSO    = 1;
            obj.N0LDSOatEsN0 = 1;
            % Preamble parameters
            obj.MAX_TX_BUFF_SIZE = 1032;
            obj.N_ZPAD_PRE       = 64;
            obj.TX_SCALE         = 1;
            obj.NFFT             = 64;
            obj.hNFFT            = floor(obj.NFFT/2)+1;
            obj.preamble         = [];
            obj.lts_t            = [];
        end

        function setLDSO(obj, txD)
            if isfield(txD, 'LDSO')
                obj.LDSO = txD.LDSO;  % LDSO code set
            end
            if isfield(txD, 'LDSO_f')
                obj.LDSO_f = txD.LDSO_f;  % LDSO_f code set
            end
        end

        function setParameters(obj, rP)
            if isfield(rP, 'M')
                obj.M          = rP.M;
            end
            if isfield(rP, 'N')
                obj.N          = rP.N;
            end
            if isfield(rP, 'K')
                obj.K          = rP.K;
            end
            if isfield(rP, 'd_v')
                obj.d_v        = rP.d_v;
            end
            if isfield(rP, 'freqSep')
                obj.freqSep    = rP.freqSep;
            end
            if isfield(rP, 'MFSKSample')
                obj.MFSKSample    = rP.MFSKSample;
            end
            if isfield(rP, 'isOFDM')
                obj.isOFDM = rP.isOFDM;
            end
            if isfield(rP, 'isSynch')
                obj.isSynch = rP.isSynch;
            end
            if isfield(rP, 'Es_N0_dB_f')
                obj.Es_N0_dB_f = rP.Es_N0_dB_f;
            end
            if isfield(rP, 'dSizeAdj')
                obj.dSizeAdj = rP.dSizeAdj;
            end
            if isfield(rP, 'MAX_TX_BUFF_SIZE')
                obj.MAX_TX_BUFF_SIZE = rP.MAX_TX_BUFF_SIZE;
            end
        end

        % Configure modulation parameters
        function configureMod(obj, rP)

            M          = obj.M;
            N          = obj.N;
            K          = obj.K;
            d_v        = obj.d_v;
            freqSep    = obj.freqSep;
            MFSKSample = obj.MFSKSample;
            isOFDM     = obj.isOFDM;
            Es_N0_dB_f = obj.Es_N0_dB_f;

            %==== QAM average energy ===============
            if rP.QAM
                constel          = qammod(0:M-1,M);
                obj.AveEnQAM     = (mean(constel.*conj(constel)));
                obj.N0QAMatEsN0  = obj.AveEnQAM .* 10^(-Es_N0_dB_f/10);
            end

            %%% MFSK ===============================
            if rP.MFSK
                obj.fskMod       = comm.FSKModulator(M,freqSep, 'SamplesPerSymbol', MFSKSample);
                obj.fskDemod     = comm.FSKDemodulator(M,freqSep, 'SamplesPerSymbol', MFSKSample);
                %A_M             = de2bi(0:M-1,kBit, 'left-msb')';
                cMFSK            = step(obj.fskMod,[0:M-1]');
                obj.AveEnMFSK    = mean(cMFSK.*conj(cMFSK));
                obj.N0MFSKatEsN0 = obj.AveEnMFSK * 10^(-Es_N0_dB_f/10);
            end

            %==== LDS average energy ===============
            if rP.LDS
                S                = obj.gen_lds_array_v3(N,K,d_v);
                S                = S .* 1/sqrt(d_v);
                obj.AveEnLDS     = 1;
                obj.S            = S;
                obj.N0LDSatEsN0  = obj.AveEnLDS .* 10^(-Es_N0_dB_f/10);
            end

            %==== LDS Orthogonal average energy ===============
            if rP.LDSO
                [LDSO]           = obj.LDSUserOne(K, N, M);
                LDSO             = LDSO*diag(1./sqrt(diag(LDSO'*LDSO)));
                obj.AveEnLDSO    = 1;
                if isOFDM
                    obj.LDSO_f   = fft(LDSO);
                end
                obj.LDSO         = LDSO;
                obj.N0LDSOatEsN0 = obj.AveEnLDSO * 10^(-Es_N0_dB_f/10);
            end

        end

        % Generate Tx data for transmission
        function obj = generateTx(obj, txPathDir, rP)

            % Tx data filename
            if ispc
                txfileName = [txPathDir '\'];
                txDataName = [txPathDir '\tx_data.mat'];
            else
                txfileName = [txPathDir '/'];
                txDataName = [txPathDir '/tx_data.mat'];
            end

            if ~exist(txPathDir, 'dir')
                % Folder does not exist so create it.
                mkdir(txPathDir);
            end

            %% Parameters Setup
            K                = 8;       % Number of "users"
            N                = 6;       % Number of subcarriers
            d_v              = 2;       % Number of subcarriers per user
            M                = K;       % Modulation Order
            kBit             = log2(M); % Number of bits per symbol
            n_symb           = 2^20;    % Number of Symbols per period
            Es_N0_dB_f       = 10;      % fixed Es/N0 for the pulsed jamming
            MAX_TX_BUFF_SIZE = 1032;    % Set the preamble size

            % Save in Data file
            rP.K                = K;
            rP.N                = N;
            rP.d_v              = d_v;
            rP.M                = M;
            rP.Es_N0_dB_f       = Es_N0_dB_f;
            rP.MAX_TX_BUFF_SIZE = MAX_TX_BUFF_SIZE;
            txD.kBit            = kBit;
            txD.n_symb          = n_symb;

            obj.setParameters(rP);
            obj.configureMod(rP);
            S           = obj.S;  % LDS code set
            LDSO        = obj.LDSO;  % LDSO code set
            LDSO_f      = obj.LDSO_f;  % LDSO code set
            fskMod      = obj.fskMod;  % MFSK modulator

            % Transmitted integer values
            tx_data     = randi(M, n_symb, 1) - 1;
            txD.tx_data = tx_data;

            % create preamble
            if rP.isSynch
                [tx_preamble]   = obj.create_preamble();
            end

            if rP.QAM
                x_QAM     = qammod(tx_data,M);
                txD.x_QAM = x_QAM;
                if rP.isSynch
                    x_QAM     = [tx_preamble/max(abs(tx_preamble))*max(abs(x_QAM)); x_QAM];
                end
                obj.writeToFile([txfileName, 'QAM.bin'], x_QAM);
            end

            if rP.MFSK
                x_MFSK     = step(fskMod,tx_data);
                txD.x_MFSK = x_MFSK;
                if rP.isSynch
                    x_MFSK     = [tx_preamble/max(abs(tx_preamble))*max(abs(x_MFSK)); x_MFSK];
                end
                obj.writeToFile([txfileName, 'MFSK.bin'], x_MFSK);
            end

            if rP.LDS
                x_LDS     = obj.lds_mod(tx_data,S,"is_ofdm", rP.isOFDM);
                txD.x_LDS = x_LDS;
                if rP.isSynch
                    x_LDS     = [tx_preamble/max(abs(tx_preamble))*max(abs(x_LDS)); x_LDS];
                end
                obj.writeToFile([txfileName, 'LDS.bin'], x_LDS);
            end

            if rP.LDSO
                if rP.isOFDM
                    x_LDSO = ifft( LDSO_f(:,tx_data(:, 1)+1));
                else
                    x_LDSO = LDSO(:,tx_data(:, 1)+1);
                end
                txD.x_LDSO = x_LDSO;
                txD.LDSO = LDSO; txD.LDSO_f = LDSO_f;
                if rP.isSynch
                    x_LDSO     = [tx_preamble/max(abs(tx_preamble))*max(abs(x_LDSO(:))); x_LDSO(:)];
                end
                obj.writeToFile([txfileName, 'LDSO.bin'], x_LDSO(:));
            end

            txD.rP  = rP;
            % Saving the .mat file for txData
            save (txDataName, 'txD');

        end % generateTx

        % Process Rx data
        function obj = processRx(obj, txPathDir, rxPathDir, rP)

            if exist('txPathDir','var')
                % third parameter does exist
                if isempty(txPathDir)
                    txPathDir = [pwd];
                elseif ~exist(txPathDir,'dir')
                    msg = ['The txPathDir directory: ' txPathDir ',   not found, please verify the directory'];
                    error(msg)
                end
            else
                txPathDir = [pwd];
            end

            if exist('rxPathDir','var')
                % third parameter does exist
                if isempty(rxPathDir)
                    rxPathDir = [pwd];
                elseif ~exist(rxPathDir,'dir')
                    msg = ['The txPathDir directory: ' rxPathDir ',   not found, please verify the directory'];
                    error(msg)
                end
            else
                rxPathDir = [pwd];
            end

            % Rx data filename
            if ispc
                txDataName = [txPathDir '\tx_data.mat'];
                rxDataName = [rxPathDir '\rx_data.mat'];
            else
                txDataName = [txPathDir '/tx_data.mat'];
                rxDataName = [rxPathDir '/rx_data.mat'];
            end

            isTxData   = true;
            if  ~isfile(txDataName)
                % No tx Data is present
                isTxData = false;
                obj.setParameters(rP);
                obj.configureMod(rP);
            else
                % tx configuration and data
                load(txDataName);
                obj.setParameters(txD.rP);
                obj.configureMod(txD.rP);
                obj.setLDSO(txD);
                tx_data    = txD.tx_data;
                n_symb     = txD.n_symb;
            end
            % Set parameters
            M          = obj.M;
            %N          = obj.N;
            %K          = obj.K;
            %d_v        = obj.d_v;
            %freqSep    = obj.freqSep;
            isOFDM     = obj.isOFDM;
            %Es_N0_dB_f = obj.Es_N0_dB_f;
            S          = obj.S;  % LDS code set
            LDSO       = obj.LDSO;  % LDSO code set
            LDSO_f     = obj.LDSO_f;  % LDSO code set
            fskDemod   = obj.fskDemod;  % MFSK demodulator

            % Parameters
            Nspec = 3076;
            rP.fs = 1;

            % create preamble
            if rP.isSynch
                [tx_preamble]   = obj.create_preamble();
            end

            if rP.debugLog
                fprintf(' Start data loading process... \n');
            end

            myFiles = dir(fullfile(rxPathDir,'**','*.bin')); %gets all bin files in struct
            NmF = length(myFiles);

            if NmF > 0
                for k = 1:NmF
                    baseFileName = myFiles(k).name;
                    [ext, name] = fileparts(baseFileName);
                    if ispc
                        pathparts = strsplit(myFiles(k).folder,'\');
                        NameFile = [myFiles(k).folder '\' myFiles(k).name];
                    else
                        pathparts = strsplit(myFiles(k).folder,'/');
                        NameFile = [myFiles(k).folder '/' myFiles(k).name];
                    end
                    if rP.debugLog
                        fprintf(' file name %s \n', baseFileName);
                    end
                    data = obj.readFile( NameFile);

                    % remove preamble
                    if rP.isSynch
                        data = obj.extractPreamble (data);
                    end

                    if ~isempty(data)
                        [figSig, figSpectrogram] = obj.captureSig ( data(1:Nspec), rP);
                        set(figSpectrogram, 'Visible', 'on');
                        titleStr = sprintf(' File: %s ', baseFileName);
                        title(titleStr); xtickformat('%.1f'); %title(titleStr);

                        set(figSig, 'Visible', 'on');
                        title(titleStr);

                        % Demodulate data
                        if strcmp(name,'QAM')
                            Lqam = length(data);
                            rx_data_QAM     = qamdemod(data(1:Lqam-mod(Lqam, M)),M);
                            rxD.rx_data_QAM = rx_data_QAM;
                            if isTxData
                                Nsym1 = min(length(rx_data_QAM), n_symb);
                                errQAM     = sum(tx_data(1:Nsym1) ~= rx_data_QAM(1:Nsym1));
                                rxD.errQAM = errQAM;
                                if rP.debugLog
                                    fprintf(' QAM Error %d out of %d, or %.2f%% \n', errQAM, Nsym1, errQAM/Nsym1*100);
                                end
                            else
                                if rP.debugLog
                                    fprintf(' No Tx data cannot compute error \n');
                                end
                            end
                        end

                        if strcmp(name,'MFSK')
                            Lmfsk = length(data);
                            rx_data_MFSK     = step(fskDemod,data(1:Lmfsk-mod(Lmfsk, M)));
                            if isTxData
                                Nsym1 = min(length(rx_data_MFSK), n_symb);
                                errMFSK     = sum(tx_data(1:Nsym1) ~= rx_data_MFSK(1:Nsym1));
                                rxD.errMFSK = errMFSK;
                                if rP.debugLog
                                    fprintf(' MFSK Error %d out of %d, or %.2f%% \n', errMFSK, Nsym1, errMFSK/Nsym1*100);
                                end
                            else
                                if rP.debugLog
                                    fprintf(' No Tx data cannot compute error \n');
                                end
                            end
                            rxD.rx_data_MFSK = rx_data_MFSK;
                        end

                        if strcmp(name,'LDS')
                            rx_data_LDS     = obj.lds_demod(data, S, "is_ofdm", isOFDM);
                            if isTxData
                                Nsym1 = min(length(rx_data_LDS), n_symb);
                                errLDS     = sum(tx_data(1:Nsym1) ~= rx_data_LDS(1:Nsym1));
                                rxD.errLDS = errLDS;
                                if rP.debugLog
                                    fprintf(' LDS Error %d out of %d, or %.2f%% \n', errLDS, Nsym1, errLDS/Nsym1*100 );
                                end
                            else
                                if rP.debugLog
                                    fprintf(' No Tx data cannot compute error \n');
                                end
                            end
                            rxD.rx_data_LDS = rx_data_LDS;
                        end

                        if strcmp(name,'LDSO')
                            %ldso = reshape(data,size(S,2),[]);
                            ldso = buffer(data,size(S,2));
                            if isOFDM
                                y_LDSO = fft(ldso);
                                rx_data_LDSO = real(LDSO_f'*y_LDSO);
                            else
                                rx_data_LDSO = real(LDSO'*y_LDSO);
                            end
                            [mA, rx_data_LDSO] = max(rx_data_LDSO);
                            rx_data_LDSO = rx_data_LDSO.'-1;

                            if isTxData
                                Nsym1 = min(length(rx_data_LDSO), n_symb);
                                errLDSO     = sum(tx_data(1:Nsym1) ~= rx_data_LDSO(1:Nsym1));
                                rxD.errLDSO = errLDSO;
                                if rP.debugLog
                                    fprintf(' LDSO Error %d out of %d, or %.2f%% \n', errLDSO, Nsym1, errLDSO/Nsym1*100);
                                end
                            else
                                if rP.debugLog
                                    fprintf(' No Tx data cannot compute error \n');
                                end
                            end
                            rxD.rx_data_LDSO = rx_data_LDSO;
                        end
                    else
                        if rP.debugLog
                            fprintf(' Empty data... \n');
                        end
                    end
                end
            else
                if rP.debugLog
                    fprintf(' No Files are found... \n');
                end
            end

            rxD.rP = rP;
            % Saving the .mat file for rxData
            save (rxDataName, 'rxD');

        end % processRx

        % Process SimChannel data
        function obj = processSimChannel(obj, txPathDir, rxPathDir, snr, rP)

            if exist('txPathDir','var')
                % third parameter does exist
                if isempty(txPathDir)
                    txPathDir = [pwd];
                elseif ~exist(txPathDir,'dir')
                    msg = ['The txPathDir directory: ' txPathDir ',   not found, please verify the directory'];
                    error(msg)
                end
            else
                txPathDir = [pwd];
            end

            if exist('rxPathDir','var')
                % third parameter does exist
                if isempty(rxPathDir)
                    rxPathDir = [pwd];
                elseif ~exist(rxPathDir,'dir')
                    % Folder does not exist so create it.
                    mkdir(rxPathDir);
                end
            else
                rxPathDir = [pwd];
            end

            % Rx data filename
            if ispc
                rxfileName = [rxPathDir '\'];
                txDataName = [txPathDir '\tx_data.mat'];
                rxDataName = [rxPathDir '\tx_data.mat'];
            else
                rxfileName = [rxPathDir '/'];
                txDataName = [txPathDir '/tx_data.mat'];
                rxDataName = [rxPathDir '/tx_data.mat'];
            end

            isTxData   = true;
            if  ~isfile(txDataName)
                % No tx Data is present
                isTxData = false;
                rP.Es_N0_dB_f = snr;
                obj.setParameters(rP);
                obj.configureMod(rP);
            else
                % tx configuration and data
                load(txDataName);
                txD.rP.Es_N0_dB_f = snr;
                obj.setParameters(txD.rP);
                obj.configureMod(txD.rP);
                obj.setLDSO(txD);
                tx_data    = txD.tx_data;
                n_symb     = txD.n_symb;
            end
            % Set parameters
            M          = obj.M;
            isOFDM     = obj.isOFDM;
            %Es_N0_dB_f = obj.Es_N0_dB_f;
            S          = obj.S;  % LDS code set
            LDSO       = obj.LDSO;  % LDSO code set
            LDSO_f     = obj.LDSO_f;  % LDSO code set
            fskDemod   = obj.fskDemod;  % MFSK demodulator

            % Parameters
            Nspec = 3076;
            rP.fs = 1;

            % % create preamble
            % if rP.isSynch
            %     [tx_preamble]   = obj.create_preamble();
            % end

            if rP.debugLog
                fprintf(' Start data loading process... \n');
            end

            myFiles = dir(fullfile(txPathDir,'**','*.bin')); %gets all bin files in struct
            NmF = length(myFiles);

            if NmF > 0
                for k = 1:NmF
                    baseFileName = myFiles(k).name;
                    [ext, name] = fileparts(baseFileName);
                    if ispc
                        pathparts = strsplit(myFiles(k).folder,'\');
                        NameFile = [myFiles(k).folder '\' myFiles(k).name];
                    else
                        pathparts = strsplit(myFiles(k).folder,'/');
                        NameFile = [myFiles(k).folder '/' myFiles(k).name];
                    end
                    if rP.debugLog
                        fprintf(' file name %s \n', baseFileName);
                    end
                    data = obj.readFile( NameFile);

                    if ~isempty(data)
                        [figSig, figSpectrogram] = obj.captureSig ( data(1:Nspec), rP);
                        set(figSpectrogram, 'Visible', 'on');
                        titleStr = sprintf(' File: %s ', baseFileName);
                        title(titleStr); xtickformat('%.1f'); %title(titleStr);

                        set(figSig, 'Visible', 'on');
                        title(titleStr);

                        % Demodulate data
                        if strcmp(name,'QAM')
                            Lqam = length(data);
                            N0QAMJam = obj.AveEnQAM .* 10^(-snr/10);
                            NoiseQAMatEsN0 =sqrt(obj.N0QAMatEsN0/2).*(randn(Lqam,1) + 1i.*randn(Lqam,1));
                            %NoiseQAMJam = sqrt(N0QAMJam/2)*jam_en_v.*(randn(n_symb,1) + 1i.*randn(n_symb,1));
                            x_QAM = data + NoiseQAMatEsN0; % + 0*NoiseQAMJam;;
                            obj.writeToFile([rxfileName, 'QAM.bin'], x_QAM);
                        end

                        if strcmp(name,'MFSK')
                            Lmfsk = length(data);
                            N0MFSKJam = obj.AveEnMFSK .* 10^(-snr/10);
                            NoiseMFSKatEsN0 = sqrt(obj.N0MFSKatEsN0/2)*(randn(Lmfsk, 1) + 1i*randn(Lmfsk, 1));  % when using pulsed jamming, we fixed the AWGN noise at certain Es/N0
                            %NoiseMFSKJam = sqrt(N0MFSKJam/2)*(randn(obj.MFSKSample, n_symb) + 1i.*randn(obj.MFSKSample,n_symb)).*jam_en_v.'; % pulsed jammer Es/Nj
                            x_MFSK = data + NoiseMFSKatEsN0; % + NoiseMFSKJam(:);
                            obj.writeToFile([rxfileName, 'MFSK.bin'], x_MFSK);
                        end

                        if strcmp(name,'LDS')
                            Llds = length(data);
                            N0LDSJam = obj.AveEnLDS .* 10^(-snr/10);
                            NoiseLDSatEsN0 =sqrt(obj.N0LDSatEsN0/2).*(randn(Llds,1) + 1i.*randn(Llds,1));
                            %NoiseLDSJam = sqrt(N0LDSJam/2)*jam_en_v.*(randn(n_symb,N) + 1i.*randn(n_symb,N));
                            %NoiseLDSJam = NoiseLDSJam.';
                            x_LDS = data + NoiseLDSatEsN0(:); % + NoiseLDSJam(:);
                            obj.writeToFile([rxfileName, 'LDS.bin'], x_LDS);
                        end

                        if strcmp(name,'LDSO')
                            Lldso = length(data);
                            N0LDSOJam = obj.AveEnLDSO .* 10^(-snr/10);
                            NoiseLDSOatEsN0 = sqrt(obj.N0LDSOatEsN0/2)*(randn(Lldso,1) + 1i*randn(Lldso,1));  % when using pulsed jamming, we fixed the AWGN noise at certain Es/N0
                            %NoiseLDSOJam = sqrt(N/K)*sqrt(N0LDSOJam/2)*(randn(K, n_symb) + 1i.*randn(K,n_symb)).*jam_en_v.'; % pulsed jammer Es/Nj
                            x_LDSO = data + NoiseLDSOatEsN0; % + NoiseLDSOJam;
                            obj.writeToFile([rxfileName, 'LDSO.bin'], x_LDSO(:));
                        end
                    else
                        if rP.debugLog
                            fprintf(' Empty data... \n');
                        end
                    end
                end

                txD.rP  = rP;
                % Saving the .mat file for txData
                save (rxDataName, 'txD');
            else
                if rP.debugLog
                    fprintf(' No Files are found... \n');
                end
            end
        end % processSimChannel

        % Visualize Data
        function visualize(obj, dataPathDir, rP);

            if exist('dataPathDir','var')
                % third parameter does exist
                if isempty(dataPathDir)
                    dataPathDir = [pwd];
                elseif ~exist(dataPathDir,'dir')
                    msg = ['The data directory: ' dataPathDir ',   not found, please verify the directory'];
                    error(msg)
                end
            else
                dataPathDir = [pwd];
            end


            % Output Folder
            if ispc
                outPutFolder = [dataPathDir '\visualize'];
            else
                outPutFolder = [dataPathDir '/visualize'];
            end

            if ~exist(outPutFolder, 'dir')
                % Folder does not exist so create it.
                mkdir(outPutFolder);
            end

            % Parameters
            Nspec = 3076;
            rP.fs = 1;

            if rP.debugLog
                fprintf(' Start data loading process... \n');
            end

            myFiles = dir(fullfile(dataPathDir,'**','*.bin')); %gets all bin files in struct
            NmF = length(myFiles);

            if NmF > 0
                for k = 1:NmF
                    baseFileName = myFiles(k).name;
                    [ext, name] = fileparts(baseFileName);
                    if ispc
                        pathparts = strsplit(myFiles(k).folder,'\');
                        NameFile = [myFiles(k).folder '\' myFiles(k).name];
                    else
                        pathparts = strsplit(myFiles(k).folder,'/');
                        NameFile = [myFiles(k).folder '/' myFiles(k).name];
                    end
                    if rP.debugLog
                        fprintf(' file name %s \n', baseFileName);
                    end
                    data = obj.readFile( NameFile);

                    if ~isempty(data)
                        [figSig, figSpectrogram] = obj.captureSig ( data(1:Nspec), rP);
                        set(figSpectrogram, 'Visible', 'on');
                        titleStr = sprintf(' File: %s ', baseFileName);
                        title(titleStr); xtickformat('%.1f'); %title(titleStr);

                        set(figSig, 'Visible', 'on');
                        title(titleStr);

                        if ispc
                            saveas(figSig,[outPutFolder '\' 'SIG_' name  '.png']);
                            saveas(figSpectrogram,[outPutFolder '\' 'SPECTROGRAM_' name  '.png']);
                        else
                            saveas(figSig,[outPutFolder '/' 'SIG_'  name '.png']);
                            saveas(figSpectrogram,[outPutFolder '/' 'SPECTROGRAM_'  name '.png']);
                        end
                    else
                        if rP.debugLog
                            fprintf(' Empty data... \n');
                        end
                    end
                end

            else
                if rP.debugLog
                    fprintf(' No Files are found... \n');
                end
            end

        end % Visuzalize Data

        % Simulation
        function obj = simulate(obj, rP)

            %% Parameters Setup
            K                = 8;       % Number of "users"
            N                = 6;       % Number of subcarriers
            d_v              = 2;       % Number of subcarriers per user
            M                = K;       % Modulation Order
            kBit             = log2(M); % Number of bits per symbol
            n_symb           = 2^20;    % Number of Symbols per period
            n_symb_burst     = 256;
            n_symb_pilot     = K;
            MAX_TX_BUFF_SIZE = 1032;    % Set the preamble size

            %%% Eb/N0 =========
            Es_N0_dB_f       = 10;      % fixed Es/N0 for the pulsed jamming
            startEbN0        = 1;      % start of Eb/N0(dB)
            stepEbN0         = 2;       % increment of Eb/N0(dB)
            endEbN0          = 20;      % end of EbN0 (dB)
            Eb_N0_dB_J       = [startEbN0:stepEbN0:endEbN0]; % multiple Es/N0 values
            Es_N0_dB_J       = Eb_N0_dB_J + 10*log10(kBit);
            LEbN0J           = length(Eb_N0_dB_J);
            rho              = 0.5;   % Jamming Probability

            % Save in Data file
            rP.K                = K;
            rP.N                = N;
            rP.d_v              = d_v;
            rP.M                = M;
            rP.Es_N0_dB_f       = Es_N0_dB_f;
            rP.MAX_TX_BUFF_SIZE = MAX_TX_BUFF_SIZE;
            txD.kBit            = kBit;
            txD.n_symb          = n_symb;

            obj.setParameters(rP);
            obj.configureMod(rP);
            S           = obj.S;  % LDS code set
            LDSO        = obj.LDSO;  % LDSO code set
            LDSO_f      = obj.LDSO_f;  % LDSO code set
            fskMod      = obj.fskMod;  % MFSK modulator

            if rP.QAM
                AveEnQAM    = obj.AveEnQAM;
                N0QAMatEsN0 = obj.N0QAMatEsN0;
                ber_QAM     = zeros(1, LEbN0J);
                ser_QAM     = zeros(1, LEbN0J);
            end

            if rP.MFSK
                fskMod       = obj.fskMod;    % MFSK Modulator
                fskDemod     = obj.fskDemod;  % MFSK demodulator
                MFSKSample   = obj.MFSKSample;
                AveEnMFSK    = obj.AveEnMFSK;
                N0MFSKatEsN0 = obj.N0MFSKatEsN0;
                ber_MFSK     = zeros(1, LEbN0J);
                ser_MFSK     = zeros(1, LEbN0J);
            end

            if rP.LDS
                AveEnLDS     = obj.AveEnLDS;
                N0LDSatEsN0  = obj.N0LDSatEsN0;
                ber_LDS      = zeros(1, LEbN0J);
                ser_LDS      = zeros(1, LEbN0J);
            end

            if rP.LDSO
                AveEnLDSO    = obj.AveEnLDSO;
                N0LDSOatEsN0 = obj.N0LDSOatEsN0;
                ber_LDSO     = zeros(1, LEbN0J);
                ser_LDSO     = zeros(1, LEbN0J);
            end

            jam_en_v = rand(n_symb,1)<= rho;
            jam_opts.rho = rho;
            jam_opts.pulse_len = N;

            for ii = 1:LEbN0J

                % Transmitted integer values, insert pilot symbols  with K per burst
                ts = 1:M;
                tx_data = randi(M, n_symb, 1);
                n_burst = n_symb / n_symb_burst;
                n_symb_between_pilots = n_symb_burst / n_symb_pilot;
                for jj = 1:n_burst
                    for kk = 1:n_symb_pilot
                        tx_data(n_symb_burst*(jj-1) + n_symb_between_pilots*(kk-1) + 1) = kk;
                    end
                end
                tx_data = tx_data - 1;
                %Tx_bits = A_M(:,tx_data+1);

                % Transmitted integer values
                %tx_data     = randi(M, n_symb, 1) - 1;

                if rP.QAM
                    N0QAMJam = AveEnQAM .* 10^(-Es_N0_dB_J(ii)/10);
                    NoiseQAMatEsN0 =sqrt(N0QAMatEsN0/2).*(randn(n_symb,1) + 1i.*randn(n_symb,1));
                    NoiseQAMJam = sqrt(N0QAMJam/2)*jam_en_v.*(randn(n_symb,1) + 1i.*randn(n_symb,1));
                    x_QAM = qammod(tx_data,M);
                    y_QAM = x_QAM + NoiseQAMatEsN0 + NoiseQAMJam;

                    rx_data_QAM = qamdemod(y_QAM,M);
                    ser_QAM(ii) = sum(tx_data ~= rx_data_QAM) / n_symb;
                    [~,ber_QAM(ii)] = biterr(tx_data, rx_data_QAM);
                    fprintf('QAM EbN0 Jam: %.2f dB  |  SER: %.2E  |  BER: %.2E \n', Eb_N0_dB_J(ii), ser_QAM(ii), ber_QAM(ii));
                end

                if rP.MFSK
                    N0MFSKJam = AveEnMFSK .* 10^(-Es_N0_dB_J(ii)/10);
                    NoiseMFSKatEsN0 = sqrt(N0MFSKatEsN0/2)*(randn(MFSKSample*n_symb, 1) + 1i*randn(MFSKSample*n_symb, 1));  % when using pulsed jamming, we fixed the AWGN noise at certain Es/N0
                    NoiseMFSKJam = sqrt(N0MFSKJam/2)*(randn(MFSKSample, n_symb) + 1i.*randn(MFSKSample,n_symb)).*jam_en_v.'; % pulsed jammer Es/Nj
                    x_MFSK = step(fskMod,tx_data)/sqrt(MFSKSample);
                    y_MFSK = x_MFSK + NoiseMFSKatEsN0 + NoiseMFSKJam(:);

                    rx_data_MFSK = step(fskDemod,y_MFSK);
                    ser_MFSK(ii) = sum(tx_data ~= rx_data_MFSK) / n_symb;
                    [~,ber_MFSK(ii)] = biterr(tx_data, rx_data_MFSK);
                    fprintf('MFSK EbN0 Jam: %.2f dB  |  SER: %.2E  |  BER: %.2E \n', Eb_N0_dB_J(ii), ser_MFSK(ii), ber_MFSK(ii));
                end

                if rP.LDS
                    N0LDSJam = AveEnLDS .* 10^(-Es_N0_dB_J(ii)/10);
                    NoiseLDSatEsN0 =sqrt(N0LDSatEsN0/2).*(randn(n_symb,N) + 1i.*randn(n_symb,N));
                    NoiseLDSJam = sqrt(N0LDSJam/2)*jam_en_v.*(randn(n_symb,N) + 1i.*randn(n_symb,N));
                    NoiseLDSJam = NoiseLDSJam.';

                    x_LDS = obj.lds_mod(tx_data,S,"is_ofdm",rP.isOFDM);
                    y_LDS = x_LDS + NoiseLDSatEsN0(:) + NoiseLDSJam(:);

                    rx_data_LDS = obj.lds_demod(y_LDS,S,"is_ofdm",rP.isOFDM);
                    ser_LDS(ii)  = sum(tx_data ~= rx_data_LDS) / n_symb;
                    [~,ber_LDS(ii)] = biterr(tx_data, rx_data_LDS);
                    fprintf('LDS EbN0 Jam: %.2f dB  |  SER: %.2E  |  BER: %.2E \n', Eb_N0_dB_J(ii), ser_LDS(ii), ber_LDS(ii));

                end

                if rP.LDSO

                    N0LDSOJam = AveEnLDSO .* 10^(-Es_N0_dB_J(ii)/10);
                    NoiseLDSOatEsN0 = sqrt(N0LDSOatEsN0/2)*(randn(K,n_symb) + 1i*randn(K,n_symb));  % when using pulsed jamming, we fixed the AWGN noise at certain Es/N0
                    NoiseLDSOJam = sqrt(N0LDSOJam/2)*(randn(K, n_symb) + 1i.*randn(K,n_symb)).*jam_en_v.'; % pulsed jammer Es/Nj
                    if rP.isOFDM
                        x_LDSO = ifft( LDSO_f(:,tx_data(:, 1)+1));
                        y_LDSO = x_LDSO + NoiseLDSOatEsN0 + NoiseLDSOJam;
                        y_LDSO = fft(y_LDSO);
                        rx_data_LDSO = real(LDSO_f'*y_LDSO);
                    else
                        x_LDSO = LDSO(:,tx_data(:, 1)+1);
                        y_LDSO = x_LDSO + NoiseLDSOatEsN0 + NoiseLDSOJam;
                        rx_data_LDSO = real(LDSO'*y_LDSO);
                    end

                    [mA, rx_data_LDSO] = max(rx_data_LDSO);
                    rx_data_LDSO = rx_data_LDSO.'-1;

                    ser_LDSO(ii) = sum(tx_data ~= rx_data_LDSO) / n_symb;
                    [~,ber_LDSO(ii)] = biterr(tx_data, rx_data_LDSO);
                    fprintf('LDSO EbN0 Jam: %.2f dB  |  SER: %.2E  |  BER: %.2E \n', Eb_N0_dB_J(ii), ser_LDSO(ii), ber_LDSO(ii));
                end
            end

            f1 = figure;
            indL = 1;
            if rP.QAM
                plot(Eb_N0_dB_J, ber_QAM,'-o')
                legendInfo{indL} = sprintf('%i-QAM',K);
                indL = indL + 1;
                hold on;
            end
            if rP.MFSK
                plot(Eb_N0_dB_J, ber_MFSK,'--')
                legendInfo{indL} = sprintf('%i-MFSK',K);
                indL = indL + 1;
                hold on;
            end

            if rP.LDS
                plot(Eb_N0_dB_J, ber_LDS,'-*')
                legendInfo{indL} = 'LDS';
                indL = indL + 1;
                hold on;
            end
            if rP.LDSO
                plot(Eb_N0_dB_J, ber_LDSO,'-v')
                legendInfo{indL} = 'LDSO';
            end

            hold off ; title("BER vs Jamming EbN0")
            ylabel('BER') ;xlabel('Pulse Jamming EbN0(dB)');
            grid on ; set(gca, 'YScale', 'log'); legend(legendInfo)
            %exportgraphics(f1,fig_path+"f_ber.eps")


            f2 = figure;  indL = 1; clear legendInfo;
            if rP.QAM
                plot(Eb_N0_dB_J, ser_QAM,'-o')
                legendInfo{indL} = sprintf('%i-QAM',K);
                indL = indL + 1;
                hold on;
            end
            if rP.MFSK
                plot(Eb_N0_dB_J, ser_MFSK,'--')
                legendInfo{indL} = sprintf('%i-MFSK',K);
                indL = indL + 1;
                hold on;
            end
            if rP.LDS
                plot(Eb_N0_dB_J, ser_LDS,'-*')
                legendInfo{indL} =  'LDS';
                indL = indL + 1;
                hold on;
            end
            if rP.LDSO
                plot(Eb_N0_dB_J, ser_LDSO,'-v')
                legendInfo{indL} = 'LDSO';
            end

            hold off; title("SER vs Jamming EbN0")
            ylabel('SER') ; xlabel('Pulse Jamming EbN0(dB)')
            set(gca, 'YScale', 'log') ; grid on; legend(legendInfo);
            %exportgraphics(f1,fig_path+"f_ser.eps")


            %f_pj = figure;
            %spectrogram(NoiseLDSJam(:));
            %exportgraphics(f_pj,fig_path+"f_pf_spec.eps")

            % Saving the .mat file for txData
            %save (txDataName, 'txD');

        end % simulate

        function [fig1, fig2] = captureSig (obj, data, rP)
            C         = viridis(45);
            sig       = data;
            numXticks = 4;

            % plot signal
            fig1 = figure('visible','off');
            plot (real(sig));
            xlabel('Samples'); ylabel('Real Mag.'); title('SIgnal');

            % Plot spectrogram
            fig2 = figure('visible','off');
            NFFT = 64;
            NOVERLAP = 60;
            WINDOW = 64;
            [~,Ftx,Ttx,PPtx] = spectrogram(sig,hanning(WINDOW),NOVERLAP,NFFT, rP.fs, 'centered', 'yaxis');
            imagesc(Ttx/1e-3, Ftx/1e6, 10*log10(PPtx)); title("spectrogram");
            xticks(linspace(Ttx(1)/1e-3, Ttx(end)/1e-3, numXticks));
            xlabel('Time (ms)'); ylabel('Frequency (MHz)');
            set(gca,'YDir','normal'); colormap(C);
        end


        % Generate Preamble
        function [tx_vec_cal] = create_preamble(obj)

            MAX_TX_BUFF_SIZE = obj.MAX_TX_BUFF_SIZE;
            N_ZPAD_PRE       = obj.N_ZPAD_PRE;
            TX_SCALE         = obj.TX_SCALE;
            NFFT             = obj.NFFT;
            hNFFT            = obj.hNFFT;

            % Preamble
            % Add a preamble LTS for fine CFO and channel estimation
            lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
            lts_t        = ifft(lts_f, NFFT);                    % time domain
            lts          = [lts_t(hNFFT:NFFT) lts_t lts_t];      % 2.5 LTS
            preamble     = lts;
            preamble     = preamble/max(abs(preamble));
            obj.preamble = preamble;
            obj.lts_t    = lts_t;

            n_samp_cal      = MAX_TX_BUFF_SIZE;
            remaining_samps = n_samp_cal - length(preamble) - N_ZPAD_PRE;
            tx_vec_cal      = [zeros(1,N_ZPAD_PRE) preamble zeros(1,remaining_samps)];
            % Scale the Tx vector to +/- 1
            tx_vec_cal      = TX_SCALE .* tx_vec_cal ./ max(abs(tx_vec_cal));
            tx_vec_cal      = tx_vec_cal.';
        end

        % Removes the preamble of the received signal (Syncronization)
        function [rx] = extractPreamble (obj, rx)

            preamble         = obj.preamble;
            lts_t            = obj.lts_t;
            MAX_TX_BUFF_SIZE = obj.MAX_TX_BUFF_SIZE;
            N_ZPAD_PRE       = obj.N_ZPAD_PRE;

            % percent of max value
            percM    = 0.8;
            % Offset values
            lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx)));
            lts_peaks = find(lts_corr > percM*max(lts_corr));
            [LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
            [lts_second_peak_index,y] = find(LTS2-LTS1 == length(lts_t));  % use size of lts+cp (80)
            % Stop if no valid correlation peak was found
            if(isempty(lts_second_peak_index))
                fprintf('SAMP OFFSET CAL: No LTS Correlation Peaks Found Exit now!\n');
                return;
            end
            peak_ind   = lts_peaks(lts_second_peak_index(1));  % Get second peak
            cal_offset = peak_ind + MAX_TX_BUFF_SIZE - length(preamble) - N_ZPAD_PRE +1;
            rx         = rx(cal_offset:end);
        end

        % LDS modulation
        function [tx_sig] = lds_mod(obj, varargin)
            p = inputParser;
            p.addRequired("tx_data",@isnumeric);
            p.addRequired("S");
            p.addParameter("is_ofdm",0,@islogical);
            p.addParameter("is_2d",0,@isnumeric);
            p.parse(varargin{:});
            Opt = p.Results;

            if Opt.is_ofdm
                N = size(Opt.S,1);
                Opt.S = ifft(Opt.S) * sqrt(N); % Compensate for IFFT gain
            end

            X = Opt.S(:,Opt.tx_data+1);
            if(Opt.is_2d)
                tx_sig = X;
            else
                tx_sig = X(:);
            end
        end

        % LDS demodulation
        function rx_data = lds_demod(obj, varargin)
            p = inputParser;
            p.addRequired("rx_sig",@isnumeric);
            p.addRequired("S");
            p.addParameter("is_ofdm",0,@islogical);
            p.parse(varargin{:});
            Opt = p.Results;


            %X = reshape(Opt.rx_sig,size(Opt.S,1),[]);
            X = buffer(Opt.rx_sig,size(Opt.S,1));

            if Opt.is_ofdm
                X = fft(X);
            end

            out = Opt.S'*X;

            [~,rx_data] = max(real(out));

            if iscolumn(Opt.rx_sig)
                rx_data = rx_data.'-1;
            else
                rx_data = rx_data - 1;
            end

        end

        % Generate LDS User One
        function [C] = LDSUserOne(obj, L, N, M)

            im = sqrt(-1);         % imaginary
            FLDS_O = zeros(L,M);

            for m = 1 : floor(M/2)
                c1 = randn + im*randn;
                c2 = randn + im*randn;
                FLDS_O([1 2]+(m-1)*2, [1 2]+(m-1)*2) = [[c1 c2].' [c2' -c1'].'];
                %     c1 = randn + im*randn;
                %     c2 = randn + im*randn;
                %     FLDS_O([1 2]+(m)*2, [1 2]+(m)*2) = [[c1 c2].' [c2' -c1'].'];
            end
            C= FLDS_O;
        end

        % Generate LDS
        function [S] = gen_lds_array_v3(obj, N,K,d_v)

            S = zeros(N,K);
            n_el = d_v*K;

            for idx = 0:(n_el - 1)
                S(mod(idx,N)+1,mod(idx,K)+1) = exp(1i*2*pi*idx/K);
            end

            if K == 8
                S = S(:,[1,2,3,4,8,7,6,5]);
            end

            my_corr = zeros(K);
            for ii = 1:K
                for jj = 1:K
                    my_corr(ii,jj) = sum(real(S(:,ii).*conj(S(:,jj))));
                end
            end
            [~,my_corr_rank] = sort(my_corr);
            %
            % % Use pskmod to optimize graygen coding
            % symgray = pskmod(0:K-1,K,0,'gray');
            %
            % [~,psk_corr_rank] = sort(real(symgray'*symgray));


        end

        % read .bin file
        function [data] = readFile(obj, filename)

            if (0)
                s = dir(filename);
                fsize = s.bytes;
                fID = fopen(filename,'r');
                data = fread(fID, fsize/obj.dSizeAdj,'float');
                fclose(fID);
                data = data(1:2:end) + 1i*data(2:2:end);
            end

            skip = 4;
            fidr=fopen(filename); fidi=fopen(filename);  % handle for real, complex parts
            fseek(fidi,skip,'bof');                      % position file pointer for first complex
            data = complex(fread(fidr,inf,'single',skip), fread(fidi,inf,'single',skip));
            fclose all; clear fidr fidi
        end

        % write to .bin file
        function writeToFile(obj, filename, x)

            if (0)
                fID = fopen(filename,'w');
                fwrite(fID,x,'float');
                fclose(fID);
            end
            L = length(x);
            x_real = real(x);
            x_imag = imag(x);
            z(1:2:2*L) =  single(x_real);
            z(2:2:2*L) = single(x_imag);
            fidw=fopen(filename,'w' );
            fwrite(fidw,z,'single');
            fclose all; clear fidw

        end
    end
end
