# Antijamming
# Matlab files
1. Place matlab files inside a folder (e.g., /space/code/ , etc.). The main file is Antijamming.m

# To run Antijamming.m as TX - generator
1. You need to specify your txPathDir path directory. If not specified, it will take 'pwd'. As an example, you can run matlab command as below.
2. Also, you need to specify runType   = "TX"; (e.g., Antijamming(txPathDir,[], [], [],'TX') )
3. matlab -nodisplay -nosplash -nodesktop -r "txPathDir='/space/TX'; rxPathDir=''; dataPathDir='';  snr=1; runType='TX'; run /space/code/Antijamming(txPathDir, rxPathDir, dataPathDir, runType);exit;"
4. The .bin files will be saved in txPathDir directory. 

# To run Antijamming.m as RX - receiver
1. You need to specify your txPathDir and rxPathDir path directories. If not specified, it will take 'pwd'. As an example, you can run matlab command as below.
2. Also, you need to specify runType   = "RX"; (e.g., Antijamming(txPathDir,rxPathDir, [], [],'RX') )
3. matlab -nodisplay -nosplash -nodesktop -r "txPathDir='/space/TX'; rxPathDir='/space/RX'; dataPathDir='';  snr=1; runType='RX'; run /space/code/Antijamming(txPathDir, rxPathDir, dataPathDir, runType);exit;"
4. It will read the .bin from the RX folder and will demodulate it and compare with the TX folder files. Error will be displayed.

# To run Antijamming.m as SIM_CHANNEL 
1. You need to specify your txPathDir and rxPathDir path directories. If not specified, it will take 'pwd'. As an example, you can run matlab command as below.
2. Also, you need to specify runType   = "SIM_CHANNEL"; (e.g., Antijamming(txPathDir,rxPathDir, [], snr, 'SIM_CHANNEL') )
3. matlab -nodisplay -nosplash -nodesktop -r "txPathDir='/space/TX'; rxPathDir='/space/RX'; dataPathDir=''; snr=1; runType='SIM_CHANNEL'; run /space/code/Antijamming(txPathDir, rxPathDir, dataPathDir, runType);exit;"
4. It will read the .bin from the TX folder and will add AWGN channel and save it to RX folder.


# To run Antijamming.m as a DATA VISUALIZATION. 
1. You need to specify your dataPathDir path directory. If not specified, it will take 'pwd'. As an example, you can run matlab command as below.
3. Also, you need to specify  runType   = "DATA_VISUALIZATION"; (e.g., Antijamming([],[], dataPathDir, [], 'DATA_VISUALIZATION') )
4. matlab -nodisplay -nosplash -nodesktop -r "txPathDir=''; rxPathDir=''; dataPathDir='/space/RX';  snr=1; runType='DATA_VISUALIZATION'; run /space/code/Antijamming(txPathDir, rxPathDir, dataPathDir, runType);exit;"
5. The collective spectrogram result will be saved in dataPathDir/visualize folder. 

# To run Antijamming.m as a SIMULATE. 
1. You donnot need to specify your any director. 
3. Also, you need to specify  runType   = "SIMULATE"; (e.g., Antijamming([],[], [], [], 'SIMULATE') )
4. matlab -nodisplay -nosplash -nodesktop -r "txPathDir=''; rxPathDir=''; dataPathDir='';  snr=1; runType='SIMULATE'; run /space/code/Antijamming(txPathDir, rxPathDir, dataPathDir, runType);exit;"
5. The BER and SER plots will be generated.
