function [DAQ, SUP] = importData(DAQFileName, SUPERVISORFileName)    
    %% Import data
    DAQ = table();
    SUP = table();
    
    % Detect table import options
    TDAQopts = detectImportOptions(DAQFileName);
    TSUPopts = detectImportOptions(SUPERVISORFileName);
    
    % Read data file into MATLAB table format
    TDAQ = readtable(DAQFileName,TDAQopts);
    TSUP = readtable(SUPERVISORFileName, TSUPopts);

    %% Import DAQ Data
    DAQ.CTAHFreq        = TDAQ.CTAH_Frequency;
    DAQ.TCHXFreq        = TDAQ.TCHX_Frequency;
    DAQ.PumpFreq        = TDAQ.Pump_Frequency;
    DAQ.DesiredPower    = TDAQ.Desired_Power;
    DAQ.InputPowerArray = TDAQ.Power_Input_Signal;

    % Read result Power Output data
    DAQ.PowerOut = TDAQ.Power_Output;

    % Read system time 
    DAQ.Time = TDAQ.Time;

    DAQ.CX101   = TDAQ.CX_10_1;
    DAQ.CX102   = TDAQ.CX_10_2;
    DAQ.CX103   = TDAQ.CX_10_3;
    DAQ.CX104   = TDAQ.CX_10_4;
    DAQ.CX111   = TDAQ.CX_11_1;
    DAQ.CX112   = TDAQ.CX_11_2;
    DAQ.CX113   = TDAQ.CX_11_3;
    DAQ.CX114   = TDAQ.CX_11_4; 
    DAQ.WT10    = TDAQ.WT_10;
    DAQ.BT11    = TDAQ.BT_11;
    DAQ.BT12    = TDAQ.BT_12;
    DAQ.WT13    = TDAQ.WT_13;
    DAQ.ST14E   = TDAQ.ST_14_E;
    DAQ.ST14W   = TDAQ.ST_14_W;
    DAQ.ST14N   = TDAQ.ST_14_N;
    DAQ.WT40    = TDAQ.WT_40;
    DAQ.BT41    = TDAQ.BT_41;
    DAQ.WT42    = TDAQ.WT_42;
    DAQ.BT43    = TDAQ.BT_43;
    DAQ.AT01    = TDAQ.AT_01;
    DAQ.AT02    = TDAQ.AT_02;
    DAQ.WT62    = TDAQ.WT_62;
    DAQ.BT63    = TDAQ.BT_63;
    DAQ.WT64    = TDAQ.WT_64;
    DAQ.BT65    = TDAQ.BT_65;
    DAQ.FM40    = TDAQ.FM_40;
     
    %% Import Supervisor Data
    SUP.Time = TSUP.Time_R;
    SUP.WT10 = TSUP.WT_10_R;
    SUP.BT11    = TSUP.BT_11_R;
    SUP.CX101   = TSUP.CX_10_1_R;
    SUP.CX102   = TSUP.CX_10_2_R;
    SUP.CX103   = TSUP.CX_10_3_R;
    SUP.CX104   = TSUP.CX_10_4_R;
    SUP.CX111   = TSUP.CX_11_1_R;
    SUP.CX112   = TSUP.CX_11_2_R;
    SUP.CX113   = TSUP.CX_11_3_R;
    SUP.CX114   = TSUP.CX_11_4_R;
    SUP.BT12    = TSUP.BT_12_R;
    SUP.WT13    = TSUP.WT_13_R;
    SUP.WT20 = TSUP.WT_20_R;
    SUP.BT21 = TSUP.BT_21_R;
    SUP.WT22 = TSUP.WT_22_R;
    SUP.BT23 = TSUP.BT_23_R;
    SUP.WT24 = TSUP.WT_24_R;
    SUP.BT25 = TSUP.BT_25_R;
    SUP.WT26 = TSUP.WT_26_R;
    SUP.BT27 = TSUP.BT_27_R;
    SUP.WT28 = TSUP.WT_28_R;
    SUP.BT29 = TSUP.BT_29_R;
    SUP.BT30 = TSUP.BT_30_R;
    SUP.WT31 = TSUP.WT_31_R;
    SUP.BT32 = TSUP.BT_32_R;
    SUP.WT33 = TSUP.WT_33_R;
    SUP.BT34 = TSUP.BT_34_R;
    SUP.WT35 = TSUP.WT_35_R;
    SUP.WT40 = TSUP.WT_40_R;
    SUP.BT41 = TSUP.BT_41_R;
    SUP.WT42 = TSUP.WT_42_R;
    SUP.BT43 = TSUP.BT_43_R;
    SUP.BT60 = TSUP.BT_60_R;
    SUP.WT61 = TSUP.WT_61_R;
    SUP.WT62 = TSUP.WT_62_R;
    SUP.BT63 = TSUP.BT_63_R;
    SUP.WT64 = TSUP.WT_64_R;
    SUP.BT65 = TSUP.BT_65_R;
    SUP.BT66 = TSUP.BT_66_R;
    SUP.WT67 = TSUP.WT_67_R;
    SUP.ST10 = TSUP.ST_10_R;
    SUP.ST11 = TSUP.ST_11_R;
    SUP.ST12N   = TSUP.ST_12_N_R;
    SUP.ST12SE  = TSUP.ST_12_SE_R;
    SUP.ST12SW  = TSUP.ST_12_SW_R;
    SUP.ST13    = TSUP.ST_13_R;
    SUP.ST14N   = TSUP.ST_14_N_R;
    SUP.ST14E   = TSUP.ST_14_E_R;
    SUP.ST14S   = TSUP.ST_14_S_R;
    SUP.ST14W   = TSUP.ST_14_W_R;
    SUP.AT01    = TSUP.AT_01_R;
    SUP.AT02    = TSUP.AT_02_R;
    SUP.FM20    = TSUP.FM_20_R;
    SUP.FM30    = TSUP.FM_30_R;
    SUP.FM40    = TSUP.FM_40_R;
    SUP.FM60    = TSUP.FM_60_R;
    % In the future, SUP.Desired power will be taken from the instructor
    % log
    SUP.DesiredPower = ones(height(SUP),1).*8000;
    
    
    % Data below only used in analysis
    SUP.PowerOut = TSUP.PowerOutput_R;
    % Read Spoof Enabled data
    SUP.Spoof = TSUP.SpoofEnabled;

    SUP.PumpFrequency = TSUP.PumpFrequency_R;
    SUP.PowerInputSignal = TSUP.PowerInputSignal_R;
    SUP.CTAHFrequency = TSUP.CTAHFrequency_R;
    
    %% Special Formatting options
    SUP = SUP(1:32376,:);
    
    % Align the data in DAQ with the data in SUP
    for i = 2:height(DAQ)
        if DAQ.Time(i) == SUP.Time(1)
            if (i+height(SUP)) < height(DAQ)
                DAQ = DAQ(i:(i-1)+height(SUP), :);
            else
                disp("Second Portion");
                DAQ = DAQ(i:end, :);
                SUP = SUP(1:height(DAQ), :);
            end
            break
        end  
    end
    
end