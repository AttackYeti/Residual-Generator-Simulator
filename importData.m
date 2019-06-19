function [DATA] = importData(filename)

    % Import data
    DATA = table();
    % Detect table import options
    Topts = detectImportOptions(filename);
    
    % Read data file into MATLAB table format
    T = readtable(filename,Topts);
    %T = table2struct(T);
    % Extract data

    % Read thermocouples
    DATA.CX101   = T.CX_10_1;
    DATA.CX102   = T.CX_10_2;
    DATA.CX103   = T.CX_10_3;
    DATA.CX104   = T.CX_10_4;
    DATA.CX111   = T.CX_11_1;
    
    DATA.CX112   = T.CX_11_2;
    DATA.CX113   = T.CX_11_3;
    DATA.CX114   = T.CX_11_4;
    DATA.WT10    = T.WT_10;
    DATA.BT11    = T.BT_11;
    
    DATA.BT12    = T.BT_12;
    DATA.WT13    = T.WT_13;
    DATA.ST14E   = T.ST_14_E;
    DATA.ST14W   = T.ST_14_W;
    DATA.ST14N   = T.ST_14_N;
    
    DATA.WT40    = T.WT_40;
    DATA.BT41    = T.BT_41;
    DATA.WT42    = T.WT_42;
    DATA.BT43    = T.BT_43;
    DATA.AT01    = T.AT_01;
    
    DATA.AT02    = T.AT_02;
    DATA.WT62    = T.WT_62;
    DATA.BT63    = T.BT_63;
    DATA.WT64    = T.WT_64;
    DATA.BT65    = T.BT_65;

    % Read result Power Output data
    DATA.PowerOut = T.Power_Output;

    % Read desired power data
    DATA.DesiredPower = T.Desired_Power;

    % Read Spoof Enabled data
    %DATA.Spoof = T.SpoofEnabled;

    % Read primary loop flow data
    DATA.FM40 = T.FM_40;

    % Read CTAH Frequency data
    DATA.CTAHFreq = T.CTAH_Frequency;

    % % Read time stamps
    % TimeArray = T{:,find(DataName == 'Time')};
    % Time = strings(length(TimeArray),1);
    % for i =1:length(TimeArray)
    %     Time(i) = strcat('07-May-2019',{' '},convertCharsToStrings(TimeArray{i}));
    % end
    % Time = datetime(Time,'InputFormat','dd-MMM-yyyy h:mm:ss.SSS a');

    % Read system time
    DATA.Time = T.Time;

end