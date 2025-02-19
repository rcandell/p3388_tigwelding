% ANALYZE_MEASUREMENTS top level script to process measurements
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
close all;
fclose all;
clc;
makeTSPlots = true;

meta_path = './MeasurementParameters.xlsx';
measdir = 'data';
fext = '.float32';

% output files
path_to_plots = 'figs';
fidResults = nan;

% write header to power stats file
fidResults = fopen("tigwelding_power_results.txt",'w');
fprintf(fidResults, 'Weld Type\t Base Metal\t Gas\t Distance (m)\t Fc (GHz)\t HF (ON/OFF)\t Current (A)\t Int Pwr Arc (dBm)\t input file name\n');
fprintf(            'Weld Type\t Base Metal\t Gas\t Distance (m)\t Fc (GHz)\t HF (ON/OFF)\t Current (A)\t Int Pwr Arc (dBm)\t input file name\n');

% open the meta data file
meta_data_tbl = tigwelding.importMeta(meta_path);

% loop through each row of the meta table
for jj = 1:height(meta_data_tbl)

    status = char(table2array(meta_data_tbl(jj, 'HF')));
    if strcmp(status, 'OFF')
        continue
    end

    % get title info
    weldtype = char(table2array(meta_data_tbl(jj, 'WeldType')));
    setnum = mat2str(table2array(meta_data_tbl(jj, 'Number')));
    runnum = mat2str(table2array(meta_data_tbl(jj, 'Run')));
    titlestr = [weldtype ' Set ' setnum ' Run ' runnum];

    % list files within specified directory
    dpath = char(table2array(meta_data_tbl(jj, 'Directory')));
    flist = dir([measdir '/' dpath '/*' fext]);
    if isempty(flist)
        continue
    end
    fpath = [flist(1).folder '/' flist(1).name];

    % construct an analysis object
    A = tigwelding(meta_data_tbl, jj, fpath, path_to_plots);
    A.loadCData();
    %A.resample(8,25); % resample to 100 Msps

    pThreshold = NaN;
    fc = table2array(meta_data_tbl(jj, 'Fc'));
    if fc < 1
        pThreshold = -48;
    else % fc >= 1
        pThreshold = -62;
    end

    % compute the power stats
    if true
        [PeakPower_dBm, AvgPower_dBm, Pleak_dBm, PInt_dBm, ...
            PMaxArc_dBm, PMaxInt_dBm] = A.powerStats(pThreshold, titlestr);
        % fidResults = fopen("tigwelding_power_results.txt",'a+');
        if jj~=1 
            fseek(fidResults, 0, 'eof'); 
        end
        fprintf(fidResults,'%-s\t %s\t %s\t %3.2f\t %3.2f\t %s\t %3.2f\t %3.2f\t %s\n', ...
            table2array(meta_data_tbl(jj,'WeldType')),  ...
            table2array(meta_data_tbl(jj,'BaseMetal')),  ...
            table2array(meta_data_tbl(jj,'GasUsed')),  ...
            table2array(meta_data_tbl(jj,'Distance')),  ...
            table2array(meta_data_tbl(jj,'Fc')),  ...
            table2array(meta_data_tbl(jj,'HF')),  ...
            table2array(meta_data_tbl(jj,'Current')),  ...
            PMaxArc_dBm,...
            table2array(meta_data_tbl(jj,'FileName')));

        fprintf('%-s\t %s\t %s\t %3.2f\t %3.2f\t %s\t %3.2f\t %3.2f\t %s\n', ...
            table2array(meta_data_tbl(jj,'WeldType')),  ...
            table2array(meta_data_tbl(jj,'BaseMetal')),  ...
            table2array(meta_data_tbl(jj,'GasUsed')),  ...
            table2array(meta_data_tbl(jj,'Distance')),  ...
            table2array(meta_data_tbl(jj,'Fc')),  ...
            table2array(meta_data_tbl(jj,'HF')),  ...
            table2array(meta_data_tbl(jj,'Current')),  ...
            PMaxArc_dBm, ...
            table2array(meta_data_tbl(jj,'FileName')));
    end

    if makeTSPlots
        A.plotTimeSeries(titlestr);
    end

    % save the workspace
    sPathToWs = strcat(measdir, '/', dpath, '/workspace.mat');
    save(sPathToWs);

end

fclose(fidResults);
fclose all;






