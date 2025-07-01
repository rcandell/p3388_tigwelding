% ANALYZE_MEASUREMENTS top level script to process measurements
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
close all;
fclose all;
clc;

meta_path = './MeasurementParameters.xlsx';

% output files
path_to_plots = 'figs';
regions_dir = 'regions';
fidResults = nan;

% write header to power stats file
fidResults = fopen("tigwelding_power_results.txt",'w');
fprintf(fidResults, 'Weld Type\t Base Metal\t Gas\t Distance (m)\t Fc (GHz)\t HF (ON/OFF)\t Current (A)\t Intf Pwr Arc (dBm)\t input file name\n');
fprintf(            'Weld Type\t Base Metal\t Gas\t Distance (m)\t Fc (GHz)\t HF (ON/OFF)\t Current (A)\t Intf Pwr Arc (dBm)\t input file name\n');

% open the meta data file
meta_data_tbl = tigwelding.importMeta(meta_path);

% what to enable
BASELINE_ON = false; % when TRUE turns off everything else
PSTATS_ON = true;
PMAP_ON = true;
PREGIONS_ON = true;
DBSCAN_ON = true;
SAVEWS_ON = true;
TSPLOT_ON = true;

% scaling factors
tscale = 1;
fscale = 1e6;

% loop through each row of the meta table
for jj = 1:height(meta_data_tbl)

    status = char(table2array(meta_data_tbl(jj, 'HF')));
    weldtype = char(table2array(meta_data_tbl(jj, 'WeldType')));
    if BASELINE_ON
        if ~strcmp(weldtype, 'BASELINE')
            continue
        end        
    else
        if ~strcmp(status, 'ON')
            continue
        end
    end

    % get title info
    setnum = mat2str(table2array(meta_data_tbl(jj, 'Number')));
    runnum = mat2str(table2array(meta_data_tbl(jj, 'Run')));
    centerfreq = mat2str(table2array(meta_data_tbl(jj, 'Fc')));
    titlestr = [weldtype ' Set ' setnum ' Run ' runnum ' Fc ' centerfreq];

    % construct an analysis object
    A = tigwelding(meta_data_tbl, path_to_plots);

    A.loadCData(jj);
    meas_name = table2array(A.meta_data_tbl(jj, 'Directory'));
    %A.resample(8,25); % resample to 100 Msps

    pThreshold = NaN;
    fc = table2array(meta_data_tbl(jj, 'Fc'));
    if fc < 1
        pThreshold = -48;
    else % fc >= 1
        pThreshold = -62;
    end

    % compute the power stats
    if PSTATS_ON
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

    if TSPLOT_ON
        A.plotTimeSeries(titlestr);
    end    

    % Power map
    if PMAP_ON
        A.computePowerSpectrum();

        figure()
        A.plotIntensityMap(tscale,fscale);
        title("Spectrogram for " + titlestr) 
        [~, pngPathFull] = tigwelding.savePlotTo(gcf, path_to_plots, meas_name, 'pmap'); 

        % figure()
        % A.plotIntensityWaterfall(tscale,fscale);
        % title("Waterfall for " + titlestr) 
        % [~, pngPathFull] = tigwelding.savePlotTo(gcf, path_to_plots, meas_name, 'waterfall');         
        % % R.addSubSubSection('Power Spectrum Vs. Time');
        % % R.addPngFigure(latexreport.FIG_FLOAT, pngPathFull, 'Power intensity in dBm versus time');        
    end
    
    if PREGIONS_ON                      
        A.computePowerSpectrum();
        preg_thresh = mean(tigwelding.UtilizationThresholds);
        thresh_str = string(preg_thresh);

        if 0
        figure()
        A.computeBandwidthDurationsAboveThresholds(preg_thresh, 1, 1, tscale, fscale);
        title('All Regions for ' + meas_name + ' above thresh ' + thresh_str + ' dBm');
        [~, pngPathFull] = tigwelding.savePlotTo(gcf, path_to_plots, meas_name, 'allregions');
        % R.addSubSubSection('Clustering by Interference Regions');
        % R.addPngFigure(latexreport.FIG_FLOAT, pngPathFull, ...
        %     'All interference regions found without deselection applied.');

        % write to a text file
        regions_fname = regions_dir + "/" + meas_name + " - allregs.csv";
        fidRegions = fopen(regions_fname,'w');
        fprintf(fidRegions, 'Start Freq (MHz), Start Time (s), Duration, Bandwidth (MHz)'); 
        fclose(fidRegions);
        writetable(A.connectedRegionBoxesScaled, regions_fname, 'WriteMode','append')       

        end

        figure()
        A.computeBandwidthDurationsAboveThresholds(preg_thresh, 25, 0.999, tscale, fscale);
        title('Large Regions for ' + meas_name + ' above thresh ' + thresh_str + ' dBm');
        [~, pngPathFull] = tigwelding.savePlotTo(gcf, path_to_plots, meas_name, 'filteredregions');
        % R.addPngFigure(latexreport.FIG_FLOAT, pngPathFull, 'Large interference regions only.');

        % write to a text file
        regions_fname = regions_dir + "/" + meas_name + " - lgregs.csv";
        fidLgRegions = fopen(regions_fname,'w');
        fprintf(fidLgRegions, 'Start Time (s), Start Freq (MHz), Duration, Bandwidth (MHz)'); 
        fclose(fidLgRegions);
        if ~isempty(A.connectedRegionBoxesScaled)
            writetable(A.connectedRegionBoxesScaled, regions_fname, 'WriteMode','append')            
        end
    end

    if DBSCAN_ON
        A.computePowerSpectrum();
        clusterThreshold = mean(tigwelding.UtilizationThresholds);

        figure()
        A.computeClustersDBScan(clusterThreshold, 2, 32, fscale);
        title('DBSCAN Clusters for ' + meas_name);    
        [~, pngPathFull] = tigwelding.savePlotTo(gcf, path_to_plots, meas_name, 'dbscan');  
    end

    if SAVEWS_ON 
        % save the workspace
        sPathToWs = strcat(A.measdir, '/', A.dpath, '/workspace.mat');
        save(sPathToWs);        
    end    

    close all

end

fclose(fidResults);
fclose all;






