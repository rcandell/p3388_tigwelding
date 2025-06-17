classdef tigwelding < handle
    %tigwelding Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % meta data spreadsheet
        meta_data_tbl = [];
        meta_row_index = [];

        % input file information
        path_to_cdata = '';
        dpath;

        % output file information
        path_to_plots = '';

        % complex data loaded from binary file
        cData = nan;

        % resampler filter
        Fs = 625e6;

        % scaling parameters
        target_power = 0;  % in dB, default to unit average power

        % power calculations over time and frequency
        tt = [];               % durations from beginning in seconds
        meas_pwr_dBm = [];      % measured powers
        ff = [];                % the frequencies        

        % for clustering 
        channelUtilization = []; % channel utilization above thresholds
        bandwidthDurations = []; % array of bandwidths and durations above thresholds
        connectedRegions = [];   % contains connected regions
        connectedRegionStats = [];  % contains bounded box information on regions
        connectedRegionBoxesScaled = []; % contains all of the spectral activity events
        dbscanClusters  = [];  % results from dbscan clustering        

    end

    properties (Constant)
        
        % file extension
        fext = '.float32';
        measdir = 'data';

        C = 299792458;  % speed of light m/s
        Fc = 2450e6;    % Hz
        FreqLimit = 50e6; % Hz
        Rload = 50; % Load in Ohms

        % Frequency resolution of the RTSA
        FrequencyResolution = 160e3;

        UtilizationThresholds = [-79 -85 -91];  % for utilization metrics based on backoff
        UtilizationThresholdsStrings = {'Above -79 dBm at 5 m';'Above -85 dBm at 5 m';'Above -91 dBm at 5 m'};

        % wifi channel map
        PsuedoWiFiChannelBounds = [ ...
            2412-10, 2412+10;
            2437-10, 2437+10;
            2462-10, 2462+10
            ];     

        % parameters for DBSCAN clustering
        maxDistanceDBSCAN_Hz = Inf; %20e6;
        dbscanEpsilon_dB = 1.5; % +/- 1.5 dB neighborhood        
    end


    methods(Static)

        microwave = importMeta(obj, workbookFile)
        L = freeSpaceGain(dist, fc)

        % function savePlotTo(hFig, sPathToFolder, sRootName, sPlotTypeName)
        %     fpath = strcat(sPathToFolder, '/', sRootName, '_', sPlotTypeName);
        %     % savefig(hFig, strcat(fpath, '.fig'));
        %     print(hFig, strcat(fpath, '.png'), '-dpng', '-r300');
        % end
        function [figPathFull, pngPathFull] = savePlotTo(hFig, sPathToFolder, sRootName, sPlotTypeName)
            fpath = strcat(sPathToFolder, '/', sRootName, ' - ', sPlotTypeName);
            figPathFull = strcat(fpath, '.fig');
            savefig(hFig, figPathFull);
            pngPathFull = strcat(fpath, '.png');
            print(hFig, pngPathFull, '-dpng', '-r300');
        end        

        function dist = dbscan_Distance(ZI, ZJ)

            % by power
            % assuming the out of power range already set to -Inf
            dist = abs(ZI(3) - ZJ(:,3));  

            % now by distance in frequency
            if ~isinf(tigwelding.maxDistanceDBSCAN_Hz)
                indsFar = abs(ZJ(:,1)-ZI(1)) > tigwelding.maxDistanceDBSCAN_Hz;
                dist(indsFar) = nan;
            end

        end        

    end

    methods
        function obj = tigwelding(h_meta_tbl, path_to_plots)
            obj.meta_data_tbl = h_meta_tbl;          
            obj.path_to_plots = path_to_plots;
        end

        function obj = loadCData(obj, meta_row_index)

            obj.meta_row_index = meta_row_index;

            % list files within specified directory
            obj.dpath = char(table2array(obj.meta_data_tbl(meta_row_index, 'Directory')));
            flist = dir([tigwelding.measdir '/' obj.dpath '/*' tigwelding.fext]);
            obj.path_to_cdata = [flist(1).folder '/' flist(1).name];    

            % open the file
            fid = fopen(obj.path_to_cdata,"r");
            x = fread(fid, 'single');
            x = complex(x(1:2:end-1),x(2:2:end));
            obj.cData = x;
            fclose(fid);
        end

        function obj = rescale(obj, target_power)
            % rescale to specified average power
            obj.target_power = target_power;
            linpwr = 10^(target_power/10);
            xpwr = mean(abs(obj.cData).^2);
            pscale = linpwr/xpwr;
            y = obj.cData*sqrt(pscale);
            % ypwr = 10*log10(mean(abs(y).^2));
            % fprintf('Power dBm is %f\n', ypwr);
            obj.cData = y;
        end

        function obj = writeCData(obj)

            % write real data 
            rData = real(obj.cData);
            fpath = strcat(obj.path_to_cdata,'.real');
            fid = fopen(fpath,'w');
            fwrite(fid,single(rData),'float32');
            fclose(fid);

            % write the imaginary data
            rData = imag(obj.cData);
            fpath = strcat(obj.path_to_cdata,'.imag');
            fid = fopen(fpath,'w');
            fwrite(fid,single(rData),'float32');
            fclose(fid);

        end

        function resample(obj, M, N)
            obj.cData = resample(obj.cData, M, N);
            obj.Fs = M/N*obj.Fs;
        end

        function scalePowerToReferenceDistance(obj)
            obj.computePowerSpectrum();
            measGain = intfmeas.freeSpaceGain(obj.distance, mean(obj.ff));
            refGain = intfmeas.freeSpaceGain(intfmeas.ReferenceDistance, mean(obj.ff));
            obj.meas_pwr_dBm = obj.meas_pwr_dBm + (refGain-measGain);
        end

        function computeUtilizationAboveThresholds(obj)
            utilization = zeros(length(obj.UtilizationThresholds),size(obj.meas_pwr_dBm,2));
            for utii = 1:length(obj.UtilizationThresholds)
                ut = obj.UtilizationThresholds(utii);
                derUtilization = obj.meas_pwr_dBm > ut;
                derUtilization = sum(derUtilization)/size(derUtilization,1);
                utilization(utii,:) = derUtilization;
            end
            obj.channelUtilization = utilization;
        end

        function plotIntensityMap(obj, t_scale, f_scale)
            obj.computePowerSpectrum();
            s=pcolor(obj.tt/t_scale,obj.ff/f_scale,obj.meas_pwr_dBm);
            ylabel('frequency (MHz)');
            xlabel('time (s)');
            set(s, 'EdgeColor', 'none');
            s.FaceColor = 'interp';
            colormap turbo
            colorbar
        end

        function plotIntensityWaterfall(obj, t_scale, f_scale)
            obj.computePowerSpectrum();
            [T,F] = meshgrid(obj.tt/t_scale, obj.ff/f_scale);
            waterfall(T, F, obj.meas_pwr_dBm);
            ylabel('frequency (MHz)');
            xlabel('time (s)');
            % set(s, 'EdgeColor', 'none');
            % s.FaceColor = 'interp';
            colormap turbo
            colorbar
        end        

        function computePowerSpectrum(obj) 
            if isempty(obj.meas_pwr_dBm)
                [p,f,t] = pspectrum(obj.cData, obj.Fs, 'spectrogram', 'OverlapPercent', 3);
                obj.meas_pwr_dBm = 10*log10(p);
                obj.ff = f;
                obj.tt = t;
            end
        end

        function plotPowerSpectrum(obj, f_scale)
            obj.computePowerSpectrum();
            % plot max and average power vs frequency
            plot(obj.ff/f_scale, [obj.max_meas_pwrs(:), obj.avg_meas_pwrs(:)]);
            xlabel('frequency (MHz)');
            ylabel('Power (dBm)');
            legend(['max';'avg'],'AutoUpdate','off')
        end

        function plotUtilizationAboveThresholds(obj, f_scale)
            % plot channel utilization above the thresholds
            plot(obj.ff/f_scale, obj.channelUtilization);
            xlabel('frequency (MHz)');
            ylabel('Pr. Channel Utilization');
            legend(obj.UtilizationThresholdsStrings,'AutoUpdate','off');
        end        

        function plotTimeSeries(obj, titlestr)
            x = 20*log10(abs(obj.cData));
            t = 1000*(0:length(x)-1)/obj.Fs;
            plot(t,x)
            xlabel('time (ms)')
            ylabel('Measured Interference Amplitude (dB-V)')
            title(titlestr);
            drawnow
            sroot = table2array(obj.meta_data_tbl(obj.meta_row_index, 'Directory'));
            tigwelding.savePlotTo(gcf, obj.path_to_plots, sroot, 'timeseries');
            close(gcf)
        end

        function computeBandwidthDurationsAboveThresholds(obj, pwrthresh, eventarea, solidity, t_scale, f_scale)

            obj.computePowerSpectrum();
            X = obj.meas_pwr_dBm;
            above = X > pwrthresh;
            imagesc(obj.tt/t_scale, obj.ff/f_scale, above);
            set(gca,'YDir','normal')
            ylabel('frequency (MHz)');
            xlabel('time (s)');            
            colormap gray

            obj.connectedRegions = bwconncomp(above, 8);
            obj.connectedRegionStats = regionprops(obj.connectedRegions, ...
                'BoundingBox','Area','Solidity');    

            % created table for properties of connected regions
            connectedRegionBoxesScaledColumnNames = {'Start Time (s)', ...
                'Start Freq (MHz)', 'Duration', 'Bandwidth (MHz)'};
            BoxesScaled = [];

            % overlay bounding boxes
            hold on;
            dt = obj.tt(end)/size(above,2);
            df = (obj.ff(end)-obj.ff(1))/size(above,1);
            ff0 = obj.ff(1);
            tt0 = obj.tt(1);
            
            % Loop through each connected component and plot the bounding boxes
            kk = 1;
            for ii = 1:length(obj.connectedRegionStats)
                % Get the bounding box for the current component
                bb = obj.connectedRegionStats(ii).BoundingBox;

                % ignore small and unfilled events
                % basically want to keep events connected in frequency and
                % time that are significant in impact
                if obj.connectedRegionStats(ii).Area < eventarea || ...
                    obj.connectedRegionStats(ii).Solidity > solidity
                    continue
                end

                % scale and bias each bounding box
                bb(1) = bb(1) - 0.5;
                bb(2) = bb(2) - 0.5;
                % bbs = [ ...
                %     (bb(1)*df+ff0)/f_scale, ...
                %     (bb(2)*dt+tt0)/t_scale, ...
                %     bb(3)*df/f_scale, ...
                %     bb(4)*dt/t_scale];              
                bbs = [ ...
                    (bb(1)*dt+tt0)/t_scale, ...
                    (bb(2)*df+ff0)/f_scale, ...
                    bb(3)*dt/t_scale, ...
                    bb(4)*df/f_scale];              

                BoxesScaled = [BoxesScaled; bbs];
                kk = kk + 1;
                
                % Plot the bounding box
                rectangle('Position', bbs, 'EdgeColor', 'r', 'LineWidth', 2);
            end
            hold off

            % Add data to table for connected regions
            if ~isempty(BoxesScaled)
                obj.connectedRegionBoxesScaled = array2table(BoxesScaled);
                obj.connectedRegionBoxesScaled.Properties.VariableNames = ...
                    connectedRegionBoxesScaledColumnNames;
            end
        end        

        function computeClustersDBScan(obj, pwrthresh, minpts, decrate, f_scale)
            powers = obj.meas_pwr_dBm;
            % imagesc(powers);

            if size(powers,1) > 800
                % Define a Gaussian smoothing filter
                % filter = fspecial('gaussian', [3, 3], 0.5);
                filter = fspecial('gaussian', decrate, 1/decrate);
                % filter = filter/sum(filter(:));
                
                % Smooth the matrix using 2D convolution
                powers2 = conv2(powers, filter, 'same');
                
                % Downsample the matrix (though it's already small, we'll use 'imresize')
                powers2 = imresize(powers2, 1/decrate, 'bilinear');

                % Downsample the freq and time coordinates
                ff2 = obj.ff(1:decrate:end, 1:decrate:end);
                tt2 = obj.tt(1:decrate:end, 1:decrate:end);
                
                imagesc(ff2/f_scale,tt2,powers2);
                set(gca,'YDir','normal')
                colorbar
            else
                powers2 = powers;
                % Downsample the freq and time coordinates
                ff2 = obj.ff;
                tt2 = obj.tt;
                imagesc(ff2/f_scale,tt2,powers2);
                set(gca,'YDir','normal')
                colorbar
            end

            % form the matrix used for dbscan
            [row, col] = size(powers2);
            [X, Y] = meshgrid(obj.ff(1:col), obj.tt(1:row));
            data = [X(:), Y(:), powers2(:)];

            % clean out out of range measurements
            rowsUnder = data(:,3)<pwrthresh;
            data(rowsUnder, 3) = -Inf;  % will cause nan distance

            % Apply DBSCAN clustering
            [idx, ~] = dbscan(data, tigwelding.dbscanEpsilon_dB, minpts, 'Distance', @tigwelding.dbscan_Distance);
            
            % Reshape the cluster indices back into the original matrix shape
            clusters = reshape(idx, row, col);
            obj.dbscanClusters = clusters;
            figure()
            imagesc(tt2, ff2/f_scale, clusters);
            set(gca,'YDir','normal')
            ylabel('Frequency (MHz)');
            xlabel('Time (s)');
            % grid on;         
        end

    end
end






