classdef tigwelding < handle
    %tigwelding Summary of this class goes here
    %   Detailed explanation goes here

    properties

        % meta data spreadsheet
        meta_data_tbl = [];
        meta_row_index = 0;

        % input file information
        path_to_cdata = '';

        % output file information
        path_to_plots = '';

        % complex data loaded from binary file
        cData = nan;

        % resampler filter
        Fs = 625e6;

        % scaling parameters
        target_power = 0;  % in dB, default to unit average power

    end

    properties (Constant)
        C = 299792458;  % speed of light m/s
        Fc = 2450e6;    % Hz
        FreqLimit = 50e6; % Hz
        Rload = 50; % Load in Ohms
    end


    methods(Static)

        microwave = importMeta(obj, workbookFile)
        L = freeSpaceGain(dist, fc)

        function savePlotTo(hFig, sPathToFolder, sRootName, sPlotTypeName)
            fpath = strcat(sPathToFolder, '/', sRootName, '_', sPlotTypeName);
            % savefig(hFig, strcat(fpath, '.fig'));
            print(hFig, strcat(fpath, '.png'), '-dpng', '-r300');
        end

    end

    methods
        function obj = tigwelding(h_meta_tbl, meta_row_index, path_to_cdata, path_to_plots)
            obj.meta_data_tbl = h_meta_tbl;
            obj.meta_row_index = meta_row_index;
            obj.path_to_cdata = path_to_cdata;
            obj.path_to_plots = path_to_plots;
        end

        function obj = loadCData(obj)
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

    end
end






