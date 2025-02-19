function obj = show_pspectrum(obj, plotType, doSave)
% SHOW_PSPECTRUM compute, display, and save power spectrum
%
% Inputs:
%
%   plotType: the type of plot to produce 'contour' or 'waterfall'
%   doSave:  TRUE to save plot to file
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

    if strcmp(plotType, 'contour')

        pspectrum(obj.cData, obj.Fs, 'spectrogram', 'Leakage',0.5, ...
            'FrequencyLimits', [-obj.FreqLimit obj.FreqLimit]);
        drawnow
        
        if doSave
            sroot = table2array(obj.meta_data_tbl(obj.meta_row_index, 'Directory'));
            tigwelding.savePlotTo(gcf, obj.path_to_plots, sroot, 'contour');
        end

    elseif strcmp(plotType, 'waterfall')

        [p,f,t] = pspectrum(obj.cData, obj.Fs, 'spectrogram');
        waterfall(f/1e6, t*1000, p')
        xlim([-obj.FreqLimit obj.FreqLimit]/1e6)
        xlabel('Frequency (MHz)')
        ylabel('Time (ms)')
        wtf = gca;
        wtf.XDir = 'reverse';
        view([30 35])
        drawnow

        if doSave
            sroot = table2array(obj.meta_data_tbl(obj.meta_row_index, 'Directory'));
            tigwelding.savePlotTo(gcf, obj.path_to_plots, sroot, 'waterfall');
        end

    else

        pspectrum(obj.cData, obj.Fs, 'spectrogram', ...
            'Leakage',0.5);
        drawnow

    end

end