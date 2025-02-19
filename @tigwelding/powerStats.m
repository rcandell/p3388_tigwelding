function [...
    PeakDataPower_dBm, ...
    AvgDataPower_dBm, ...
    PAvgLeak_dBm, ...
    PAvgAnt_dBm, ...
    PMaxLeak_dBm, ...
    PMaxAnt_dBm] = powerStats(obj, pThreshold, titlestr) 
% POWERSTATS Compute the power statistics for the measurement object
%
% pThreshold: the dBm threshold for sample selection
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

% compute peak power
X = obj.cData;

pp = 10*log10((abs(X).^2)/obj.Rload);
pp = pp(pp>pThreshold);  % select powers above pThreshold
PMAX = max(pp);
ppl = pp(pp>(PMAX-15));  %includes low and high peaks
pph = pp(pp>(PMAX-6));   % includes only high peaks
PPEAK = mean(pph); % measured average peak power excludes noise floor
PAVG = mean(ppl); % measured avg power

% PPEAK = 10*log10(max(abs(X).^2)/obj.Rload); % measured peak power
% PAVG = 10*log10(mean(abs(X).^2)/obj.Rload); % measured avg power

% extract the meta data
FC = obj.meta_data_tbl(obj.meta_row_index, 'Fc');
FC = table2array(FC)*1e9;
DIST = table2array(obj.meta_data_tbl(obj.meta_row_index, 'Distance'));
GANT = table2array(obj.meta_data_tbl(obj.meta_row_index, 'AntennaGain'));
LCAB = table2array(obj.meta_data_tbl(obj.meta_row_index, 'CableLoss'));

% compute the power leakage at the arc
FSGAIN = tigwelding.freeSpaceGain(DIST,FC);
PAVGARC = PAVG - GANT + LCAB - FSGAIN;  % this is the avg power at the arc
PMAXARC = PPEAK - GANT + LCAB - FSGAIN;  % this is the peak power at the arc

% compute average power at 1 m, 10, m, and 100 m
% power is based on average power estimated at the arc
PAVGINT(1) = PAVGARC + tigwelding.freeSpaceGain(1,FC);
PAVGINT(2) = PAVGARC + tigwelding.freeSpaceGain(10,FC);
PAVGINT(3) = PAVGARC + tigwelding.freeSpaceGain(100,FC);

% compute peak power at 1 m, 10, m, and 100 m
% power is based on maximum power estimated at the arc
PMAXINT(1) = PMAXARC + tigwelding.freeSpaceGain(1,FC);
PMAXINT(2) = PMAXARC + tigwelding.freeSpaceGain(10,FC);
PMAXINT(3) = PMAXARC + tigwelding.freeSpaceGain(100,FC);

% plot average interference interference curve by distance
dPts = [1 10 100];
semilogx(dPts,PAVGINT);
xlabel('Distance (m)');
ylabel('Interference Power (dBm)')
grid on
hold on
semilogx(dPts,PMAXINT,'--');
xlabel('Distance (m)');
ylabel('Measured Interference Power (dBm)')
legend(['15 dB from max ' num2str(PMAX,3)], [' 6 dB from max ' num2str(PMAX,3)])
hold off 
title(titlestr);
drawnow

sroot = table2array(obj.meta_data_tbl(obj.meta_row_index, 'Directory'));
tigwelding.savePlotTo(gcf, obj.path_to_plots, sroot, 'intpwr');
close(gcf)

% return values
PeakDataPower_dBm = PPEAK;
AvgDataPower_dBm = PAVG;
PAvgLeak_dBm = PAVGARC;
PAvgAnt_dBm = PAVGINT;
PMaxLeak_dBm = PMAXARC;
PMaxAnt_dBm = PMAXINT;

end

