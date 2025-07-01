% GENERATE_TESTVECTORS_IEEE3388 top level script to make test vectors for IEEE 3388
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

clear;
close all;
fclose all;
clc;

meta_path = './MeasurementParameters.xlsx';
measdir = 'data';

% PSTATS output table
pstats_tbl_colnames = { ...
    'TVID','Num', 'Run', 'Fc_GHz', 'PowerIgnition', ...
    'HFOn', 'PulserOn', 'CurrentAmps', 'Fs', 'TVFilename', 'OriginalFilename' };
catalog_tbl = table('Size', [0, 11], ...
    'VariableNames', pstats_tbl_colnames, ...
    'VariableTypes', { 'double', 'double', 'double', 'double', 'string', ...
                       'string', 'string', 'double', 'double', 'string', 'string'});

% output files
testvec_path = 'testvec3388';

% open the meta data file
meta_data_tbl = tigwelding.importMeta(meta_path);

% loop through each row of the meta table
ii = 1;
for jj = 1:height(meta_data_tbl)

    num_ = meta_data_tbl{jj,"Number"};
    run_ = meta_data_tbl{jj,"Run"};
    Fc = table2array(meta_data_tbl(jj, 'Fc'));
    HF = char(table2array(meta_data_tbl(jj, 'HF')));
    OriginalFilename = char(table2array(meta_data_tbl(jj, 'FileName')));
    if ~strcmp(HF, 'ON')
        continue;
    end

    % construct a test vector object
    TV = tigwelding_testvector(meta_data_tbl, testvec_path);

    % % create test vector
    iq = TV.make_testvector(jj);
    if iq == -1
        disp("Skipping " + jj + " no data");
        continue
    end

    % write to file
    tv_filename = testvec_path + "/tv" + ii + ".csv";
    disp("Writing test vector " + ii + " Fc " + Fc + "GHz");
    writematrix(iq, tv_filename);

    % create the catalog
    piopt = meta_data_tbl{jj,"PowerAndIgnitionOption"};
    pulseron = meta_data_tbl{jj,"Pulser"};
    current_amps = meta_data_tbl{jj,"Current"};
    catalog = { ii, num_, run_, Fc, piopt, HF, pulseron, current_amps, TV.SampleRate, tv_filename, OriginalFilename };
    catalog_tbl(height(catalog_tbl)+1,:) = catalog;  

    t = 1000*(0:length(iq)-1)/TV.SampleRate;
    subplot(2,1,1), plot(t, real(iq)), ylabel('real(V)'), xlabel('t (ms)')
    subplot(2,1,2), plot(t, imag(iq)), ylabel('imag(V)'), xlabel('t (ms)')
    fnc = char(tv_filename);
    fnc = string([fnc(1:end-3) 'svg']);
    saveas(gcf, fnc);

    % increment main counter
    ii = ii + 1;

end

% write the power stats
writetable(catalog_tbl, testvec_path + "/catalog.xlsx");

fclose all;
close all






