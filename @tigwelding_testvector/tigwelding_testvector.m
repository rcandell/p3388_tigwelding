classdef tigwelding_testvector < tigwelding
    
    properties
        path_to_tv
        co_tbl
    end

    properties (Constant)
        SampleRate = 625e6; % Hz
        mag_lower_bound = sqrt(power(10,-60/10));  % dBm
        tv_pwr_cutoff = [ ...
            1, 1, -45; 1, 2, -45; 1, 3, -65; 
            2, 1, -55; 2, 2, -62; 
            3, 1, -60; 3, 2, -62;  
            4, 1, -60; 4, 2, -62; 
            5, 1, -55; 5, 2, -62; 
            6, 1, -55; 6, 2, -62; 
            11, 1, -55; 11, 2, -62             
            ];
        norm_pwr = false;
    end
    
    methods
        function obj = tigwelding_testvector(h_meta_tbl, path_to_tv)
            obj = obj@tigwelding(h_meta_tbl, nan);
            obj.path_to_tv = path_to_tv;
            obj.co_tbl = array2table(obj.tv_pwr_cutoff, 'VariableNames',{'Num','Run','Cutoff'});
        end       
        
        function iq = make_testvector(obj, jj)

            % get measurement number and run 
            num_ = obj.meta_data_tbl{jj,"Number"};
            run_ = obj.meta_data_tbl{jj,"Run"};

            % get the raw iq data
            obj.loadCData(jj);
            iq = obj.cData;

            % eliminate the noise floor
            co_tbl = obj.co_tbl; %#ok<*PROPLC>
            logicalIndex = co_tbl.Num == num_;
            co_tbl = co_tbl(logicalIndex, :);
            logicalIndex = co_tbl.Run == run_;
            if sum(logicalIndex) > 0
                co_tbl = co_tbl(logicalIndex, :);
                cutoff = table2array(co_tbl(1,"Cutoff"));
                iq(10*log10((abs(iq).^2))<cutoff) = 0;
            else
                iq = -1;
                return
            end

            % normalize the data to unity power
            if obj.norm_pwr
                p = mean(abs(iq).^2);
                iq = iq/sqrt(p);
            else
                a = max(abs(iq));
                iq = iq/a;
            end

        end
    end
end

