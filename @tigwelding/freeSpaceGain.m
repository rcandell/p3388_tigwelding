function gain_dB = freeSpaceGain(DIST_m, FC_Hz)
% FREESPACEGAIN Compute the free-space path loss per Frii's equations
%
% Outputs:
%   gain_dB the path gain in dB
%
% Inputs:
%   DIST_m distance in meters
%   FC_Hz center frequency in Hz
%
% Author: Rick Candell
% Organization: National Institute of Standards and Technology
% Email: rick.candell@nist.gov

gain = (tigwelding.C/(4*pi*DIST_m*FC_Hz));
gain_dB = 20*log10(gain);

end
