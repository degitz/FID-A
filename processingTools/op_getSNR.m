% op_getSNR.m
% Jamie Near, McGill University 2014.
% Edits from
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% [SNR]=op_getSNR(in,signalppmmin,signalppmmax,noiseppmmin,noiseppmmax,suppressPlots);
% 
% DESCRIPTION:
% Find the SNR of the main peak in a spectrum. This is typically NAA for 1H.
% 
% INPUTS:
% in             = input data in matlab structure format
% signalppmmin      = min of frequncy range in which to search for main peak.
%                  (Optional.  Default = 1.8 ppm);
% signalppmmax      = max of frequncy range in which to search for main peak.
%                  (Optional.  Default = 2.2 ppm);
% noiseppmmin    = min of frequency range in which to measure noise.
%                  (Optional.  Default = -2 ppm);
% noiseppmmax    = max of frequency range in which to measure noise.
%                  (Optional.  Default = 0 ppm);
% suppressPlots  = Boolean to suppress plotting results.
%                  (Optional.  Default = false);
%
% OUTPUTS:
% SNR            = Estimated SNR of the input spectrum.
% signal         = The measured raw signal intensity
% noisesd        = The measured noise standard deviation

function [SNR,signal,noisesd]=op_getSNR(in,signalppmmin,signalppmmax,noiseppmmin,noiseppmmax,suppressPlots);

if nargin<6
    suppressPlots=false;
    if nargin<5
        noiseppmmax=0;
        if nargin<4
            noiseppmmin=-2;
            if nargin<3
                signalppmmax=2.2;
                if nargin<2
                    signalppmmin=1.8;
                end
            end
        end
    end
end

% Average data if not already done
if ~in.flags.averaged
    in = op_averaging(in);
end

%FIRST FIND THE NAA SIGNAL INTENSITY.  USE THE MAX PEAK HEIGHT OF THE 
%MAGNITUDE SPECTRUM INSIDE THE DESIRED SPECTRAL RANGE:
signalwindow=in.specs(in.ppm>signalppmmin & in.ppm<signalppmmax);
ppmwindow=in.ppm(in.ppm>signalppmmin & in.ppm<signalppmmax);

maxsignal_index=find(abs(signalwindow)==max(abs((signalwindow))));
maxsignal=abs(signalwindow(maxsignal_index));

if ~suppressPlots
    figure;
    plot(ppmwindow,abs(real(signalwindow)));

    figure
    plot(in.ppm,real(in.specs));
    if nargin < 5
        noiseppmmax=input('input upper ppm limit for noise: ');
        if nargin < 4
            noiseppmmin=input('input lower ppm limit for noise: ');
        end
    end
end

%NOW FIND THE STANDARD DEVIATION OF THE NOISE:
noisewindow=in.specs(in.ppm>noiseppmmin & in.ppm<noiseppmmax);
ppmwindow2=in.ppm(in.ppm>noiseppmmin & in.ppm<noiseppmmax)';

P=polyfit(ppmwindow2,noisewindow,2);

% Calculate SNR
noise=noisewindow-polyval(P,ppmwindow2);
noisesd=std(real(noise));
signal=(maxsignal-mean(real(noisewindow))); %Removes DC offset
%SNR=maxsignal/noisesd
SNR=signal/noisesd;

if ~suppressPlots
    figure;
    plot(ppmwindow,abs(real(signalwindow)));

    figure
    plot(in.ppm,real(in.specs));
    figure
    plot(ppmwindow2,real(noisewindow),...
        ppmwindow2,real(polyval(P,ppmwindow2)),...
        ppmwindow2,real(noise));
    disp(['The calculated signal-to-noise ratio is:  ' num2str(SNR) '.']);
end