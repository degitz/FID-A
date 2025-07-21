% op_getLW.m
% Jamie Near, McGill University 2014.
%Edits from
%   Jacob Degitz, Texas A&M University 2025.
% 
% USAGE:
% [FWHM]=op_getLW(in,Refppmmin,Refppmmax,zpfactor,suppressPlots);
% 
% DESCRIPTION:
% Estimates the linewidth of a reference peak in the spectrum.  By default, 
% the reference peak is either: between 4.4 and 5.0 ppm for 1H, or between
% -0.5 and 0.5 ppm for all other nuclei.  Two methods are used to estimate
% the linewidth:  1.  FWHM is measured by simply taking the full width at
% half max of the reference peak.  2.  The FWHM is measured by fitting the
% reference peak to a lorentzian lineshape and determine the FWHM of the
% best fit.  The output FWHM is given by the average of these two measures.
% 
% INPUTS:
% in            = input spectrum in structure format.
% Refppmmin     = Min of frequency range (ppm) in which to search for reference peak.
%                  (Optional.  Default = 4.4 or -0.5 ppm, depending on gamma);
% Refppmmax     = Max of frequency range to (ppm) in which search for reference peak
%                  (Optional.  Default = 5.0 or 0.5 ppm, depending on gamma);
% zpfactor      = zero-padding factor (used for method 1.)
%                  (Optional.  Default = 8);
% suppressPlots = Boolean to suppress plotting results. 
%                  (Optional.  Default = false)
%
% OUTPUTS:
% FWHM       = Estimated linewidth of the input spectrum (in Hz).


function [FWHM]=op_getLW(in,Refppmmin,Refppmmax,zpfactor,suppressPlots);

if nargin<5
    suppressPlots=false;
    if nargin<4
        zpfactor=8;
        if nargin<3
            Refppmmax=5.0;
            if nargin<2
                Refppmmin=4.4;
            end
        end
    end
end

%in=io_readlcmraw(filestring,'dat');
in=op_zeropad(in,zpfactor);

%FIRST FIND FWHM USING TWO METHODS:

%METHOD 1:  ACTUALLY MEAUSURE FWHM OF WATER PEAK
Refwindow=in.specs(in.ppm>Refppmmin & in.ppm<Refppmmax);
ppmwindow=in.ppm(in.ppm>Refppmmin & in.ppm<Refppmmax);

maxRef_index=find(abs(real(Refwindow))==max(abs(real((Refwindow)))));
maxRef=real(Refwindow(maxRef_index));

if ~suppressPlots
    plot(ppmwindow,abs(real(Refwindow)),'.');
end

gtHalfMax=find(abs(real(Refwindow)) >= 0.5*abs(maxRef));

FWHM1=ppmwindow(gtHalfMax(1)) - ppmwindow(gtHalfMax(end));
FWHM1=FWHM1*(in.gamma*in.Bo);


%METHOD 2:  FIT WATER PEAK TO DETERMINE FWHM PARAM
sat='n';
centFreq=ppmwindow(maxRef_index);
while sat=='n'
    parsGuess=zeros(1,5);
    parsGuess(1)=maxRef; %AMPLITUDE
    parsGuess(2)=(5*in.Bo/3)/(in.gamma*in.Bo); %FWHM. LW = 5/3 Hz/T.
    parsGuess(3)=centFreq; %FREQUENCY
    parsGuess(4)=0; %Baseline Offset
    parsGuess(5)=0; %Phase
    
    yGuess=op_lorentz(parsGuess,ppmwindow);
    parsFit=nlinfit(ppmwindow,real(Refwindow'),@op_lorentz,parsGuess);
    yFit=op_lorentz(parsFit,ppmwindow);
    
    if ~suppressPlots
        figure;
        plot(ppmwindow,real(Refwindow),'.',ppmwindow,real(yGuess),':',ppmwindow,yFit);
        legend('data','guess','fit');
    
        sat=input('are you satisfied with fit? y/n [y] ','s');
        if isempty(sat)
            sat='y';
        end
        if sat=='n';
            centFreq=input('input new center frequency guess: ');
        end
    else
        sat='y';
    end

end


FWHM2=abs(parsFit(2));
FWHM2=FWHM2*(in.gamma*in.Bo);

FWHM=mean([FWHM1 FWHM2]);  

if ~suppressPlots
    disp(['The calculated linewidth is:  ' num2str(FWHM) ' Hz.' ]);
end

