% op_addphase.m
% Jamie Near, McGill University 2014.
% Edits from
%   Edith Touchet-Valle, Texas A&M University 2024.
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% out=op_addphase(in,ph0,ph1,ppm0,suppressPlot);
% 
% DESCRIPTION:
% add zero and first order phase to a spectrum.
% 
% INPUTS:
% in             = input spectrum in matlab structure format
% ph0            = zero order phase to add (degrees)
% ph1            = 1st order phase to add (in seconds);
% ppm0           = (optional) frequency reference point.  Default = 4.65;
% suppressPlot   = (optional) Boolean to suppress plots.  Default = 0;
%
% OUTPUTS:
% out            = Phase adjusted output spectrum.


function out=op_addphase(in,ph0,ph1,ppm0,suppressPlot);

if nargin<5
    suppressPlot=1; % 01/05/23 modified by EV to 1
    if nargin<4
        gamma = (in.txfrq/1e6)/in.Bo;
        if gamma > 42 % For proton
            ppm0 = 4.65;
        else % For all other nuclei
            ppm0 = 0;
        end
        if nargin<3
            ph1=0;
        end
    end
end

%Add Zero-order phase
fids=in.fids.*ones(size(in.fids))*exp(1i*ph0*pi/180);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%Now add 1st-order phase
specs=addphase1(specs,in.ppm,ph1,in.Bo,in.txfrq,ppm0);

%re-calculate Fids using fft
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(size(specs,in.dims.t),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids=fft(fftshift(specs,in.dims.t),[],in.dims.t);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids=fft(circshift(fftshift(specs,in.dims.t),1),[],in.dims.t);
end

if length(size(specs))<3 && ~suppressPlot
    plot(in.ppm,real(specs));
    set(gca,'XDir','reverse');
end

    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;