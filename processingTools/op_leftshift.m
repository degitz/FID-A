% op_leftshift.m
% Jamie Near, McGill University 2014.
% Edits from
%   Edith Touchet-Valle, Texas A&M University, 2024.
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% out=op_leftshift(in,ls);
% 
% DESCRIPTION:
% Remove leading datapoints from the fid to get rid of 1st order phase
% errors.
% 
% INPUTS:
% in     = input data in matlab strucure format.
% ls     = number of points to remove from the beginning of the fid.
%
% OUTPUTS:
% out    = Output following left shifting.  

function out=op_leftshift(in,ls);

if in.flags.leftshifted
    cont=input('WARNING:  Left shifting has already been performed!  Continue anyway?  (y or n)','s');
    if cont=='y'
        %continue;
    else
        error('STOPPING');
    end
end

fids=in.fids;

sz=in.sz;
sz(1)=sz(1)-ls;
fids=reshape(fids(ls+1:end,:),sz);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%Calculate t and ppm arrays using the calculated parameters:
f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):in.spectralwidth/(sz(1)):(in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];
ppm = -f/(in.txfrq/1e6); % to account for multiple nuclei - ETV 01/16/24
gamma = (in.txfrq/1e6)/in.Bo;
if gamma > 42 % For proton
    ppm = ppm + 4.65;
end

%t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];
t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];

    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.leftshifted=1;
