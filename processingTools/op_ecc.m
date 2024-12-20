% op_ecc.m
% Jamie Near, McGill University 2014.
% Edits from
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% [out,outw]=op_ecc(in,inw);
% 
% DESCRIPTION:
% Perform an eddy current correction by estimating any non-linearity in the
% phase of the water unsuppressed data in the time domain and applying the
% appropriate correction to both the water suppressed and water unsuppressed 
% data.
% 
% INPUTS:
% in     = water suppressed input data in matlab structure format.
% inw    = water unsuppressed input data in matlab structure format.
% caFLAG = OPTIONAL. Boolean flag that signifies averages are not combined - JND 9/2/2024
%
% OUTPUTS:
% out    = Water suppressed output following eddy current correction  
% outw   = Water unsuppressed output following eddy current correction

function [out,outw]=op_ecc(in,inw,caFLAG);

if inw.dims.coils~=0 || inw.dims.subSpecs~=0
    error('ERROR:  Must combine receivers and subspecs prior to running ecc!! ABORTING!!');
end
if inw.dims.averages~=0
    if nargin < 3
        error('ERROR:  Must combine averages prior to running ecc!! ABORTING!!');
    end
elseif ~caFLAG
    error('ERROR:  Must combine averages prior to running ecc!! ABORTING!!');
else
    warning('WARNING: Averages will be temporarily combined to perform ecc.')
end

% Temporarily combine averages - JND 9/2/2024
if ~caFLAG
    in_ca = op_averaging(in);
    inw_ca = op_averaging(inw);
else
    in_ca = in;
    inw_ca = inw;
end

%save the phase as a vector of hard numbers.
inph=double(phase(inw_ca.fids));
figure; plot(inw_ca.t,inph);

%choose the part of the phase function that is most linear
tmin=input('input min t value: ');
tmax=input('input max t value: ');
figure;

%now fit a straight line to the linear part of the phase function
p=polyfit(inw_ca.t(inw_ca.t>tmin & inw_ca.t<tmax), inph(inw_ca.t>tmin & inw_ca.t<tmax)',1);

%now fit a spline to approximate a smooth version of the phase function
pp=splinefit(inw_ca.t,inph,150);

%now subtract the line from the spline to get the eddy current related
%phase offset:
ecphase=ppval(pp,inw_ca.t)'-polyval(p,inw_ca.t)';
sz=size(in_ca.fids);
ecphase_rep=repmat(ecphase,[1 sz(2:end)]);
figure;
plot(inw_ca.t,ecphase);


%Now apply the eddy current correction to both the water suppressed and the
%water unsuppressed data:
out=in;
out.fids=out.fids.*exp(1i*-ecphase_rep);
out.specs=fftshift(ifft(out.fids,[],out.dims.t),out.dims.t);
out=op_addphase(out,180*ecphase_rep(1)/pi);

outw=inw;
outw.fids=outw.fids.*exp(1i*-ecphase);
outw.specs=fftshift(ifft(outw.fids,[],out.dims.t),out.dims.t);
outw=op_addphase(outw,180*ecphase(1)/pi);

% Temporarily combine averages - JND 9/2/2024
if ~caFLAG
    out_ca = op_averaging(out);
    outw_ca = op_averaging(outw);
else
    out_ca = out;
    outw_ca = outw;
end

% Plot data
figure;
plot(out_ca.t,phase(out_ca.fids)); title('Phase of Water suppressed data') 
figure;
plot(outw_ca.t,phase(outw_ca.fids));title('Phase of Water unsuppressed data') 