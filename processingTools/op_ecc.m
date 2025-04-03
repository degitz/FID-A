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
%
% OUTPUTS:
% out    = Water suppressed output following eddy current correction  
% outw   = Water unsuppressed output following eddy current correction

function [out,outw]=op_ecc(in,inw,caFLAG);

if inw.dims.subSpecs~=0
    error('ERROR:  Must combine subspecs prior to running ecc!! ABORTING!!');
end

% Temporarily combine channels - JND 4/1/2025
if inw.dims.coils~=0
    warning('WARNING: Channels will be temporarily combined to perform ecc.')
    CHcombinedFLAG = true;
    [in_ca,inw_ca]=op_combineRcvrs(in,inw);
else
    CHcombinedFLAG = false;
    in_ca = in;
    inw_ca = inw;
end

% Temporarily combine averages - JND 9/2/2024
if inw.dims.averages~=0
    warning('WARNING: Averages will be temporarily combined to perform ecc.')
    AVcombinedFLAG = true;
    in_ca = op_averaging(in_ca);
    inw_ca = op_averaging(inw_ca);
else
    AVcombinedFLAG = false;
end

%save the phase as a vector of hard numbers.
inph=double(phase(inw_ca.fids));
figure; plot(inw_ca.t,inph); title('Phase of Water Suppressed Data')
axis tight; ylim([-180 180]); xlabel('Time (sec)'); ylabel('Phase (degrees)');

%choose the part of the phase function that is most linear
disp('Identify a section of linear values')
tmin=input('input min t value: ');
tmax=input('input max t value: ');

%now fit a straight line to the linear part of the phase function
p=polyfit(inw_ca.t(inw_ca.t>tmin & inw_ca.t<tmax), inph(inw_ca.t>tmin & inw_ca.t<tmax)',1);

%now fit a spline to approximate a smooth version of the phase function
pp=splinefit(inw_ca.t,inph,150);

%now subtract the line from the spline to get the eddy current related
%phase offset:
ecphase=ppval(pp,inw_ca.t)'-polyval(p,inw_ca.t)';
sz=size(in_ca.fids);
ecphase_rep=repmat(ecphase,[1 sz(2:end)]);
figure; plot(inw_ca.t,ecphase); title('Eddy Current Related Phase Offset')
axis tight; ylim([-180 180]); xlabel('Time (sec)'); ylabel('Phase (degrees)');

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

% Temporarily combine channels - JND 4/1/2025
if CHcombinedFLAG
    [out_ca,outw_ca]=op_combineRcvrs(out,outw);
else
    out_ca = out;
    outw_ca = outw;
end

% Temporarily combine averages - JND 9/2/2024
if AVcombinedFLAG
    out_ca = op_averaging(out_ca);
    outw_ca = op_averaging(outw_ca);
end

% Plot data
figure; plot(out_ca.t,phase(out_ca.fids)); title('Phase of Corrected Water Suppressed Data') 
axis tight; ylim([-180 180]); xlabel('Time (sec)'); ylabel('Phase (degrees)');
figure; plot(outw_ca.t,phase(outw_ca.fids));title('Phase of Corrected Water Unsuppressed Data') 
axis tight; ylim([-180 180]); xlabel('Time (sec)'); ylabel('Phase (degrees)');
plot(outw_ca.t,phase(outw_ca.fids));title('Phase of Water unsuppressed data') 
