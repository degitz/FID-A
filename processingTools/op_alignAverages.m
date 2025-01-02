% op_alignAverages.m
% Jamie Near, McGill University 2014.
% Edits from
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% [out,fs,phs]=op_alignAverages(in,tmax,med,ref);
% 
% DESCRIPTION:
% Perform spectral registration in the time domain to correct frequency and
% phase drifts.  As described in Near et al.  Frequency and phase drift 
% correction of magnetic resonance spectroscopy data by spectral 
% registration in the time domain. Magn Reson Med 2015; 73(1):44-50.
%
% June 15th 2017:  Made the tmax and med arguments optional.  If tmax is 
% not specified, the value is determined automatically by finding the time
% at which the SNR of the FID drops permanently below 5.  This idea
% was suggested by Mark Mikkelsen.  Thanks Mark!!
% 
% INPUTS:
% in        = Input data structure.
% tmax      = Maximum time (s) in time domain to use for alignment.
%             (Optional. Default is the time at which SNR drops below 5)
% med       = Align averages to the median of the averages? ('y','n', 'a', or 
%             'r'). (Optional.  Default = 'n').  
%               - If you select 'n', all averages will be aligned to a 
%                 single average. The average chosen as the reference 
%                 average will be the one with the lowest 'unlikeness' 
%                 metric (see 'op_rmbadaverages.m'). 
%               - If select 'y', all averages will be aligned to the median
%                 of the averages.
%               - If you select 'a', all averages will be aligned to the 
%                 average of the averages.
%               - If you select 'r', all averages will be aligned to an
%                 externally provided reference spectrum.
% ref       = An externally provided reference spectrum that you would like
%             to align everything to (Required only if med = 'r').  
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.


function [out,fs,phs,tmax_est]=op_alignAverages(in,varargin)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.averages==0
    %DO NOTHING
    disp('WARNING:  No averages found.  Returning input without modification!');
    out=in;
    fs=0;
    phs=0;
    return
end

%%% Parse inputs - JND 9/2/2024
[tmax, tmax_est, med, ref] = parseInputs(in, varargin{:});

%%% Initialize data
dwelltime = in.dwelltime;
parsFit=[0,0];

if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

fs=zeros(in.sz(in.dims.averages),B);
phs=zeros(in.sz(in.dims.averages),B);
fids=zeros(in.sz(in.dims.t),1,B);

%%% Iterate through subspecs
for m=1:B
    % Set base
    switch med
        case {'y','Y'}
            disp('Aligning all averages to the median of the averages.');
            base=op_median(in);
        case {'a','A'}
            disp('Aligning all averages to the average of the averages.');
            base=op_averaging(in);
        case {'n','N'}
            %First find the average that is most similar to the total average:
            inavg=op_median(in);
            for k=1:in.sz(in.dims.averages)
                for l=1:B
                    metric(k,l)=sum((real(in.fids(in.t>=0 & in.t<=tmax,k,l))-(real(inavg.fids(inavg.t>=0 & inavg.t<=tmax,l)))).^2);
                end
            end
            [temp,ind_min]=min(metric(:,m));

            %Now set the base function using the index of the most similar average:
            disp(['Aligning all averages to average number ' num2str(ind_min) '.']);
            base=in;
            fids(:,ind_min,m)=in.fids(:,ind_min,m);
        case {'r','R'}
            disp('Aligning all averages to an externally provided reference spectrum.');
            base=ref;
    end

    % Isolate tmax points
    switch med
        case {'y','Y','a','A','r','R'}
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        case {'n','N'}
            base=[real(base.fids(in.t>=0 & in.t<tmax,ind_min,m));imag(base.fids(in.t>=0 & in.t<tmax,ind_min,m))];
    end

    % Change max iteration warning to error temporarily - JND 12/6/24
    warning('error','stats:nlinfit:IterationLimitExceeded');
    opts = statset('nlinfit'); opts.MaxIter = 400; maxitFLAG = false;

    %%% Fit data
    for n=1:in.sz(in.dims.averages)
        if n~=ind_min
            parsGuess=parsFit;
            %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
            start=in.fids(in.t>=0 & in.t<tmax,n,m);

            % Increase max iterations if it is reached - JND 12/6/24
            if ~maxitFLAG
                try
                    parsFit=nlinfit(start,base,@op_freqPhaseShiftComplexNest,parsGuess);
                catch ME
                    switch ME.identifier
                        case 'stats:nlinfit:IterationLimitExceeded'
                            maxitFLAG = true;
                            warning('on','stats:nlinfit:IterationLimitExceeded'); % Restore the warnings back to their previous (non-error) state - JND 12/6/24
                            parsFit=nlinfit(start,base,@op_freqPhaseShiftComplexNest,parsGuess,opts);
                    end
                end
            else
                parsFit=nlinfit(start,base,@op_freqPhaseShiftComplexNest,parsGuess,opts);
            end

            fids(:,n,m)=op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
            fs(n,m)=parsFit(1);
            phs(n,m)=parsFit(2);
            %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
        end
    end
end

warning('on','stats:nlinfit:IterationLimitExceeded'); % Restore the warnings back to their previous (non-error) state - JND 12/6/24

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqcorrected=1;


    function y=op_freqPhaseShiftComplexNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        shifted=addphase(fid.*exp(1i*t'*f*2*pi),p);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        y=addphase(fid.*exp(1i*t'*f*2*pi),p);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function [tmax, tmax_est, med, ref] = parseInputs(in, varargin)
        CHECK_tmx = true;
        CHECK_med = true;
        CHECK_ref = true;
        for i = 1:nargin-1
            % Check for tmax input
            if strcmpi(varargin{i}, 'tmax') && CHECK_tmx
                % Remove input from options
                CHECK_tmx = false;

                % Save variable
                tmax = varargin{i+1};
                tmax_est = NaN;

            % Check for ref input
            elseif strcmpi(varargin{i}, 'ref') && CHECK_ref
                % Remove input from options
                CHECK_ref = false;

                % Save variable
                ref = varargin{i+1};
            
            % Check for med input
            elseif strcmpi(varargin{i}, 'med') && CHECK_med
                % Remove input from options
                CHECK_med = false;

                % Save variable
                med = varargin{i+1};
            end
        end

        % Add default options based on which inputs weren't grabbed
        if CHECK_ref, ref = struct(); end
        if CHECK_med, med = 'n'; end
        if CHECK_tmx % Find the time at which the SNR drops below 5
            disp('tmax not supplied.  Calculating tmax....');

            % Calculate SNR
            sig = abs(in.fids);
            noise = std(real(in.fids(ceil(0.75*end):end,:,:)),[]);
            noise = mean(mean(mean(noise,2),3),4);
            snr = sig/noise;

            % Find tmax
            tmax_est = zeros(in.sz(in.dims.averages),1);
            for nn = 1:in.sz(in.dims.averages) % Changed for readability - JND 9/2/2024
                if any(snr(:,nn)>5)
                    N = find(snr(:,nn)>5);
                elseif any(snr(:,nn)>4)
                    N = find(snr(:,nn)>4);
                else % End loop if SNR is less than four - ETV 2024
                    % error('SNR is too low to perform frequency drift correction. Try again with a different file or without performing the correction')
                    [~, N] = max(snr(:,nn));
                end
                tmax_est(nn) = in.t(N(end));
            end
            tmax = median(tmax_est);
            disp(['tmax = ' num2str(tmax*1000) 'ms']);
        end
        if strcmpi(med,'r') % Combined two conditions into case-insensitive condition - JND 8/27/2024
            if CHECK_ref
                error('ERROR:  If using the ''r'' option for input variable ''med'', then a 4th input argument must be provided');
            end
        end
    end
end
