function [out,params] = channelImpairments(tdFrameOut,rf_outFreq,ChanParams)

out                     = tdFrameOut;
params.timesamplesOff   = 0;
params.snrVal           = 0;
params.CFO              = 0;
params.phaseErr         = 0;

if ChanParams.timingOffsetEn
    timesamplesOff          = randi([1000,3000],1,1);
    out                     = [complex(rand(1,timesamplesOff)/50,rand(1,timesamplesOff)/50) out];
    params.timesamplesOff   = timesamplesOff;
end

if ChanParams.awgnNoiseEn
    snrVal          = 5 + (10 - 5) * rand(1,1); % random snr value between 5dB and 40dB
    out             = awgn(out,snrVal,'measured');
    params.snrVal   = snrVal;
end

if  ChanParams.CFOEn
    df                  = 10 + (40 - 10) * rand(1,1); % random frequency between 50 and 400Hz
    out                 = out.*exp(2i*pi*df/(rf_outFreq)*(0:length(out)-1));
    params.CFO          = df;
end

if  ChanParams.PhaseErrEn
    theta               = pi/32 + (2*pi - pi/32) * rand(1,1); % random Phase between pi/32 and 2pi
    out                 = out*exp(2i*theta);
    params.phaseErr     = theta;
end

end
