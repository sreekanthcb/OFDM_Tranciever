function [freqErrVec, PhaseErrVec] = freqErrEstimation(rxData,refData,dmrsSymbolIdx,scs)

if size(rxData) ~= size(refData)
    error('Recived Data and Reference data are not of same size')
end
nSymbols = size(rxData,2);

freqErrVec = [];
PhaseErrVec = [];
for ii = 1:nSymbols-1
    c_out1 = sum(refData(:,ii).*conj(rxData(:,ii)));
    phase1 = angle(c_out1);
    PhaseErrVec = [PhaseErrVec c_out1]; %#ok<AGROW>

    c_out2 = sum(refData(:,ii+1).*conj(rxData(:,ii+1)));
    phase2 = angle(c_out2);
    PhaseErrVec = [PhaseErrVec c_out2]; %#ok<AGROW>

    nSymbolsBetween = dmrsSymbolIdx(ii+1)- dmrsSymbolIdx(ii);
    timeDuration    = nSymbolsBetween/scs;
    freqErrVec = [freqErrVec (phase2-phase1)/(2*pi*timeDuration)]; %#ok<AGROW>
end
PhaseErrVec = unique(PhaseErrVec);
freqErrVec  = -freqErrVec;
end

