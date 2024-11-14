%% Inputs
figs_en             = 0;

%% Resource Grid parameters
scs                 = 15e3;
bw                  = 40e6;
n_symbols           = 14;
preambleSymbolIdx   = 1;
dmrsSymbolIdx       = [3 6 9];
dmrsFreqGap         = 2;

%% Modulation Parameters
mod_order       = 16;

%% BB parameters
upSamplingFactor    = 2;

%% OFDM Parameters
cpSize              = 1/32;

%% FEC parameters
FEC.cr                      = 1/2;
FEC.viterbiTailBits         = 16;
FEC.scramblerInitialState   = 'CEA63';

%% RF parameters
Fc                          = 996e6;

%% Channel parameters
disableAllErrors          = 0; % if '1' no channel impairments

ChanParams.timingOffsetEn = true;
ChanParams.awgnNoiseEn    = true;
ChanParams.PhaseErrEn     = false;
ChanParams.CFOEn          = true;

%% Derived Parameters (Dont Edit anything from here onwards casually)
n_RBs               = n_RBCalc(scs,bw);
nFFT                = n_FFTCalc(n_RBs);
bb_Fs               = nFFT*scs;

bitsperConstPoint   = log2(mod_order);

if Fc < 2*bb_Fs*upSamplingFactor
    error('RF Sampling Rate is insufficient, Atleast required is %0.2f MHz\n',(2*bb_Fs*upSamplingFactor)/1e6);
end

if disableAllErrors
    ChanParams.timingOffsetEn = false; %#ok<UNRCH>
    ChanParams.awgnNoiseEn    = false;
    ChanParams.PhaseErrEn     = false;
    ChanParams.CFOEn          = false;
end

%% Bit processing
n_dataPos               = n_RBs*12*n_symbols - (numel(preambleSymbolIdx)*n_RBs*12 + numel(dmrsSymbolIdx)*n_RBs*12/dmrsFreqGap);
n_dataBits              = n_dataPos*bitsperConstPoint;
n_origBits              = n_dataBits*FEC.cr-FEC.viterbiTailBits;
origBits                = round(rand(1,n_origBits));
FEC.interleaveLength    = max(factor(n_dataBits));
bitProcessedData        = bitprocessingTx(origBits,FEC);

%% Modulation
bitProcessedData_1  = reshape(bitProcessedData,numel(bitProcessedData)/bitsperConstPoint,bitsperConstPoint);
X                   = bi2de(bitProcessedData_1)';
modData             = qammod(X,mod_order)/bitsperConstPoint;

preambleLength      = n_RBs*12;
n_preambleBits      = preambleLength*2;
preambleBits        = round(rand(1,n_preambleBits));
preambleBits_1      = reshape(preambleBits,numel(preambleBits)/2,2);
X                   = bi2de(preambleBits_1)';
preambleData        = qammod(X,4)/2; % Preamble is QPSK modulated

dmrsLength          = numel(dmrsSymbolIdx)*n_RBs*12/dmrsFreqGap;
n_dmrsBits          = dmrsLength*2;
dmrsBits            = round(rand(1,n_dmrsBits));
dmrsBits_1          = reshape(dmrsBits,numel(dmrsBits)/2,2);
X                   = bi2de(dmrsBits_1)';
dmrsData            = qammod(X,4)/2; % Preamble is QPSK modulated

%% Grid construction
gridBoundaries      = [1:n_RBs*12:n_RBs*12*n_symbols;n_RBs*12:n_RBs*12:n_RBs*12*n_symbols];
gridIndices         = 1:n_RBs*12*n_symbols;

preambleIndices     = gridBoundaries(1,preambleSymbolIdx):gridBoundaries(2,preambleSymbolIdx);

dmrsIndices         = [];
for ii = dmrsSymbolIdx
    dmrsIndices     = [dmrsIndices gridBoundaries(1,ii):dmrsFreqGap:gridBoundaries(2,ii)]; %#ok<AGROW>
end

dataIndices         = setdiff(gridIndices,[preambleIndices dmrsIndices]);

refGrid                     = zeros(1,n_RBs*12*n_symbols);
refGrid(preambleIndices)    = preambleData;
refGrid(dmrsIndices)        = dmrsData;
refGrid(dataIndices)        = modData;

refGrid                     = reshape(refGrid,n_RBs*12,n_symbols);

if figs_en
    grid2Plot = refGrid;
    grid2Plot(dataIndices) = 10+10i;
    grid2Plot(preambleIndices) = 100+100i;
    grid2Plot(dmrsIndices) = 5000+5000i;
    image(real(grid2Plot));
    title('baseband freq grid before OFDM mod');
end

%% OFDM modulation
% Zero padding
nRowsInGrid                 = size(refGrid,1);
nColsInGrid                 = size(refGrid,2);

nZerosUpperSide             = (nFFT - nRowsInGrid)/2;
nZerosDownSide              = (nFFT - nRowsInGrid)/2;

ZP_Frame                    = [zeros(nZerosUpperSide,nColsInGrid);refGrid;zeros(nZerosDownSide,nColsInGrid)];
% IFFT Shift
ShiftedFrame      = ifftshift(ZP_Frame,1);

% IFFT Operation
tdFrame           = ifft(ShiftedFrame,[],1);

% CP addition
cpVector          = repmat(nFFT*cpSize,1,nColsInGrid);
tdFramewithCP     = zeros((1+cpSize)*nFFT,nColsInGrid);
for ii = 1:nColsInGrid
    symFrame            = tdFrame(:,ii);
    symFrameWithCP      = [symFrame(end-cpVector(ii)+1:end);symFrame];
    tdFramewithCP(:,ii) = symFrameWithCP;
end

% final baseband time signal
finalTDFrame        = tdFramewithCP(:).';

%% Upsampling
tdFrameOut          = resample(finalTDFrame,upSamplingFactor,1); % upconversion
if figs_en
    figure()
    subplot 221
    plot(real(finalTDFrame));
    title('Time domain Baseband output ');
    subplot 222
    spectrumVisualizer(finalTDFrame,bb_Fs/mod_order)

    subplot 223
    plot(real(tdFrameOut));
    title('Output after Upsampling');
    subplot 224
    spectrumVisualizer(tdFrameOut,bb_Fs*upSamplingFactor)
end

preambleIndices     = 1:(1+cpSize)*nFFT*upSamplingFactor;
refTDPreamble       = tdFrameOut(preambleIndices);

%% RF conversion
tdFrameOutRF        = tdFrameOut.*exp(2i*pi*Fc/(bb_Fs*upSamplingFactor)*(0:length(tdFrameOut)-1));

%% Channel impairments
[chan_out, ChannelParams]    = channelImpairments(tdFrameOutRF,bb_Fs*upSamplingFactor,ChanParams);

%% RF Downconversion
rxIn                = chan_out.*exp(2i*pi*(-Fc)/(bb_Fs*upSamplingFactor)*(0:length(chan_out)-1));

%% Timing offset Estimation
corrOut             = xcorr(refTDPreamble,rxIn);
c1                  = corrOut.*conj(corrOut);
if figs_en
    figure()
    plot(c1,'-o');
    title('Correlation plot for timing sync')
end
[~,p]               = max(c1);
StartIdxOfPreamble  = length(rxIn)-p+1;
timeSyncedFrame     = rxIn(StartIdxOfPreamble:end);

%% Down Sampling
timeSyncedFrameBB   = resample(timeSyncedFrame,1,upSamplingFactor);

%% Grid processing
frameInGrid         = reshape(timeSyncedFrameBB,(1+cpSize)*nFFT,n_symbols);

%% OFDM Demod
% CP removal
cpVector          = repmat(nFFT*cpSize,1,n_symbols);
tdFramewithoutCP  = zeros(nFFT,n_symbols);
for ii = 1:n_symbols
    symFrame2               = frameInGrid(:,ii);
    symFrameWithoutCP       = symFrame2(cpVector(ii)+1:end);
    tdFramewithoutCP(:,ii)  = symFrameWithoutCP;
end

% FFT
fdFrame           = fft(tdFramewithoutCP,[],1);

% fft shifting
ShiftedFrame      = ifftshift(fdFrame,1);

% Zero padding removal
gridRx            = ShiftedFrame(nZerosUpperSide+1:end-nZerosDownSide,:);

%% Extracting Data and Ref Signals from Grid
rxData            = gridRx(dataIndices);
rxDmrsData        = gridRx(dmrsIndices);

rxDmrsData        = reshape(rxDmrsData,[],numel(dmrsSymbolIdx));
refDmrsData       = reshape(dmrsData,[],numel(dmrsSymbolIdx));

[freqErrVec, PhaseErrVec]  = freqErrEstimation(rxDmrsData,refDmrsData,dmrsSymbolIdx,scs);
freqErr                    = mean(freqErrVec);
phaseErr                   = mean(PhaseErrVec);


%% Demodulation
demodX            = qamdemod(rxData,mod_order);
demodBits         = de2bi(demodX);
demodBits         = demodBits(:)';
%% Bit processing
decBits           = bitprocessingRx(demodBits,FEC);

%% Performance Evaluation
bitsErr           = sum(origBits~=decBits);
fprintf('BER is %f (%d out of %d Bits)\n',bitsErr/numel(origBits),bitsErr,numel(origBits));
fprintf('\n##Channel Introduced Errors##\nTiming Delay = %d Samples\nSNR          = %.2f\nCFO          = %.2f\nPhase Err    = %.2f\n',ChannelParams.timesamplesOff,ChannelParams.snrVal,ChannelParams.CFO,ChannelParams.phaseErr);

fprintf('\nEstimation##\nTime Offset Estimated = %d Samples\nFreq Estimated        = %0.2f\n',StartIdxOfPreamble-1,freqErr);

