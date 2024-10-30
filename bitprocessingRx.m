function out = bitprocessingRx(inpt,FEC)

cr                      = FEC.cr;
viterbiTailBits         = FEC.viterbiTailBits;
scramblerInitialState   = FEC.scramblerInitialState;
interleaveLength        = FEC.interleaveLength;

tr              = poly2trellis(7,[171,133]);
if cr == 1/2
    puncPat = [1 1];
elseif cr == 2/3
    puncPat = [1 1 0 1];
elseif cr == 3/4
    puncPat = [1 1 0 1 1 0];
elseif cr == 5/6
    puncPat = [1 1 0 1 1 0 0 1 1 0];
end

deintData       = RowCol_DeInterleaver(inpt,interleaveLength);
vitOut          = vitdec(deintData,tr,64,'term','hard',puncPat);
out             = descrambler(vitOut(1:end-viterbiTailBits),scramblerInitialState);
end
