function out = bitprocessingTx(inpt,FEC)

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

bitsScrambled   = scrambler(inpt,scramblerInitialState);
convOut         = convenc([bitsScrambled zeros(1,viterbiTailBits)],tr,puncPat);
out             = RowCol_Interleaver(convOut,interleaveLength);
end
