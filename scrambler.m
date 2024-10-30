function [data_out] = scrambler(data_in,scramblerInitialState)

scrambler_state             = de2bi(hex2dec(scramblerInitialState),23,'left-msb');
data_out                    = zeros(size(data_in));
for ii = 1:length(data_in)
    data_out(ii)            = xor(data_in(ii),xor(scrambler_state(18),scrambler_state(23)));
    scrambler_state(2:end)  = scrambler_state(1:end-1);
    scrambler_state(1)      = data_out(ii);
end

end
