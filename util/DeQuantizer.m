%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantization interval: l
% number of bits: n

function Decoded = DeQuantizer(n, l, coded)

Decoded = double(coded).*l/(2^n-2);

end