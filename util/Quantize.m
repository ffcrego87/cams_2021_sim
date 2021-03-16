%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantization interval: l
% number of bits: n

function coded = Quantize(n, l, x)

delta = l./(2^n-2);

coded = int64(sign(x).*min(floor(1/2+abs(x)./delta),2^(n-1)-1));

end