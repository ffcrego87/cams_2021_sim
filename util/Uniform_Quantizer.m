%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform quantizer: Q
% quantization interval: l
% mid-value: x_bar
% number of bits: n

function [x_quantized,sat] = Uniform_Quantizer(n, l, x_bar, x)

delta = l/(2^n-2);

sat = max(floor(1/2+abs(x-x_bar)./delta)>=2^(n-1)-1);

x_quantized = x_bar + delta*sign(x-x_bar).*min(floor(1/2+abs(x-x_bar)./delta),2^(n-1)-1);
%disp(min(floor(1/2+abs(x-x_bar)./delta),2^(n-1)-1))
end