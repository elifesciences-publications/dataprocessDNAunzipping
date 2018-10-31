function force = getWLCapproxF(zL,lpkT,K)
% get approximate force for a given extension of stretchable WLC
% using Marko-Siggia interpolation formula
% end extension: z/L
% lpkT = lp/kT
% K = stretch constant (in units of force)

%%
A = 4*lpkT/K^2;
B = 1/K^2*(1-4*zL) + 8*lpkT/K.*(1-zL);
C = 2/K*(1-zL).*(1-4*zL) + 4*lpkT*(1-zL).^2;

for cc = 1:length(zL)
    coeff = [A,B(cc),C(cc),-1];
    tmp = roots(coeff);
    
    % get real root
    [~,b] = min(abs(imag(tmp)));
    force(cc) = real(tmp(b));
end
