function [coh_f, coh] = calCoherence(C,P1,P2,deltaF,...
                                     T,w1w2bar,w1w2squaredbar,mask)

% calcluates the coherence and integrated coherence:
%
%           |C|^2      w1w2bar^2
% coh(f) = ------- * --------------
%           P1 P2    w1w2squaredbar
%
% coh    = T int_{-infty}^{infty} coh(f) df
%
% Assumes that C,P1,P2 are all of the same length and units (strain^2/Hz).
% In particular C has to be already corrected by 1/w1w2bar.
% This was NOT the case for old stochastic.m and C=s1^*s2,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh_f = (C.*conj(C)./(P1.*P2)) .* mask *(w1w2bar^2 / w1w2squaredbar);
coh   = 2*T*sum(coh_f*deltaF); % factor of 2 to include +/- freqs

end
