function [cspks,edges] = convolve_spiketrain(spks, binsize, duration, kern)
% Convolves spike trains by a specified kernel.
% Default is a gaussian kernel, specified by.
% Spks should be in seconds
% Joao Couto 2015

if ~exist('kern','var')
    sigma = 0.02; % in seconds!
    tt = -(2*sigma)/binsize:binsize:(2*sigma)/binsize;
    kern = normpdf(tt,0,sigma);
    kern = kern/sum(kern);
end
edges = (0:binsize:duration);
if ~iscell(spks)
    spks = {spks};
end
cspks = zeros(length(edges),length(spks));


for trial = 1:length(spks)
    % This mess is just to deal with empty spiketrains because of matlab.
    tmp = histc(spks{trial}, edges)./binsize;
    if ~isempty(tmp)
        cspks(:,trial) = conv(tmp,kern,'same');
    end
end
