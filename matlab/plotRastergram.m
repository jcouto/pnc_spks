function [spikes,ln]=plotRastergram(raster,trial_index,offset,color)
% Plot a rastergram from spike trains
% Joao Couto 2015
if ~iscell(raster)
    trials = unique(trial_index);
    spikes = cell(1,length(trials));
    for ii = 1:length(trials)
        spikes{ii} = raster(trial_index==trials(ii));
    end
else
    spikes = raster;
end
if ~exist('color','var')
    color = 'k';
end

if ~exist('offset','var')
    offset = 0;
end
%
linesize = .2;
minimum=-.5;
maximum=.5;
ln = [];

for ii=1:length(spikes)
    for jj=1:length(spikes{ii})
        ln = [ln,line(spikes{ii}(jj).*[1,1],...
            [minimum,maximum]+offset,...
            'color',color,'linewidth',linesize)];
    end
    offset=offset+1;
end

%xlabel('Time (s)','fontname','Arial','fontsize',12,'fontweight','bold')
%ylabel('Trials','fontname','Arial','fontsize'
