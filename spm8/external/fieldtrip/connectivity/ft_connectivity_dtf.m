function [dtf, dtfvar, n] = ft_connectivity_dtf(input, varargin)

% FIXME build in proper documentation

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_connectivity_dtf.m 7123 2012-12-06 21:21:38Z roboos $

hasjack = ft_getopt(varargin, 'hasjack', 0);
powindx = ft_getopt(varargin, 'powindx');
% FIXME build in feedback
% FIXME build in proper documentation
% FIXME build in dDTF etc

siz    = [size(input) 1];
n      = siz(1);
ncmb   = siz(2);
outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));

if isempty(powindx)
  % data are represented as chan_chan_therest
  for j = 1:n
    tmph   = reshape(input(j,:,:,:,:), siz(2:end));
    den    = sum(abs(tmph).^2,2);
    tmpdtf = abs(tmph)./sqrt(repmat(den, [1 siz(2) 1 1 1]));
    outsum = outsum + tmpdtf;
    outssq = outssq + tmpdtf.^2;
    %tmp    = outsum; tmp(2,1,:,:) = outsum(1,2,:,:); tmp(1,2,:,:) = outsum(2,1,:,:); outsum = tmp;
    %tmp    = outssq; tmp(2,1,:,:) = outssq(1,2,:,:); tmp(1,2,:,:) = outssq(2,1,:,:); outssq = tmp;
    % FIXME swap the order of the cross-terms to achieve the convention such that
    % labelcmb {'a' 'b'} represents: a->b
  end
else
  error('linearly indexed data for dtf computation is at the moment not supported');
%FIXME this needs to be thought through -> also, as a multivariate measure
%a pairwise decomposition does not make sense, should this be dealt with by
%ft_connectivityanalysis?
  %   % data are linearly indexed
%   sortindx = [0 0 0 0];
%   for k = 1:ncmb
%     iauto1  = find(sum(cfg.powindx==cfg.powindx(k,1),2)==2);
%     iauto2  = find(sum(cfg.powindx==cfg.powindx(k,2),2)==2);
%     icross1 = k;
%     icross2 = find(sum(cfg.powindx==cfg.powindx(ones(ncmb,1)*k,[2 1]),2)==2);
%     indx    = [iauto1 icross2 icross1 iauto2];
%     
%     if isempty(intersect(sortindx, sort(indx), 'rows')),
%       sortindx = [sortindx; sort(indx)];
%       for j = 1:n
%         tmph    = reshape(input(j,indx,:,:), [2 2 siz(3:end)]);
%         den     = sum(abs(tmph).^2,2);
%         tmpdtf  = reshape(abs(tmph)./sqrt(repmat(den, [1 2 1 1])), [4 siz(3:end)]);
%         outsum(indx,:) = outsum(indx,:) + tmpdtf([1 3 2 4],:);
%         outssq(indx,:) = outssq(indx,:) + tmpdtf([1 3 2 4],:).^2;
%         % FIXME swap the order of the cross-terms to achieve the convention such that
%         % labelcmb {'a' 'b'} represents: a->b
%       end
%     end
%   end
end
dtf = outsum./n;

if n>1, %FIXME this is strictly only true for jackknife, otherwise other bias is needed
  if hasjack,
    bias = (n - 1).^2;
  else
    bias = 1;
  end
  dtfvar = bias.*(outssq - (outsum.^2)/n)./(n-1);
else
  dtfvar = [];
end

% this is to achieve the same convention for all connectivity metrics:
% row -> column
for k = 1:prod(siz(4:end))
  dtf(:,:,k)    = transpose(dtf(:,:,k));
  if ~isempty(dtfvar)
    dtfvar(:,:,k) = transpose(dtfvar(:,:,k));
  end
end
