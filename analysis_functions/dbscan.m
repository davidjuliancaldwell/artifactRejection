% DBSCAN DBSCAN clustering algorithm
%
% Usage:  [C, ptsC, centres] = dbscan(P, E, minPts)
%
% Arguments:
%         P - dim x Npts array of points.
%         E - Distance threshold.
%    minPts - Minimum number of points required to form a cluster.
%
% Returns:
%         C - Cell array of length Nc listing indices of points associated with
%             each cluster.
%      ptsC - Array of length Npts listing the cluster number associated with
%             each point.  If a point is denoted as noise (not enough nearby
%             elements to form a cluster) its cluster number is 0.
%   centres - dim x Nc array of the average centre of each cluster.

% Reference:
% Martin Ester, Hans-Peter Kriegel, J�rg Sander, Xiaowei Xu (1996). "A
% density-based algorithm for discovering clusters in large spatial databases
% with noise".  Proceedings of the Second International Conference on Knowledge
% Discovery and Data Mining (KDD-96). AAAI Press. pp. 226-231.
% Also see: http://en.wikipedia.org/wiki/DBSCAN

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% PK January 2013

function [C, ptsC, centres] = dbscan(P, E, minPts,distanceMetric) % DJC - add distance

if ~exist('distanceMetric','var')
    distanceMetric = 'euclidean';
end

[dim, Npts] = size(P);

ptsC  = zeros(Npts,1);
C     = {};
Nc    = 0;               % Cluster counter.
Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.

for n = 1:Npts
    if ~Pvisit(n)                            % If this point not visited yet
        Pvisit(n) = 1;                       % mark as visited
        neighbourPts = regionQuery(P, n, E,distanceMetric); % and find its neighbours
        
        if length(neighbourPts) < minPts-1  % Not enough points to form a cluster
            ptsC(n) = 0;                    % Mark point n as noise.
            
        else                % Form a cluster...
            Nc = Nc + 1;    % Increment number of clusters and process
            % neighbourhood.
            
            C{Nc} = [n];    % Initialise cluster Nc with point n
            ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
            
            ind = 1;        % Initialise index into neighbourPts array.
            
            % For each point P' in neighbourPts ...
            while ind <= length(neighbourPts)
                
                nb = neighbourPts(ind);
                
                if ~Pvisit(nb)        % If this neighbour has not been visited
                    Pvisit(nb) = 1;   % mark it as visited.
                    
                    % Find the neighbours of this neighbour and if it has
                    % enough neighbours add them to the neighbourPts list
                    neighbourPtsP = regionQuery(P, nb, E,distanceMetric);
                    if length(neighbourPtsP) >= minPts
                        neighbourPts = [neighbourPts  neighbourPtsP];
                    end
                end
                
                % If this neighbour nb not yet a member of any cluster add it
                % to this cluster.
                if ~ptsC(nb)
                    C{Nc} = [C{Nc} nb];
                    ptsC(nb) = Nc;
                end
                
                ind = ind + 1;  % Increment neighbour point index and process
                % next neighbour
            end
        end
    end
end

% Find centres of each cluster
centres = zeros(dim,length(C));
for n = 1:length(C)
    for k = 1:length(C{n})
        centres(:,n) = centres(:,n) + P(:,C{n}(k));
    end
    centres(:,n) = centres(:,n)/length(C{n});
end

end % of dbscan

%------------------------------------------------------------------------
% Find indices of all points within distance E of point with index n
% This function could make use of a precomputed distance table to avoid
% repeated distance calculations, however this would require N^2 storage.
% Not a big problem either way if the number of points being clustered is
% small.   For large datasets this function will need to be optimised.

% Arguments:
%              P - the dim x Npts array of data points
%              n - Index of point of interest
%              E - Distance threshold

function neighbours = regionQuery(P, n, E,distanceMetric)

E2 = E^2;
[dim, Npts] = size(P);
neighbours = [];
vector_pts = 1:size(P,2);

switch distanceMetric
        
    case 'euclidean'
        eucl = zeros(1,size(P,2));
        eucl(n) = nan;
        v = P(:,vector_pts(vector_pts~=n)) - repmat(P(:,n),1,size(P(:,vector_pts(vector_pts~=n)),2));
        eucl(vector_pts(vector_pts~=n)) = sum(v.*v);
       [~,neighbours] = find(eucl<E2);

        
    case 'corr'
        cor = zeros(1,size(P,2));
        cor(vector_pts(vector_pts~=n)) = corr(P(:,n),P(:,vector_pts(vector_pts~=n)));
        [~,neighbours] = find(cor>E);
        
    case 'cosine'
        cosTheta = zeros(size(P,2),1);
        % cosTheta = dot(P(:,n),P(:,vector_pts(vector_pts~=n)),2)/(norm(P(:,n))*norm(P(:,vector_pts(vector_pts~=n))));
        denominator = sqrt(sum(P(:,n).*P(:,n))).*sqrt(sum(P(:,vector_pts(vector_pts~=n)).*P(:,vector_pts(vector_pts~=n))));
        numerator = (P(:,n)'*P(:,vector_pts(vector_pts~=n)));
        cosTheta(vector_pts(vector_pts~=n)) = numerator./denominator;
        
     %   cosTheta(vector_pts(vector_pts~=n)) = (P(:,n)'*P(:,vector_pts(vector_pts~=n)))'./((sqrt(P(:,n)'*...
      %      P(:,n))).*sqrt(diag(P(:,vector_pts(vector_pts~=n))'*P(:,vector_pts(vector_pts~=n)))));
        [~,neighbours] = find(cosTheta'>E);
        
end



end % of regionQuery
