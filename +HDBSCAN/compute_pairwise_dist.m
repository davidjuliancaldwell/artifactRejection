function D = compute_pairwise_dist( X,varargin )
% D = compute_pairwise_dist( X,(Y) )
%
% computes a (fast!) pairwise euclidean distance matrix of the inputs.
% If more than one matrix is provided, will compute pairwise distances
% between matrices. Else, computes a symmetric pairwise distance of X
%
% Both X and the optional matrix Y must have points along the rows, and
% dimensions of points along columns.
%
% Dynamic time warping is an additional distance metric


if ~isempty(varargin{2})
    metric = varargin{2};
else
    metric = 'eucl';
end

if strcmp(metric,'eucl')
    
    if isempty(varargin{1})
        d = sum( X.^2,2 );
        D = real( sqrt( bsxfun( @plus,d,d' ) - 2*(X * X') ) );
    else
        Y = varargin{2};
        if size( Y,2 ) ~= size( X,2 )
            error( 'X and Y must have equal number of dimensions' );
        end
        
        D = real( sqrt( bsxfun( @plus, sum( X.^2,2 ),...
            bsxfun( @minus, sum( Y.^2,2 )', 2*(X * Y') ) ) ) );
    end
    
elseif strcmp(metric,'corr')
    
    if isempty(varargin{1})
        D = corr(X');
        D = sqrt(1-D);
                D = (D + D.')/2;

    else
        Y = varargin{2};
        if size( Y,2 ) ~= size( X,2 )
            error( 'X and Y must have equal number of dimensions' );
        end
        
        D =  corr(X',Y') ;
        D = sqrt(1-D);
        D = (D + D.')/2;


    end
    
elseif strcmp(metric,'dtw')
    
    if isempty(varargin{1})
        % dynamic time warping
        sX = size(X,1);
        D = zeros(sX,sX);
        for indexOut = 1:sX
            for indexIn = indexOut:sX
                D(indexOut,indexIn) = dtw.dtw_c(X(indexOut,:)',X(indexIn,:)');
            end
        end
    else
        Y = varargin{2};
        
        sX = size(X,1);
        D = zeros(sX,sX);
        for indexOut = 1:sX
            for indexIn = indexOut:sX
                D(indexOut,indexIn) = dtw.dtw_c(X(indexOut,:)',Y(indexIn,:)');
            end
        end
    end
    
    %
    D = triu(D)+triu(D,1)';
    
end

end