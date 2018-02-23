function a = highlight(dataAxes, xSpan, ySpan, hColor)
% function highlight(dataAxes, xSpan, ySpan, hColor)
%
% changelog:
%   8/15/11 - originally written JDW
%
% highlights a section of datAxes with a transparent block.
% dataAxes - axes to be highlighted
% xSpan (optional) - the span on the x axis to highlight (defaults to XLim
% if empty)
% ySpan (optional) - the span on the y axis to highlight (defaults to YLim
% if empty)
% hColor (optional) - either a string or rgb color value (defaults to red
% if empty)
% DJC - 10-6-2017 - added return of figure handle for highlight

    fig = get(dataAxes, 'Parent');
    
    set(fig, 'currentaxes', dataAxes);
    
    daXLim = get(dataAxes, 'XLim');
    daYLim = get(dataAxes, 'YLim');

    if (nargin < 2 || isempty(xSpan))
        xSpan = daXLim;
    end
    if (nargin < 3 || isempty(ySpan))
        ySpan = daYLim;
    end
    if (nargin < 4 || isempty(hColor))
        hColor = 'r';
    end
    
    a = patch([xSpan(1); xSpan(1); xSpan(2); xSpan(2)], ...
              [ySpan(1); ySpan(2); ySpan(2); ySpan(1)], ...
              hColor);
          
    set(a, 'EdgeColor', hColor);
    set(a, 'FaceColor', hColor);
    
    set(a, 'EdgeAlpha', 'flat');
    set(a, 'FaceAlpha', 'flat');
    set(a, 'FaceVertexAlphaData', 0.5);

end