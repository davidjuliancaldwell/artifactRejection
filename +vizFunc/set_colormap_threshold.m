function set_colormap_threshold(fi,tr,bounds,col)
% function set_colormap_threshold(fi,tr,bounds,col)
% fi = figure handle
% tr = [lower_bound upper_bound]
% bounds = [min max]
% col =[r g b];
cmin=bounds(1);
cmax=bounds(2);
dc=cmax-cmin;

caxis([cmin cmax]);
figure(fi);
ax=get(fi,'Children');
% cm=colormap(jet(256));
cm = colormap;

cm = interp1(1:size(cm, 1), cm, linspace(1, size(cm, 1), 256));

ix1=round((tr(1)-cmin)/dc*256);
ix2=round((tr(2)-cmin)/dc*256);
if(ix1<1)
    ix1=1;
end
if(ix2>256)
    ix2=256;
end
if(ix2>ix1)
    cm(ix1:ix2,:)=ones(length(ix1:ix2),1)*col;
end
for i=1:length(ax)
colormap(cm);
end
end

