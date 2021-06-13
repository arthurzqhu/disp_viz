function bubblepie(xlist,ylist,slist,graph_data)
% A Bubble Plot, Pie Chart Combination
% bubblepie(xlist,ylist,slist,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
% 
% Creates a plot with pie charts at (xlist, ylist) using graph_data, having
% size of the pie scaled by slist.  Graph_labels contains a title for each
% pie chart (it can be empty [] to avoiid labeling), Graph_legend indicates
% the contents of graph_data, and lab is a binary value indicating whether
% pie chart labels are displayed. Pie colormaps can be modified afterwards
% with the colormap function or by acessing indivdual axes children.
% 
% Example:
% x = -pi:1:pi;
% x = x';
% y = sin(x);
% s = 1.1+cos(x);
% graph_labels = mat2cell(x,ones(1,length(x)),1);
% graph_data = 10*rand(length(x),3);
% graph_legend = {'one','two','three'};
% xlab = 'radians';
% ylab = 'sin(x)';
% lab = 1;
% 
% bubblepie(x,y,s,graph_data,graph_labels,graph_legend,xlab,ylab,lab)
% 
% title('BubblePie Plot')
% 
% To avoid titling pie charts pass in empty vector:
% bubblepie(x,y,s,graph_data,[],graph_legend,xlab,ylab,lab)
% 
% To indicate no pie labels are wanted by setting the parameter lab = 0:
% bubblepie(x,y,s,graph_data,graph_labels,graph_legend,xlab,ylab,0)
% 
%   Abraham Anderson
%   July 30, 2007, updated Oct 27, 2017

graph_max_size = 0.02;
graph_min_size = 0.01;
graph_range = graph_max_size-graph_min_size;

canvas_max = .9;
canvas_min = 0.1;
canvas_range = canvas_max-canvas_min;

maxx = max(xlist);
maxy = max(ylist);
minx = min(xlist);
miny = min(ylist);
maxs = max(slist);

if maxy==miny
    maxy = 1;
    miny = 0;
end

if maxx==minx
    maxx = 1;
    minx = 0;
end

% figure
% h0 = axes('position',[canvas_min,canvas_min,canvas_range,canvas_range], ...
%     'xlim',[minx maxx]);

% set(get(gca,'XLabel'),'String',xlab)
% set(get(gca,'YLabel'),'String',ylab)
ax_pos = get(gca,'Position');
ax_xspan = ax_pos(3);
ax_yspan = ax_pos(4);
% text(0.1*maxx+minx,maxy-0.1*maxy,{'PieChart Groups (CCW):' graph_legend{:}},'verticalalignment','top');
for i = 1:size(graph_data,1)
        s = slist(i)/maxs*graph_range+graph_min_size;
        x = (xlist(i)-minx)/(maxx-minx)*ax_xspan+ax_pos(1)-s/2;
        y = (ylist(i)-miny)/(maxy-miny)*ax_yspan+ax_pos(2)-s/2;
        d = graph_data(i,:);
        if sum(d)==0
            continue
        end
%         d(d==0) = 0.1; im sorry why would you do this?!
        axes('position',[x y s s])
        pie(d,repmat({''},length(d),1))
end

set(gcf,'Currentaxes',gca)
