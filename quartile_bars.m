function quartile_bars(dat, axis, fig, transparency, marker_size, linewidth, inst_pos)
% input: an array of data; which axis ('x', or 'y') you want to put the quartile bars;
% 'fig' as the figure handle is used to determine which figure to put the quartile bars.
% defaults: marker_size = 20, linewidth = 3, inst_pos (position of instruction) is upper 
% right. need to specify the relative position of the instruction box on the figure; e.g.,
% default is [0.725, 1, 0.6, 1].
%
% output: quartile bars on the existing axes of the plot

figure(fig)
    

dat_quartl = num2cell(prctile(dat,0:25:100));
[dat_min, dat_25, dat_med, dat_75, dat_max] = deal(dat_quartl{:});

xax_lim = get(gca,'xlim');
yax_lim = get(gca,'ylim');

if nargin < 7
    inst_pos = [0.725, 1, 0.6, 1];
    
    if nargin < 6
        linewidth = 3;
        
        if nargin < 5
            marker_size = 20;
            if nargin < 4
                transparency = 0.5; %havent implement this into the code below
            end
        end
    end
end

bar_array_cell = {[dat_min, dat_25]; dat_med; [dat_75, dat_max]};

% determine where on the plot to put the quartiles
if axis == 'x'
    [x1, x2, x3] = deal(bar_array_cell{:});
    y1 = [yax_lim(1) yax_lim(1)]; 
    y2 = yax_lim(1); y3 = y1;
elseif axis == 'y'
    [y1, y2, y3] = deal(bar_array_cell{:});
    x1 = [xax_lim(1) xax_lim(1)]; 
    x2 = xax_lim(1); x3 = x1;
else
    error('only 2D plots are currently supported, the second input needs to be either "x" or "y"');
end

% plotting those quartiles
line(x1, y1, 'linewidth', linewidth)
plot(x2, y2, '.', 'MarkerSize', marker_size)
line(x3, y3, 'linewidth', linewidth)

% instruction
box_minx = (xax_lim(2)-xax_lim(1))*inst_pos(1)+xax_lim(1);
box_maxx = (xax_lim(2)-xax_lim(1))*inst_pos(2)+xax_lim(1);
box_miny = (yax_lim(2)-yax_lim(1))*inst_pos(3)+yax_lim(1);
box_maxy = (yax_lim(2)-yax_lim(1))*inst_pos(4)+yax_lim(1);
box_wid = box_maxx-box_minx;
box_len = box_maxy-box_miny;

line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.15 box_miny+box_len*.4],...
    'linewidth',linewidth)
plot(box_minx+box_wid*.2,box_miny+box_len*.5,'.','MarkerSize',marker_size)
line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.6 box_miny+box_len*.85],...
    'linewidth',linewidth)

text(box_minx+box_wid*.3,box_miny+box_len*.15, 'min')
text(box_minx+box_wid*.3,box_miny+box_len*.4, '1st quartile')
text(box_minx+box_wid*.3,box_miny+box_len*.5, 'median')
text(box_minx+box_wid*.3,box_miny+box_len*.6, '3rd quartile')
text(box_minx+box_wid*.3,box_miny+box_len*.85, 'max')

rectangle('Position',[box_minx box_miny box_wid box_len],'LineWidth',linewidth/2)

end