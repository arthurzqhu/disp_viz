function hAx = createOverlayAxis(ax_ll,ax_ur)

% create an overlay axis so that you can put a common x and/or y label for
% all subplots. do this after everything and right before you write
% xlabe = ... ylabel = ...
% 
% createOverlayAxis(ax_ll,ax_ur), where ax_ll = lower left axis, ax_ur =
% upper right axis. use gca in a for loop to get the axis handle.
% e.g.,
% for i = 1:4
%     subplot(2,2,i)
%     ax_handle(i) = gca;
% %     do whatever you want
% end
% hAx = createOverlayAxis(ax_handle(3),ax_handle(2));

hAx_ll = ax_ll.Position;
hAx_ur = ax_ur.Position;

height = hAx_ur(2)+hAx_ur(4)-hAx_ll(2);
width = hAx_ur(1)+hAx_ur(3)-hAx_ll(1);

hAx=axes('position',[hAx_ll(1) hAx_ll(2) width height],'visible','off'); 
hAx.XLabel.Visible='on';
hAx.YLabel.Visible='on';
axes(hAx)

end