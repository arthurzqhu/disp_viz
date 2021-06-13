function cbar_pos = getCbarPos(ax_above,ax_below)

ax_a_pos = get(ax_above,'position');
y_top = ax_a_pos(2) + ax_a_pos(4);

ax_b_pos = get(ax_below,'position');
y_bot = ax_b_pos(2);

x_left = ax_a_pos(1)+ax_a_pos(3)*1.05;

cbar_pos = [x_left, y_bot, 0.0143, y_top-y_bot];

end