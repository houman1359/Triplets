
function plot_grid(GRIDs,h,m,n,p)
figure(h);
subplot(m,n,p)
% gridlines ---------------------------
hold on
g_y=GRIDs(1,:); % user defined grid Y [start:spaces:end]
g_x=GRIDs(2,:); % user defined grid X [start:spaces:end]
for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'k:') %y grid lines
   hold on    
end
for i=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'k:') %x grid lines
   hold on    
end
% print(1,'-dpng','-r300','K1') %save plot as png (looks better)
title('grid')