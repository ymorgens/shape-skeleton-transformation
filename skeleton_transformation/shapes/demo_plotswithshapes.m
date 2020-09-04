NumShapesOnAxis = 7;
x = 1:NumShapesOnAxis;
y = 1:NumShapesOnAxis;
h = plot(x,y,'ro');h.LineWidth=12; hold on;
plot(x,y,'k:')
ylabel('Dependent','FontSize',18);
xlabel('Independent','FontSize',18);
title('QUEENS OF THE STONE AGE!!!!');
hold on;
axis square;
ax = gca;
YTick = ax.YTick;
XTick = ax.XTick;
dy = YTick(2) - YTick(1);
dx =  XTick(2) - XTick(1);
ax.FontSize = 25;

load('novobjs','obj');

% add stimului below
for n = 1:NumShapesOnAxis
points = 0.25*adjustContourPos(obj{n}.shape);
points(:,1) = dx*points(:,1) + n;
points(:,2) = dy*points(:,2)-dy/1.5;% - d/2;
hold on;
patch(points(:,1), points(:,2),'k');

end