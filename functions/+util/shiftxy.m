function [xshift,yshift] = shiftxy(xy,centxy,beta,order,dims)
if nargin<4
    order = 1;
end

if centxy(2)<=eps % no curvature
    weightx = ((xy(:,1)-dims(1)/2).^order);
else
    weightx = ((xy(:,1)-centxy(2)).^order);
end
if centxy(1)<=eps % no curvature
    weighty = (sign((xy(:,2)-dims(2)/2)).*abs(xy(:,2)-dims(2)/2).^1);
else
    weighty = (sign((xy(:,2)-centxy(1))).*abs(xy(:,2)-centxy(1)).^1);
end
xshift = beta(1)*weightx.*((xy(:,2)-centxy(1)).^2);
yshift = beta(2)*weighty.*((xy(:,1)-centxy(2)).^2);
if nargin>4
    xshift = reshape(xshift,dims([2 1]));
    yshift = reshape(yshift,dims([2 1]));
end
