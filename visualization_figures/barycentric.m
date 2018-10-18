function barycentric()
    clf;
    
    subplot(1,3,1);
    xp=[1];  yp=[1];  dx=[0];  dy=[0];
    draw_bary2D(0,0,'b-',xp,yp,dx,dy);
    draw_bary2D(0,1,'y-',xp,yp,dx,dy);
    draw_bary2D(1,0,'g-',xp,yp,dx,dy);
    draw_bary2D(1,1,'w-',xp,yp,dx,dy);
    
%%
    subplot(1,3,2);cla
    xp=[1,1,2,1];  yp=[1,2,1,0];  dx=[-0.25,0,-0.25,-0.25];  dy=[-0.25,-0.25,0.25,0];
%     xp=[1,1,2,1];  yp=[1,2,1,0];  
%     dx=randn(1,4); dy = randn(1,4)
    draw_bary2D(0,0,'b-',xp,yp,dx,dy);
    draw_bary2D(0,1,'y-',xp,yp,dx,dy);
    draw_bary2D(1,0,'g-',xp,yp,dx,dy);
    draw_bary2D(1,1,'w-',xp,yp,dx,dy);
%%
    
    subplot(1,3,3);
    draw_bary3D(0,0,0,'k');
    draw_bary3D(0,1,0,'r');
    draw_bary3D(1,0,0,'b');
    draw_bary3D(1,1,0,'g');
    draw_bary3D(0,0,1,'c');
    draw_bary3D(0,1,1,'m');
    draw_bary3D(1,0,1,'y');
    draw_bary3D(1,1,1,'w');
end

function [xs,ys] = move_coord(xs,ys,xp,yp,dx,dy)
    for is=1:length(xs)
        for ip=1:length(xp)
            if xs(is)==xp(ip) && ys(is)==yp(ip)
                xs(is)=xs(is)+dx(ip);
                ys(is)=ys(is)+dy(ip);
            end
        end
    end
end

function draw_bary2D(x0,y0,c,xp,yp,dx,dy)
    hold on;
    axis square;  axis off;

    % square
    [xm,ym] = move_coord([x0,x0+1],[y0,y0],xp,yp,dx,dy);
    plot(xm,ym,'k-','linewidth',1);
    [xm,ym] = move_coord([x0+1,x0+1],[y0,y0+1],xp,yp,dx,dy);
    plot(xm,ym,'k-','linewidth',1);
    [xm,ym] = move_coord([x0+1,x0],[y0+1,y0+1],xp,yp,dx,dy);
    plot(xm,ym,'k-','linewidth',1);
    [xm,ym] = move_coord([x0,x0],[y0+1,y0],xp,yp,dx,dy);
    plot(xm,ym,'k-','linewidth',1);

    % diagonal
    if mod(x0+y0,2)==1
        [xm,ym] = move_coord([x0+1,x0],[y0,y0+1],xp,yp,dx,dy);
        plot(xm,ym,c,'linewidth',1)
    else
        [xm,ym] = move_coord([x0,x0+1],[y0,y0+1],xp,yp,dx,dy);
        plot(xm,ym,c,'linewidth',1)
    end
end

function draw_bary3D(x0,y0,z0,c)
    hold on;
    axis square;  axis off;

    % vertical posts
    plot3([x0,x0],[y0,y0],[z0,z0+1],'k-');
    plot3([x0,x0],[y0+1,y0+1],[z0,z0+1],'k-');
    plot3([x0+1,x0+1],[y0,y0],[z0,z0+1],'k-');
    plot3([x0+1,x0+1],[y0+1,y0+1],[z0,z0+1],'k-');

    % bottom square
    plot3([x0,x0+1],[y0,y0],[z0,z0],'k-');
    plot3([x0+1,x0+1],[y0,y0+1],[z0,z0],'k-');
    plot3([x0+1,x0],[y0+1,y0+1],[z0,z0],'k-');
    plot3([x0,x0],[y0+1,y0],[z0,z0],'k-');

    % top square
    plot3([x0,x0+1],[y0,y0],[z0+1,z0+1],'k-');
    plot3([x0+1,x0+1],[y0,y0+1],[z0+1,z0+1],'k-');
    plot3([x0+1,x0],[y0+1,y0+1],[z0+1,z0+1],'k-');
    plot3([x0,x0],[y0+1,y0],[z0+1,z0+1],'k-');

    % inner tetrahedron
    if mod(x0+y0+z0,2)==1
        fill3([x0+1,x0,x0],[y0,y0+1,y0],[z0,z0,z0+1],c)
        fill3([x0+1,x0,x0+1],[y0,y0,y0+1],[z0,z0+1,z0+1],c)
        fill3([x0,x0,x0+1],[y0+1,y0,y0+1],[z0,z0+1,z0+1],c)
        fill3([x0+1,x0+1,x0],[y0,y0+1,y0+1],[z0,z0+1,z0],c)
    else
        fill3([x0,x0+1,x0],[y0,y0,y0+1],[z0,z0+1,z0+1],c)
        fill3([x0,x0+1,x0],[y0,y0+1,y0+1],[z0,z0,z0+1],c)
        fill3([x0,x0+1,x0+1],[y0,y0+1,y0],[z0,z0,z0+1],c)
        fill3([x0+1,x0+1,x0],[y0+1,y0,y0+1],[z0,z0+1,z0+1],c)
    end
end