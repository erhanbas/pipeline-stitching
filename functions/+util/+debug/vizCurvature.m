
figure(304),
if iadj==1
    % range
    ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
    yl = [0 dims(setdiff([1 2],iadj))];
    xlim(xl), ylim(yl)
    
    subplot(2,3,[1 3])
    cla
    plot(y,x,'+')
    hold on
    plot(y(outliers),x(outliers),'ro')
    plot(y_range,x_range,'m-','LineWidth',4)
    plot(y_range_init,x_range,'g-','LineWidth',2)
    daspect([1 50 1])
    legend('matched feats','outliers','estimated model','initial model')
    
elseif iadj==2
    % range
    ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
    yl = [0 dims(setdiff([1 2],iadj))];
    xlim(yl), ylim(xl)
    
    subplot(2,3,[2 3])
    cla
    plot(x,y,'+')
    hold on
    plot(x(outliers),y(outliers),'ro')
    plot(x_range,y_range,'m-','LineWidth',4)
    plot(x_range,y_range_init,'g-','LineWidth',2)
    daspect([30 1 1])
    legend('matched feats','outliers','estimated model','initial model')
end