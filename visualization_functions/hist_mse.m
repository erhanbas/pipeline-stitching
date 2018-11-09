function hist_mse(C_mse,A_mse,S_mse,finterior)

clc
figure(123), clf

ax = gca;
h_ctrl = histogram(C_mse(finterior));
hold on
h_aff = histogram(A_mse(finterior));
h_st = histogram(S_mse(finterior));


C_st = [mean(C_mse(finterior),'omitnan') std(C_mse(finterior),'omitnan')];
A_st = [mean(A_mse(finterior),'omitnan') std(A_mse(finterior),'omitnan')];
S_st = [mean(S_mse(finterior),'omitnan') std(S_mse(finterior),'omitnan')];
leg1 = sprintf('Control points: %1.3f \\pm %1.3f',C_st(1),C_st(2));
leg2 = sprintf('Affine: %1.3f \\pm %1.3f',A_st(1),A_st(2));
leg3 = sprintf('Stage: %1.3f \\pm %1.3f',S_st(1),S_st(2));

leg = legend(leg1,leg2,leg3);
% legend('test \pm')
leg.FontSize = 32;
h_ctrl.BinWidth = .2;
h_aff.BinWidth = .2;
h_st.BinWidth = .2;
xlim([0 20])

title('Residual out of 11657 tiles', 'FontSize', 30)
xlabel('Residual magnitude in \mum', 'FontSize', 30)
ax.YAxis.FontSize = 24;
ax.XAxis.FontSize = 24;