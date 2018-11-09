classdef resStats < handle
    properties
        Ctrl_res
        Aff_res
        Stg_res
        Aff
        Stg
        Residual_onx
        Residual_ony
        Residual_onz
        numtile
    end
    methods
        function init(obj)
            obj.numtile = length(obj.Ctrl_res);
        end
%         function estimateStats(obj)
%         end
        function estimateMSE(obj)
            [S_mse,A_mse,C_mse] = deal(nan(1,obj.numtile));
            for it = 1:obj.numtile
                if ~isnan(obj.Aff{it}); S_mse(it) = mean(sqrt(sum(obj.Stg_res{it}.^2,2))); end
                if ~isnan(obj.Aff{it}); A_mse(it) = mean(sqrt(sum(obj.Aff_res{it}.^2,2))); end
                if ~isnan(obj.Aff{it}); C_mse(it) = mean(sqrt(sum(obj.Ctrl_res{it}.^2,2)),'omitnan'); end
            end
            
            
        end
        function histmod = getMod(obj,inhist)
            [aa,bb] = max(inhist.Values);
            histmod = mean(inhist.BinEdges(bb:bb+1));
        end
        function hist_mse(obj,fignum,C_mse,A_mse,S_mse,finterior)
            figure(fignum), 
            clf
            ax = gca;
            binsp = .025;
            h_ctrl = histogram(C_mse(finterior),'BinWidth',binsp);
            C_st = [obj.getMod(h_ctrl) median(C_mse(finterior),'omitnan') mean(C_mse(finterior),'omitnan') std(C_mse(finterior),'omitnan')];
            hold on
            leg1 = sprintf('Control points x: %1.3f',C_st(1));
            legs{1} = leg1;
            % leg1 = sprintf('Control points: %1.3f \\pm %1.3f',C_st(1),C_st(2));
            if ~isempty(A_mse)
                h_aff = histogram(A_mse(finterior),'BinWidth',binsp);
                A_st = [obj.getMod(h_aff) median(A_mse(finterior),'omitnan') mean(A_mse(finterior),'omitnan') std(A_mse(finterior),'omitnan')];
                leg2 = sprintf('Control points y: %1.3f',A_st(1));
                % leg2 = sprintf('Affine: %1.3f \\pm %1.3f',A_st(1),A_st(2));
                legs{end+1} = leg2;
            end
            if ~isempty(S_mse)
                h_st = histogram(S_mse(finterior),'BinWidth',binsp);
                S_st = [obj.getMod(h_st) median(S_mse(finterior),'omitnan') mean(S_mse(finterior),'omitnan') std(S_mse(finterior),'omitnan')];
                leg3 = sprintf('Control points z: %1.3f',S_st(1));
                %leg3 = sprintf('Stage: %1.3f \\pm %1.3f',S_st(1),S_st(2));
                legs{end+1} = leg3;
            end
            
            leg = legend(legs);
            % legend('test \pm')
            leg.FontSize = 32;
            xlim([0 5])
            
            title('Residual out of 11657 tiles', 'FontSize', 30)
            xlabel('Residual magnitude in \mum', 'FontSize', 30)
            ax.YAxis.FontSize = 24;
            ax.XAxis.FontSize = 24;
        end
        function boxplt(obj,fignum,C_mse,A_mse,S_mse,finterior)
            len = length(finterior);
            C_fint = C_mse(finterior);
            S_fint = S_mse(finterior);
            A_fint = A_mse(finterior);
            X_fint = [C_fint,A_fint,S_fint];
            
            Origin = [repmat({'Ctrl'},len,1);repmat({'Affine'},len,1);repmat({'Stage'},len,1)];
            leg1 = sprintf('%1.3f \\pm %1.3f',C_st(1),C_st(2));
            leg2 = sprintf('%1.3f \\pm %1.3f',A_st(1),A_st(2));
            leg3 = sprintf('%1.3f \\pm %1.3f',S_st(1),S_st(2));
            
            figure(fignum),
            clf,
            ax = cla;
            hold on
            boxplot(X_fint,Origin)
            text(1-.1,median(C_fint,'omitnan'), leg1,'FontSize',12,'FontWeight','bold')
            text(2-.1,median(A_fint,'omitnan'), leg2,'FontSize',12,'FontWeight','bold')
            text(3-.1,median(S_fint,'omitnan'), leg3,'FontSize',12,'FontWeight','bold')
            
            ylim([-1 30])
            ylabel('Residual in \mum')
            
            ax.YAxis.FontSize = 34;
            ax.XAxis.FontSize = 44;
            set(ax,'Ytick',0:5:32)
        end
    end
end
