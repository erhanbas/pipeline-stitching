classdef ControlPoints
    properties
        Vecfield
        Scopeloc
        Neigs
        Regpts
        pixres
    end
    methods
        function xyz = gridLocations(obj,idx)
            x_lim_cntrl = obj.Vecfield.xlim_cntrl;
            y_lim_cntrl = obj.Vecfield.ylim_cntrl;
            z_lim_cntrl = obj.Vecfield.zlim_cntrl(:,idx);
            [xx,yy,zz] = ndgrid(x_lim_cntrl,y_lim_cntrl,z_lim_cntrl);
            xyz = [xx(:),yy(:),zz(:)];
        end
        function res = residualStage(obj,idx)
            
        end
        function res = residualStageOnZ(obj,idx_center,idx_target)
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            
            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            
            tile_descs_center_um = tile_descs_center.*obj.pixres + tile_loc_um_center;
            tile_descs_target_um = tile_descs_target.*obj.pixres + tile_loc_um_target;
            
            % tile_loc_um_center - tile_loc_um_target
            res = tile_descs_center_um - tile_descs_target_um;
        end
        function res = residualAffineFCOnZ(obj,idx_center,idx_target)
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            aff_center = obj.Vecfield.afftile(:,:,idx_center)/1e3; %in um
            aff_target = obj.Vecfield.afftile(:,:,idx_target)/1e3; %in um
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;

            tile_descs_center_aff_um = aff_center * [tile_descs_center ones(size(tile_descs_center,1),1)]';
            tile_descs_target_aff_um = aff_target * [tile_descs_target ones(size(tile_descs_target,1),1)]';

            % tile_loc_um_center - tile_loc_um_target
            res = (tile_descs_center_aff_um - tile_descs_target_aff_um)';

        end
        
        function r = roundOff(obj)
            r = round([obj.Value],2);
        end
        function r = multiplyBy(obj,n)
            r = [obj.Value] * n;
        end
    end
end