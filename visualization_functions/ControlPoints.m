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
        function [res,tile_descs_center_um,tile_descs_target_um] = residualStageOnZ(obj,idx_center,idx_target)
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            
            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            
            tile_descs_center_um = tile_descs_center.*obj.pixres.*[-1 -1 1] + tile_loc_um_center;
            tile_descs_target_um = tile_descs_target.*obj.pixres.*[-1 -1 1] + tile_loc_um_target;
            
            % tile_loc_um_center - tile_loc_um_target
            res = tile_descs_center_um - tile_descs_target_um;
        end
        function [res,tile_descs_center_aff_um,tile_descs_target_aff_um] = residualAffineFCOnZ(obj,idx_center,idx_target)
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            %%
            aff_center = obj.Vecfield.afftile(:,:,idx_center)/1e3; %in um
            aff_target = obj.Vecfield.afftile(:,:,idx_target)/1e3; %in um
            tile_descs_center_aff_um = [tile_descs_center ones(size(tile_descs_center,1),1)]*aff_center';
            tile_descs_target_aff_um = [tile_descs_target ones(size(tile_descs_target,1),1)]*aff_target';
            %% tile_loc_um_center - tile_loc_um_target
            res = (tile_descs_center_aff_um - tile_descs_target_aff_um);
            
            %% estimate affine from point match
            if 0
                desc_dif = tile_descs_center-tile_descs_target;
                scope_dif = tile_loc_um_center-tile_loc_um_target;
                estimated_affine = desc_dif\[ones(size(desc_dif,1),1)*scope_dif];
                res = desc_dif * estimated_affine - scope_dif;
                tile_descs_center_aff_um = tile_descs_center * estimated_affine + tile_descs_center;
                tile_descs_target_aff_um = tile_descs_target * estimated_affine + tile_descs_target;
            elseif 0
                %% estimate affine from point match
                desc_dif = tile_descs_center-tile_descs_target;
                mean_desc_dif = mean(desc_dif);
                scope_dif = tile_loc_um_center-tile_loc_um_target;
                estimated_affine = desc_dif\[ones(size(desc_dif,1),1)*scope_dif];
                tile_descs_center_aff_um = (tile_descs_center-mean(tile_descs_center)) * estimated_affine + tile_descs_center;
                tile_descs_target_aff_um = (tile_descs_target-mean(tile_descs_target)) * estimated_affine + tile_descs_target;
                res = tile_descs_center_aff_um-tile_descs_target_aff_um-mean_desc_dif;
                %%
            end
        end
        function [res, tile_descs_center_ctrl_um, tile_descs_target_ctrl_um] = residualCtrlOnZ(obj,idx_center,idx_target)
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            xyz_center = obj.gridLocations(idx_center);
            xyz_target = obj.gridLocations(idx_target);
            control_points_center = obj.Vecfield.control(:,:,idx_center)/1e3; % in um
            control_points_target = obj.Vecfield.control(:,:,idx_target)/1e3; % in um
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
            tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
            res = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
        end
        
        function um_t = interpolatedInstance(obj,xyz,um,xyz_t)
            FxU = scatteredInterpolant(xyz,um(:,1),'linear','none');
            FxV = scatteredInterpolant(xyz,um(:,2),'linear','none');
            FxW = scatteredInterpolant(xyz,um(:,3),'linear','none');
            um_t = [FxU(xyz_t) FxV(xyz_t) FxW(xyz_t)];
        end
        
        function r = roundOff(obj)
            r = round([obj.Value],2);
        end
        function r = multiplyBy(obj,n)
            r = [obj.Value] * n;
        end
    end
end