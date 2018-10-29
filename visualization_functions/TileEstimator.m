classdef TileEstimator
    properties
        Vecfield
        Scopeloc
        Neigs
        Regpts
        pixres = [1 1 1]
        Paireddescriptor
        Scopeparams
    end
    methods
        function [Aest,residual] = estimateStage(obj,idx)
            % t1 + A1*X1 = t2 + A2*X2, where 
            % t : stage coordinate
            % A: affine of a tile
            % X1: descriptors in a tile
            % assume, A1~=A2 for adjacent tiles => A = (t1-t2)/(X2-X1)
            [stage_dif,desc_dif] = obj.getDifferences(idx);
            if isempty(stage_dif)
                [Aest,residual] = deal(nan);
                return
            end
            % scale with pixres and add flip to x&y
            A = diag(obj.pixres.*[-1 -1 1]);
            residual = ( stage_dif + desc_dif * A);
            Aest = A;
        end
        function [Aest,residual] = estimateAffine(obj,idx)
            % t1 + A1*X1 = t2 + A2*X2, where 
            % t : stage coordinate
            % A: affine of a tile
            % X1: descriptors in a tile
            % assume, A1~=A2 for adjacent tiles => A = (t1-t2)/(X2-X1)
            [stage_dif,desc_dif] = obj.getDifferences(idx);
            if isempty(stage_dif)
                [Aest,residual] = deal(nan);
                return
            end
            % solve for affine and add translation
            A = -desc_dif\stage_dif;
            residual = ( stage_dif + desc_dif * A);
            Aest = A;
        end
        function [residual] = estimateResidual4ctrl(obj,idx_center,idx_target)
            % gets a tile, and its adjacent tiles and estimates residual
            if nargin < 3
                idx_target = obj.Neigs(idx_center,end);
            end
            if isnan(idx_target); residual = nan;return;end
            neigs = obj.Neigs(idx_center,:); %[id -x -y +x +y -z +z]
            xyz_center = obj.gridLocations(idx_center);
            xyz_target = obj.gridLocations(idx_target);
            control_points_center = obj.Vecfield.control(:,:,idx_center)/1e3; % in um
            control_points_target = obj.Vecfield.control(:,:,idx_target)/1e3; % in um
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
            tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
            residual = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
%             mean(sqrt(sum(res.^2,2)),'omitnan')

%             for ineigs = 2:legth(neigs)
%                 xyz_center
%                 
%             end
            
            
        end
        %% UTILITY functions
        function um_t = interpolatedInstance(obj,xyz,um,xyz_t)
            FxU = scatteredInterpolant(xyz,um(:,1),'linear','none');
            FxV = scatteredInterpolant(xyz,um(:,2),'linear','none');
            FxW = scatteredInterpolant(xyz,um(:,3),'linear','none');
            um_t = [FxU(xyz_t) FxV(xyz_t) FxW(xyz_t)];
        end
        
        function xyz = gridLocations(obj,idx)
            x_lim_cntrl = obj.Vecfield.xlim_cntrl;
            y_lim_cntrl = obj.Vecfield.ylim_cntrl;
            z_lim_cntrl = obj.Vecfield.zlim_cntrl(:,idx);
            [xx,yy,zz] = ndgrid(x_lim_cntrl,y_lim_cntrl,z_lim_cntrl);
            xyz = [xx(:),yy(:),zz(:)];
        end
        function [accumulator_stage,accumulator_descriptor] = getDifferences(obj,idx)
            % accumulate 6 adjacency for center tile, skip edge cases for
            % simplicty as they are mostly gelatine embedding
            
            neigs = obj.Neigs(idx,:); %[id -x -y +x +y -z +z]
            % descriptor match on z
            [accumulator_stage,accumulator_descriptor] = deal([]);
            %% tile match with above
            if any(isnan(neigs))
                return
            end
            idx_center = neigs(end-1);
            idx_target = idx;
            [stage_dif,des_dif] = getDifferences_onz(obj,idx_center,idx_target);
            accumulator_stage{1} = -stage_dif; % -z direction
            accumulator_descriptor{1} = -des_dif; % -z direction

            %% tile match with below
            idx_center = idx;
            idx_target = neigs(end);

            [stage_dif,des_dif] = obj.getDifferences_onz(idx_center,idx_target);
            accumulator_stage{end+1} = stage_dif; % +z
            accumulator_descriptor{end+1} = des_dif; % +z
            
            %% tile match with left
            % to make it consistent with z calculation, first flip
            % center/target, then negate the numbers
            idx_center = neigs(2);
            idx_target = idx;
            [stage_dif,des_dif] = obj.getDifferences_onx(idx_center,idx_target);
            accumulator_stage{end+1} = -stage_dif; % -x direction
            accumulator_descriptor{end+1} = -des_dif; % -x direction

            %% tile match with right
            idx_center = idx;
            idx_target = neigs(4);
            [stage_dif,des_dif] = obj.getDifferences_onx(idx_center,idx_target);
            accumulator_stage{end+1} = stage_dif; % +x direction
            accumulator_descriptor{end+1} = des_dif; % +x direction
            
            %% tile match with top
            idx_center = neigs(3);
            idx_target = idx;
            [stage_dif,des_dif] = obj.getDifferences_ony(idx_center,idx_target);
            accumulator_stage{end+1} = -stage_dif; % -y direction
            accumulator_descriptor{end+1} = -des_dif; % -y direction

            %% tile match with top
            idx_center = idx;
            idx_target = neigs(5);
            [stage_dif,des_dif] = obj.getDifferences_ony(idx_center,idx_target);
            accumulator_stage{end+1} = stage_dif; % +y direction
            accumulator_descriptor{end+1} = des_dif; % +y direction
            
            %% flatten array
            accumulator_stage = cat(1,accumulator_stage{:});
            accumulator_descriptor = cat(1,accumulator_descriptor{:});
            
        end
        function desc = correctTiles(obj,desc,dims)
            % flip dims
            if nargin<3
                dims = [1024 1536 251];
            end
            desc(:,1:2) = dims(1:2)+1 - desc(:,1:2);
        end
        function [stage_dif,des_dif] = getDifferences_onz(obj,idx_center,idx_target)
            if obj.Regpts{idx_center}.matchrate == 0
                [stage_dif,des_dif] = deal([]);
            end
            tile_descs_center = obj.Regpts{idx_center}.X;
            tile_descs_target = obj.Regpts{idx_center}.Y;
            des_dif = tile_descs_center - tile_descs_target;

            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            stage_dif = tile_loc_um_center - tile_loc_um_target;
            stage_dif = ones(size(des_dif,1),1)*stage_dif;
        end
        function [stage_dif,des_dif] = getDifferences_onx(obj,idx_center,idx_target)
            tile_descs_center = obj.Paireddescriptor{idx_center}.onx.X;% Y, not X, as we get the neighbor in -x
            tile_descs_target = obj.Paireddescriptor{idx_center}.onx.Y;
            % correct dims with flip (our scope does reverse imaging)
            % des_dif = obj.correctTiles(tile_descs_center) - obj.correctTiles(tile_descs_target);
            des_dif = tile_descs_center - tile_descs_target;
            % this is same as tile_descs_target-tile_descs_center
            
            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            stage_dif = tile_loc_um_center - tile_loc_um_target;
            stage_dif = ones(size(des_dif,1),1)*stage_dif;
        end
        function [stage_dif,des_dif] = getDifferences_ony(obj,idx_center,idx_target)
            tile_descs_center = obj.Paireddescriptor{idx_center}.ony.X;% Y, not X, as we get the neighbor in -x
            tile_descs_target = obj.Paireddescriptor{idx_center}.ony.Y;
            % correct dims with flip (our scope does reverse imaging)
            % des_dif = obj.correctTiles(tile_descs_center) - obj.correctTiles(tile_descs_target);
            % this is same as tile_descs_target-tile_descs_center
            des_dif = tile_descs_center - tile_descs_target;
            
            tile_loc_um_center = obj.Scopeloc.loc(idx_center,:)*1e3; %in um
            tile_loc_um_target = obj.Scopeloc.loc(idx_target,:)*1e3; %in um
            stage_dif = tile_loc_um_center - tile_loc_um_target;
            stage_dif = ones(size(des_dif,1),1)*stage_dif;
        end
    end
    
end
