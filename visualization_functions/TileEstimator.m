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
        
        function [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onx(obj,ix,iy)
            xyz_center = obj.gridLocations(ix);
            control_points_center = obj.Vecfield.control(:,:,ix)/1e3; % in um
            xyz_target = obj.gridLocations(iy);
            control_points_target = obj.Vecfield.control(:,:,iy)/1e3; % in um
            tile_descs_center = obj.Paireddescriptor{ix}.onx.X;
            tile_descs_target = obj.Paireddescriptor{ix}.onx.Y; % this is on ix not iy !!
        end
        
        function [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_ony(obj,ix,iy)
            xyz_center = obj.gridLocations(ix);
            control_points_center = obj.Vecfield.control(:,:,ix)/1e3; % in um
            xyz_target = obj.gridLocations(iy);
            control_points_target = obj.Vecfield.control(:,:,iy)/1e3; % in um
            tile_descs_center = obj.Paireddescriptor{ix}.ony.X;
            tile_descs_target = obj.Paireddescriptor{ix}.ony.Y; % this is on ix not iy !!
        end
        
        
        function [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onz(obj,ix,iy)
            xyz_center = obj.gridLocations(ix);
            control_points_center = obj.Vecfield.control(:,:,ix)/1e3; % in um
            xyz_target = obj.gridLocations(iy);
            control_points_target = obj.Vecfield.control(:,:,iy)/1e3; % in um
            tile_descs_center = obj.Regpts{ix}.X;
            tile_descs_target = obj.Regpts{ix}.Y; % this is on ix not iy !!
        end
        function X = finitefilter(obj,X)
            X = X(all(isfinite(X),2),:);
        end
        function mX = meanSqrt(obj,X)
            mX = mean(sqrt(sum(X.^2,2)),'omitnan');
        end
        
        function [residual,residual_onx, residual_ony, residual_onz, stats] = estimateResidual4ctrl(obj,idx_center,idx_target)
            
            %% estimates residual error for matched descriptors
            % gets a tile, and its adjacent tiles and estimates residual
            if nargin == 3 % estimate only for target
                if isnan(idx_target); residual = nan;return;end
                % on z
                xyz_center = obj.gridLocations(idx_center);
                control_points_center = obj.Vecfield.control(:,:,idx_center)/1e3; % in um
                
                xyz_target = obj.gridLocations(idx_target);
                control_points_target = obj.Vecfield.control(:,:,idx_target)/1e3; % in um
                
                tile_descs_center = obj.Regpts{idx_center}.X;
                tile_descs_target = obj.Regpts{idx_center}.Y;
                
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
                return
            end
            
            neigs = obj.Neigs(idx_center,:); %[id -x -y +x +y -z +z]
            [residual_onxp1,residual_onxm1,residual_onyp1,residual_onym1,residual_onzp1,residual_onzm1] = deal([]);
            % on +x
            ix = idx_center;
            iy = neigs(4);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onx(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onxp1 = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            
            % on -x
            ix = neigs(2);
            iy = neigs(1);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onx(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onxm1 = -(tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            % on +y
            ix = idx_center;
            iy = neigs(5);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_ony(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onyp1 = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            
            % on -y
            ix = neigs(3);
            iy = neigs(1);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_ony(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onym1 = -(tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            
            % on +z
            ix = neigs(1);
            iy = neigs(7);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onz(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onzp1 = (tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            
            % on -z
            ix = neigs(6);
            iy = neigs(1);
            if all(isfinite([ix,iy]))
                [xyz_center,xyz_target,control_points_center,control_points_target,tile_descs_center,tile_descs_target] = getPoints_onz(obj,ix,iy);
                tile_descs_center_ctrl_um = obj.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
                tile_descs_target_ctrl_um = obj.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
                residual_onzm1 = -(tile_descs_center_ctrl_um - tile_descs_target_ctrl_um);
            end
            
            %             residual_list.onxp1 = residual_onxp1;
            %             residual_list.onxm1 = residual_onxm1;
            %             residual_list.onyp1 = residual_onyp1;
            %             residual_list.onym1 = residual_onym1;
            %             residual_list.onzp1 = residual_onzp1;
            %             residual_list.onzm1 = residual_onzm1;
            %[id -x -y +x +y -z +z]
            %%
            residual = [
                residual_onxp1;residual_onxm1;
                residual_onyp1;residual_onym1;
                residual_onzp1;residual_onzm1;
                ];
            residual = obj.finitefilter(residual);
            
            residual = obj.finitefilter(residual);
            residual_onxp1 = obj.finitefilter(residual_onxp1);
            residual_onxm1 = obj.finitefilter(residual_onxm1);
            residual_onyp1 = obj.finitefilter(residual_onyp1);
            residual_onym1 = obj.finitefilter(residual_onym1);
            residual_onzp1 = obj.finitefilter(residual_onzp1);
            residual_onzm1 = obj.finitefilter(residual_onzm1);
            
            residual_onx = [residual_onxp1;residual_onxm1];
            residual_ony = [residual_onyp1;residual_onym1];
            residual_onz = [residual_onzp1;residual_onzm1];
            residual_onx = residual_onx(all(isfinite(residual_onx),2),:);
            residual_ony = residual_ony(all(isfinite(residual_ony),2),:);
            residual_onz = residual_onz(all(isfinite(residual_onz),2),:);
            
            stats = [obj.meanSqrt(residual) obj.meanSqrt(residual_onxp1) obj.meanSqrt(residual_onxm1) obj.meanSqrt(residual_onyp1) obj.meanSqrt(residual_onym1) obj.meanSqrt(residual_onzp1) obj.meanSqrt(residual_onzm1)];
            
            
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
