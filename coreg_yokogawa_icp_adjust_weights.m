%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coreg_yokogawa_icp is a function to coregister a structural MRI with MEG data
% and associated headshape information
%
% Written by Robert Seymour Oct/Nov 2017 (some subfunctions written by Paul
% Sowman)
%
% INPUTS:
% - dir_name        = directory name for the output of your coreg
% - confile         = full path to the con file
% - mrkfile         = full path to the mrk file
% - mri_file        = full path to the NIFTI structural MRI file
% - hspfile         = full path to the hsp (polhemus headshape) file
% - elpfile         = full path to the elp file
% - hsp_points      = number of points for downsampling the headshape (try 100-200)
% - scalpthreshold  = threshold for scalp extraction (try 0.05 if unsure)
%
% VARIABLE INPUTS (if using please specify all):
% - do_vids         = save videos to file. Requires CaptureFigVid.
% - weight_number   = how strongly do you want to weight the facial points?
% - bad_coil        = is there a bad coil to take out?
%
% EXAMPLE FUNCTION CALL:
% coreg_yokogawa_icp(dir_name,confile,mrkfile,mri_file,hspfile,elpfile,...
% hsp_points, scalpthreshold,'yes')
%
% OUTPUTS:
% - grad_trans              = correctly aligned sensor layout
% - headshape_downsampled   = downsampled headshape (original variable name I know)
% - mri_realigned           = the mri realigned based on fiducial points
% - trans_matrix            = transformation matrix for accurate coregistration
% - headmodel_singleshell   = coregistered singleshell headmodel
%
% THIS IS A WORK IN PROGRESS FUNCTION - any updates or suggestions would be
% much appreciated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coreg_yokogawa_icp_adjust_weights(dir_name,confile,mrkfile,mri_file,hspfile,elpfile,hsp_points,scalpthreshold,varargin)

if isempty(varargin)
    do_vids = 'no';
    weight_number = 0.1;
    bad_coil = '';
else
    do_vids = varargin{1};
    weight_number = 1./varargin{2};
    bad_coil = varargin{3}
end

cd(dir_name); disp('CDd to the right place');

% Get Polhemus Points
[shape] = parsePolhemus(elpfile,hspfile);

% Read the grads from the con file
grad_con                    = ft_read_sens(confile); %in cm, load grads

% Read the mrk file
if isempty(bad_coil)
    mrk                         = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
    markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
    [R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
    meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
    
    grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
    grad_trans.fid              = shape; %add in the head information
    % Rotate about z-axis
    rot180mat       = rotate_about_z(180);
    grad_trans      = ft_transform_geometry(rot180mat,grad_trans);
    save grad_trans grad_trans
    
    % Else if there is a bad marker
else
    fprintf(''); disp('TAKING OUT BAD MARKER');
    
    % Identify the bad coil
    badcoilpos = find(ismember(shape.fid.label,bad_coil));
    
    % Take away the bad marker
    marker_order = [2 3 1 4 5];
    
    mrk                         = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
    markers                     = mrk.fid.pos(marker_order,:);%reorder mrk to match order in shape
    % Now take out the bad marker when you realign
    markers(badcoilpos-3,:) = [];
    fids_2_use = shape.fid.pnt(4:end,:); fids_2_use(badcoilpos-3,:) = [];
    [R,T,Yf,Err]                = rot3dfit(markers,fids_2_use);%calc rotation transform
    meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
    
    grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
    grad_trans.fid              = shape; %add in the head information
    % Rotate about z-axis
    rot180mat       = rotate_about_z(180);
    grad_trans      = ft_transform_geometry(rot180mat,grad_trans);
    save grad_trans grad_trans
end

% Get headshape downsampled to 100 points with facial info preserved
headshape_downsampled = downsample_headshape(hspfile,hsp_points,grad_trans);
% Rotate about z-axis
headshape_downsampled = ft_transform_geometry(rot180mat,headshape_downsampled);
save headshape_downsampled headshape_downsampled

% Load in MRI
mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';

% Give rough estimate of fiducial points
cfg                         = [];
cfg.method                  = 'interactive';
cfg.viewmode                = 'ortho';
cfg.coordsys                = 'neuromag';
[mri_realigned]             = ft_volumerealign(cfg, mri_orig);

save mri_realigned mri_realigned

% check that the MRI is consistent after realignment
ft_determine_coordsys(mri_realigned, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_headshape(headshape_downsampled);

%% Extract Scalp Surface
cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = scalpthreshold;
scalp  = ft_volumesegment(cfg, mri_realigned);

%% Create mesh out of scalp surface
cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 90000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'cm');
% Flip the mesh around (improves coreg)
%mesh.pos(:,2) = mesh.pos(:,2).*-1;

%% Create Figure for Quality Checking
if strcmp(do_vids,'yes')
    try
        figure;
        ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on; drawnow;
        view(0,0);
        ft_plot_headshape(headshape_downsampled); drawnow;
        OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        CaptureFigVid([0,0; 360,0], 'mesh_quality',OptionZ)
    catch
        disp('You need CaptureFigVid in your MATLAB path. Download at https://goo.gl/Qr7GXb');
        figure;
        ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on; drawnow;
        ft_plot_headshape(headshape_downsampled); drawnow;
        view(0,0);print('mesh_quality','-dpdf');
    end
else
    figure;
    ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; drawnow;
    view(90,0);
    ft_plot_headshape(headshape_downsampled); drawnow;
    title('If this looks weird you might want to adjust the cfg.scalpthreshold value');
    print('mesh_quality','-dpdf');
end

%% Perform ICP using mesh and headshape information
numiter = 50;
disp('Performing ICP fit with 50 iterations');

% Weight the facial points x10 times higher than the head points
count_facialpoints2 = find(headshape_downsampled.pos(:,3)<3); x = 1;

% If there are no facial points ignore the weighting
if isempty(count_facialpoints2)
    w = ones(size(headshape_downsampled.pos,1),1).*1;
    weights = @(x)assignweights(x,w);
    
% But if there are facial points apply weighting using weight_number
else
    w = ones(size(headshape_downsampled.pos,1),1).*weight_number;
    w(count_facialpoints2) = 1;
    weights = @(x)assignweights(x,w);
end

% Now try ICP with weights
[R, t, err, dummy, info] = icp(mesh.pos', headshape_downsampled.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true, 'Weight',weights, 'WorstRejection', 0.05);

%% Create figure to display how the ICP algorithm reduces error
clear plot;
figure; plot([1:1:51]',err,'LineWidth',8);
ylabel('Error'); xlabel('Iteration');
title('Error*Iteration');
set(gca,'FontSize',25);

%% Create transformation matrix
trans_matrix = inv([real(R) real(t);0 0 0 1]);
save trans_matrix trans_matrix

%% Create figure to assess accuracy of coregistration
mesh_spare = mesh;
mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);
c = datestr(clock); %time and date

if strcmp(do_vids,'yes')
    try
        figure;
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on;
        ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
        OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)
    catch
        figure;
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on;
        ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
        print('ICP_quality','-dpdf');
        disp('You need CaptureFigVid in your MATLAB path. Download at https://goo.gl/Qr7GXb');
    end
else
    figure;
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on;
    ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
    clear c; print('ICP_quality','-dpdf');
end

%% Apply transform to the MRI
%mri_realigned = ft_transform_geometry_PFS_hacked(trans_matrix, mri_realigned);

%% Segment
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_realigned);

%% Create singleshell headmodel
cfg = [];
cfg.method='singleshell';

headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, create headmodel

% Flip headmodel around
%headmodel_singleshell.bnd.pos(:,2) = headmodel_singleshell.bnd.pos(:,2).*-1;
% Apply transformation matrix
headmodel_singleshell.bnd.pos = ft_warp_apply(trans_matrix,headmodel_singleshell.bnd.pos);

figure;ft_plot_headshape(headshape_downsampled) %plot headshape
ft_plot_sens(grad_trans, 'style', 'k*')
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 1; camlight
view([90,0]); title('After Coreg');
print('headmodel_quality','-dpdf');
save headmodel_singleshell headmodel_singleshell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function [shape] = parsePolhemus(elpfile,hspfile)
        
        fid1 = fopen(elpfile);
        C = fscanf(fid1,'%c');
        fclose(fid1);
        
        E = regexprep(C,'\r','xx');
        E = regexprep(E,'\t','yy');
        
        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');
        sensornamesi = strfind(E,'%N');
        fiducialsstarti = strfind(E,'%F');
        lastfidendi = strfind(E(fiducialsstarti(3):fiducialsstarti(length(fiducialsstarti))+100),'xx');
        fiducialsendi = fiducialsstarti(1)+strfind(E(fiducialsstarti(1):fiducialsstarti(length(fiducialsstarti))+lastfidendi(1)),'xx');
        
        NASION = E(fiducialsstarti(1)+4:fiducialsendi(1)-2);
        NASION = regexprep(NASION,'yy','\t');
        NASION = str2num(NASION);
        
        LPA = E(fiducialsstarti(2)+4:fiducialsendi(2)-2);
        LPA = regexprep(LPA,'yy','\t');
        LPA = str2num(LPA);
        
        RPA = E(fiducialsstarti(3)+4:fiducialsendi(3)-2);
        RPA = regexprep(RPA,'yy','\t');
        RPA = str2num(RPA);
        
        LPAredstarti = strfind(E,'LPAred');
        LPAredendi = strfind(E(LPAredstarti(1):LPAredstarti(length(LPAredstarti))+45),'xx');
        LPAred = E(LPAredstarti(1)+11:LPAredstarti(1)+LPAredendi(2)-2);
        LPAred = regexprep(LPAred,'yy','\t');
        LPAred = str2num(LPAred);
        
        RPAyelstarti = strfind(E,'RPAyel');
        RPAyelendi = strfind(E(RPAyelstarti(1):RPAyelstarti(length(RPAyelstarti))+45),'xx');
        RPAyel = E(RPAyelstarti(1)+11:RPAyelstarti(1)+RPAyelendi(2)-2);
        RPAyel = regexprep(RPAyel,'yy','\t');
        RPAyel = str2num(RPAyel);
        
        PFbluestarti = strfind(E,'PFblue');
        PFblueendi = strfind(E(PFbluestarti(1):PFbluestarti(length(PFbluestarti))+45),'xx');
        PFblue = E(PFbluestarti(1)+11:PFbluestarti(1)+PFblueendi(2)-2);
        PFblue = regexprep(PFblue,'yy','\t');
        PFblue = str2num(PFblue);
        
        LPFwhstarti = strfind(E,'LPFwh');
        LPFwhendi = strfind(E(LPFwhstarti(1):LPFwhstarti(length(LPFwhstarti))+45),'xx');
        LPFwh = E(LPFwhstarti(1)+11:LPFwhstarti(1)+LPFwhendi(2)-2);
        LPFwh = regexprep(LPFwh,'yy','\t');
        LPFwh = str2num(LPFwh);
        
        RPFblackstarti = strfind(E,'RPFblack');
        RPFblackendi = strfind(E(RPFblackstarti(1):end),'xx');
        RPFblack = E(RPFblackstarti(1)+11:RPFblackstarti(1)+RPFblackendi(2)-2);
        RPFblack = regexprep(RPFblack,'yy','\t');
        RPFblack = str2num(RPFblack);
        
        allfids = [NASION;LPA;RPA;LPAred;RPAyel;PFblue;LPFwh;RPFblack];
        fidslabels = {'NASION';'LPA';'RPA';'LPAred';'RPAyel';'PFblue';'LPFwh';'RPFblack'};
        
        fid2 = fopen(hspfile);
        C = fscanf(fid2,'%c');
        fclose(fid2);
        E = regexprep(C,'\r','xx'); %replace returns with "xx"
        E = regexprep(E,'\t','yy'); %replace tabs with "yy"
        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');
        
        headshapestarti = strfind(E,'position of digitized points');
        headshapestartii = strfind(E(headshapestarti(1):end),'xx');
        headshape = E(headshapestarti(1)+headshapestartii(2)+2:end);
        headshape = regexprep(headshape,'yy','\t');
        headshape = regexprep(headshape,'xx','');
        headshape = str2num(headshape);
        
        shape.pnt = headshape;
        shape.fid.pnt = allfids;
        shape.fid.label = fidslabels;
        
        %convert to BESA style coordinates so can use the .pos file or sensor
        %config from .con
        shape.pnt = cat(2,fliplr(shape.pnt(:,1:2)),shape.pnt(:,3)).*1000;
        %shape.pnt = shape.pnt(1:length(shape.pnt)-15,:); % get rid of nose points may want to alter or comment this depending on your digitisation
        %shape.pnt = shape.pnt*1000;
        neg = shape.pnt(:,2)*-1;
        shape.pnt(:,2) = neg;
        
        shape.fid.pnt = cat(2,fliplr(shape.fid.pnt(:,1:2)),shape.fid.pnt(:,3)).*1000;
        %shape.fid.pnt = shape.fid.pnt*1000;
        neg2 = shape.fid.pnt(:,2)*-1;
        shape.fid.pnt(:,2) = neg2;
        shape.unit='mm';
        shape = ft_convert_units(shape,'cm');
        
        new_name2 = ['shape.mat'];
        save (new_name2,'shape');
    end

    function [R,T,Yf,Err] = rot3dfit(X,Y)
        %ROT3DFIT Determine least-square rigid rotation and translation.
        % [R,T,Yf] = ROT3DFIT(X,Y) permforms a least-square fit for the
        % linear form
        %
        % Y = X*R + T
        %
        % where R is a 3 x 3 orthogonal rotation matrix, T is a 1 x 3
        % translation vector, and X and Y are 3D points sets defined as
        % N x 3 matrices. Yf is the best-fit matrix.
        %
        % See also SVD, NORM.
        %
        % rot3dfit: Frank Evans, NHLBI/NIH, 30 November 2001
        %
        
        % ROT3DFIT uses the method described by K. S. Arun, T. S. Huang,and
        % S. D. Blostein, "Least-Squares Fitting of Two 3-D Point Sets",
        % IEEE Transactions on Pattern Analysis and Machine Intelligence,
        % PAMI-9(5): 698 - 700, 1987.
        %
        % A better theoretical development is found in B. K. P. Horn,
        % H. M. Hilden, and S. Negahdaripour, "Closed-form solution of
        % absolute orientation using orthonormal matrices", Journal of the
        % Optical Society of America A, 5(7): 1127 - 1135, 1988.
        %
        % Special cases, e.g. colinear and coplanar points, are not
        % implemented.
        
        %error(nargchk(2,2,nargin));
        narginchk(2,2); %PFS Change to update
        if size(X,2) ~= 3, error('X must be N x 3'); end;
        if size(Y,2) ~= 3, error('Y must be N x 3'); end;
        if size(X,1) ~= size(Y,1), error('X and Y must be the same size'); end;
        
        % mean correct
        
        Xm = mean(X,1); X1 = X - ones(size(X,1),1)*Xm;
        Ym = mean(Y,1); Y1 = Y - ones(size(Y,1),1)*Ym;
        
        % calculate best rotation using algorithm 12.4.1 from
        % G. H. Golub and C. F. van Loan, "Matrix Computations"
        % 2nd Edition, Baltimore: Johns Hopkins, 1989, p. 582.
        
        XtY = (X1')*Y1;
        [U,S,V] = svd(XtY);
        R = U*(V');
        
        % solve for the translation vector
        
        T = Ym - Xm*R;
        
        % calculate fit points
        
        Yf = X*R + ones(size(X,1),1)*T;
        
        % calculate the error
        
        dY = Y - Yf;
        Err = norm(dY,'fro'); % must use Frobenius norm
    end

    function [output] = ft_transform_geometry_PFS_hacked(transform, input)
        
        % FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to
        % a structure with geometric information, for example a volume conduction model
        % for the head, gradiometer of electrode structure containing EEG or MEG
        % sensor positions and MEG coil orientations, a head shape or a source model.
        %
        % The units in which the transformation matrix is expressed are assumed to
        % be the same units as the units in which the geometric object is
        % expressed. Depending on the input object, the homogeneous transformation
        % matrix should be limited to a rigid-body translation plus rotation
        % (MEG-gradiometer array), or to a rigid-body translation plus rotation
        % plus a global rescaling (volume conductor geometry).
        %
        % Use as
        %   output = ft_transform_geometry(transform, input)
        %
        % See also FT_WARP_APPLY, FT_HEADCOORDINATES
        
        % Copyright (C) 2011, Jan-Mathijs Schoffelen
        %
        % This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
        % for the documentation and details.
        %
        %    FieldTrip is free software: you can redistribute it and/or modify
        %    it under the terms of the GNU General Public License as published by
        %    the Free Software Foundation, either version 3 of the License, or
        %    (at your option) any later version.
        %
        %    FieldTrip is distributed in the hope that it will be useful,
        %    but WITHOUT ANY WARRANTY; without even the implied warranty of
        %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        %    GNU General Public License for more details.
        %
        %    You should have received a copy of the GNU General Public License
        %    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
        %
        % $Id: ft_transform_geometry.m$
        
        % flg rescaling check
        allowscaling = ~ft_senstype(input, 'meg');
        
        % determine the rotation matrix
        rotation = eye(4);
        rotation(1:3,1:3) = transform(1:3,1:3);
        
        if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
            error('invalid transformation matrix');
        end
        
        %%### get rid of this accuracy checking below as some of the transformation
        %%matricies will be a bit hairy###
        if ~allowscaling
            % allow for some numerical imprecision
            %if abs(det(rotation)-1)>1e-6%100*eps
            %if abs(det(rotation)-1)>100*eps  % allow for some numerical imprecision
            %error('only a rigid body transformation without rescaling is allowed');
            %end
        end
        
        if allowscaling
            % FIXME build in a check for uniform rescaling probably do svd or so
            % FIXME insert check for nonuniform scaling, should give an error
        end
        
        tfields   = {'pos' 'pnt' 'o' 'coilpos' 'chanpos' 'chanposold' 'chanposorg' 'elecpos', 'nas', 'lpa', 'rpa', 'zpoint'}; % apply rotation plus translation
        rfields   = {'ori' 'nrm'     'coilori' 'chanori' 'chanoriold' 'chanoriorg'};                                          % only apply rotation
        mfields   = {'transform'};           % plain matrix multiplication
        recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
        % the field 'r' is not included here, because it applies to a volume
        % conductor model, and scaling is not allowed, so r will not change.
        
        fnames    = fieldnames(input);
        for k = 1:numel(fnames)
            if ~isempty(input.(fnames{k}))
                if any(strcmp(fnames{k}, tfields))
                    input.(fnames{k}) = apply(transform, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, rfields))
                    input.(fnames{k}) = apply(rotation, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, mfields))
                    input.(fnames{k}) = transform*input.(fnames{k});
                elseif any(strcmp(fnames{k}, recfields))
                    for j = 1:numel(input.(fnames{k}))
                        input.(fnames{k})(j) = ft_transform_geometry(transform, input.(fnames{k})(j));
                    end
                else
                    % do nothing
                end
            end
        end
        output = input;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogeneous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [new] = apply(transform, old)
        old(:,4) = 1;
        new = old * transform';
        new = new(:,1:3);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotate_about_z - make a rotation matix for arbitrary rotation in degrees
% around z axis
%
% Written by Paul Sowman Oct 2017 (http://web.iitd.ac.in/~hegde/cad/lecture/L6_3dtrans.pdf - page 4)
%
% INPUTS:
% - deg        = degrees of rotation required
%
% OUTPUTS:
% - rmatx      = a 4*4 rotation matrix for deg degrees about z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rmatx=rotate_about_z(deg)
        
        deg   = deg2rad(deg);
        rmatx = [cos(deg) sin(deg) 0 0;-sin(deg) cos(deg) 0 0;0 0 1 0;0 0 0 1];
    end

    function [headshape_downsampled] = downsample_headshape(path_to_headshape,numvertices,sensors)
        % Get headshape
        headshape = ft_read_headshape(path_to_headshape);
        % Convert to cm
        headshape = ft_convert_units(headshape,'cm');
        % Convert to BESA co-ordinates
        headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
        headshape.pos(:,2) = headshape.pos(:,2).*-1;
        
        % Get indices of facial points (up to 4cm above nasion)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 4cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        count_facialpoints = find(headshape.pos(:,3)<4);
        if isempty(count_facialpoints)
            disp('CANNOT FIND ANY FACIAL POINTS - COREG BY ICP MAY BE INACCURATE');
        else
            facialpoints = headshape.pos(count_facialpoints,:,:);
            rrr = 1:4:length(facialpoints);
            facialpoints = facialpoints(rrr,:); clear rrr;
        end
        
        % Remove facial points for now
        headshape.pos(count_facialpoints,:) = [];
        
        % Create mesh out of headshape downsampled to x points specified in the
        % function call
        cfg.numvertices = numvertices;
        cfg.method = 'headshape';
        cfg.headshape = headshape.pos;
        mesh = ft_prepare_mesh(cfg, headshape);
        
        % Replace the headshape info with the mesh points
        headshape.pos = mesh.pos;
        
        % Create figure for quality checking
        figure; subplot(2,2,1);ft_plot_mesh(mesh); hold on;
        title('Downsampled Mesh');
        view(0,0);
        subplot(2,2,2);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 1');
        view(0,0);
        subplot(2,2,3);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 2');
        view(90,0);
        subplot(2,2,4);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 3');
        view(180,0);
        print('headshape_quality','-dpdf');
        
        % Add the facial info back in
        headshape.pos = vertcat(headshape.pos,facialpoints);
        % Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};
        
        % Convert fiducial points to BESA
        headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
        headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;
        
        % Plot for quality checking
        figure;ft_plot_sens(sensors) %plot channel position : between the 1st and 2nd coils
        ft_plot_headshape(headshape) %plot headshape
        view(0,0);
        print('headshape_quality2','-dpdf');
        
        % Export filename
        headshape_downsampled = headshape;
        
    end

    % Assign Weights Function
    function y = assignweights(x, w)
        
        % x is an indexing vector with the same number of arguments as w
        y = w(:)';
    end
end


