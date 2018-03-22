%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coreg_elelta_icp_adjust_weights is a function to coregister a structural 
% MRI with MEG data and associated headshape information
%
% Written by Robert Seymour, March 2018
%
% INPUTS:
% - dir_name        = directory name for the output of your coreg
% - fif_file         = full path to the fif file
% - mri_file        = full path to the NIFTI structural MRI file
% - scalpthreshold  = threshold for scalp extraction (try 0.05 if unsure)
%
% VARIABLE INPUTS (if using please specify all):
% - do_vids         = save videos to file. Requires CaptureFigVid.
% - weight_number   = how strongly do you want to weight the facial points?
%                     Input a number from 0 (no weighting) to 1 (full
%                     weighting)
%
% EXAMPLE FUNCTION CALL:
% coreg_elekta_icp_adjust_weights(dir_name,fif_file,mri_file,...
% scalpthreshold,'yes',0.8)
%
% OUTPUTS:
% - sens                    = correctly aligned sensor layout
% - mri_realigned           = the mri realigned based on fiducial points
% - trans_matrix            = transformation matrix for accurate coregistration
% - mri_realigned2          = the coregistered mri based on ICP algorithm
% - headmodel_singleshell   = coregistered singleshell headmodel
%
% THIS IS A WORK IN PROGRESS FUNCTION - any updates or suggestions would be
% much appreciated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coreg_elekta_icp_adjust_weights(dir_name,fif_file,mri_file,scalpthreshold,varargin)

% If user doesn't specify varargin, use defaults
if isempty(varargin)
    do_vids = 'no';
    weight_number = 0.1;
else
    do_vids = varargin{1};
    weight_number = 1./varargin{2};
end

cd(dir_name); disp('CDd to the right place');

% Get Polhemus Points
[headshape] = ft_read_headshape(fif_file);

% Read the grads from the con file
sens  = ft_read_sens(fif_file); %in cm, load grads
save sens sens

% Load in MRI
mri_orig = ft_read_mri(mri_file); % in mm, read in mri from DICOM
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
ft_plot_headshape(headshape);

%% Extract Scalp Surface
cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = scalpthreshold;
scalp  = ft_volumesegment(cfg, mri_realigned);

%% Create mesh out of scalp surface
cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 75000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'cm');

%% Create Figure for Quality Checking
if strcmp(do_vids,'yes')
    try
        figure;
        ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor',...
            'none','facealpha',0.8); hold on;
        camlight; hold on; drawnow;
        view(0,0);
        ft_plot_headshape(headshape); drawnow;
        OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        CaptureFigVid([0,0; 360,0], 'mesh_quality',OptionZ)
    catch
        disp('You need CaptureFigVid in your MATLAB path. Download at https://goo.gl/Qr7GXb');
        figure;
        ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor',...
            'none','facealpha',0.8); hold on;
        camlight; hold on; drawnow;
        ft_plot_headshape(headshape); drawnow;
        view(0,0);print('mesh_quality','-dpdf');
    end
else
    figure;
    ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; drawnow;
    view(90,0);
    ft_plot_headshape(headshape); drawnow;
    title('If this looks weird you might want to adjust the cfg.scalpthreshold value');
    print('mesh_quality','-dpdf');
end

%% Perform ICP using mesh and headshape information
numiter = 50;
disp('Performing ICP fit with 50 iterations');

% Weight the facial points x10 times higher than the head points
count_facialpoints2 = find(headshape.pos(:,3)<3); x = 1;

% If there are no facial points ignore the weighting
if isempty(count_facialpoints2)
    w = ones(size(headshape.pos,1),1).*1;
    weights = @(x)assignweights(x,w);
    
% But if there are facial points apply weighting using weight_number
else
    w = ones(size(headshape.pos,1),1).*weight_number;
    w(count_facialpoints2) = 1;
    weights = @(x)assignweights(x,w);
end

% Now try ICP with weights
[R, t, err] = icp(mesh.pos', headshape.pos', numiter, 'Minimize', 'plane',...
    'Extrapolation', true, 'Weight', weights, 'WorstRejection', 0.05);

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
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,...
            'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on;
        ft_plot_headshape(headshape); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
        OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)
    catch
        figure;
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on;
        ft_plot_headshape(headshape); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
        print('ICP_quality','-dpdf');
        disp('You need CaptureFigVid in your MATLAB path. Download at https://goo.gl/Qr7GXb');
    end
else
    figure;
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor',...
        'none','facealpha',0.8); hold on;
    camlight; hold on;
    ft_plot_headshape(headshape); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
    clear c; print('ICP_quality','-dpdf');
end

%% Apply transform to the MRI
mri_realigned2 = ft_transform_geometry(trans_matrix,mri_realigned);

save mri_realigned2 mri_realigned2

% check that the MRI is consistent after realignment
ft_determine_coordsys(mri_realigned2, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_headshape(headshape);
drawnow;

%% Segment
fprintf(' Segmenting the brain ... could take a while\n');
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_realigned2);

%% Create singleshell headmodel
cfg = [];
cfg.method='singleshell';

headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, create headmodel

figure;ft_plot_headshape(headshape) %plot headshape
ft_plot_sens(sens, 'style', 'k*')
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor',...
    'cortex'); alpha(1.0); hold on;
ft_plot_mesh(mesh_spare,'facecolor','skin'); alpha(0.2); camlight;
view([90,0]); title('After Coreg');
print('headmodel_quality','-dpdf');
save headmodel_singleshell headmodel_singleshell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign Weights Function
    function y = assignweights(x, w)
        
        % x is an indexing vector with the same number of arguments as w
        y = w(:)';
    end

end


