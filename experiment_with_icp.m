%% Set paths to necesary files
dir_name = '/Users/44737483/Documents/alien_data/2704/visual'; cd(dir_name);
confile = '/Users/44737483/Documents/mcq_data/2704/meg/run-alien/2704_AT_ME160_2017_09_27_aliens.con';
mrkfile = '/Users/44737483/Documents/mcq_data/2704/meg/run-alien/2704_AT_ME160_2017_09_27_aliens_PRE.mrk';
mri_file = '/Users/44737483/Documents/mcq_data/2704/anat/2704.nii';
elpfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.elp';
hspfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.hsp';

coreg_yokogawa_icp_adjust_weights(dir_name,confile,mrkfile,mri_file,hspfile,elpfile,100, 0.04,'yes',10);

cd(dir_name); disp('CDd to the right place');

% Get Polhemus Points
[shape] = parsePolhemus(elpfile,hspfile);

% Read the grads from the con file
grad_con                    = ft_read_sens(confile); %in cm, load grads

% Read the mrk file
mrk                         = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
[R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix

% Transform sensors based on the MRKfile
grad_trans      = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
grad_trans.fid  = shape; %add in the head information
% Rotate about z-axis
rot180mat       = rotate_about_z(180);
grad_trans      = ft_transform_geometry(rot180mat,grad_trans);
%save grad_trans grad_trans

% Get headshape downsampled to 100 points with facial info preserved
headshape_downsampled = downsample_headshape(hspfile,300,grad_trans);
% Rotate about z-axis
headshape_downsampled = ft_transform_geometry(rot180mat,headshape_downsampled);
%save headshape_downsampled headshape_downsampled

% Load in MRI
mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';

% Give rough estimate of fiducial points
cfg                         = [];
cfg.method                  = 'interactive';
cfg.viewmode                = 'ortho';
cfg.coordsys                = 'neuromag';
[mri_realigned]             = ft_volumerealign(cfg, mri_orig);

%save mri_realigned mri_realigned

% check that the MRI is consistent after realignment
ft_determine_coordsys(mri_realigned, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_headshape(headshape_downsampled);

%% Extract Scalp Surface
cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = 0.04;
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
figure;
ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on; drawnow;
view(0,0);
ft_plot_headshape(headshape_downsampled); drawnow;
%title('If this looks weird you might want to adjust the cfg.scalpthreshold value');
%print('mesh_quality','-dpdf');
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true; 
CaptureFigVid([0,0; 360,0], 'mesh_quality',OptionZ)

%% Perform ICP using mesh and headshape information
numiter = 50;
disp('Performing ICP fit with 50 iterations');

x = 1;
w = ones(size(headshape_downsampled.pos,1),1).*0.1;
count_facialpoints = find(headshape_downsampled.pos(:,3)<3);
w(count_facialpoints) = 1;
weights = @(x)assignweights(x,w);

[R, t, err, dummy, info] = icp(mesh.pos', headshape_downsampled.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true, 'Weight',weights, 'WorstRejection', 0.05);

clear plot;
figure; plot([1:1:51]',err,'LineWidth',8);
ylabel('Error'); xlabel('Iteration');
title('Error*Iteration');
set(gca,'FontSize',25);

%% Create transformation matrix
trans_matrix = inv([real(R) real(t);0 0 0 1]);
%save trans_matrix trans_matrix

%% Create figure to assess accuracy of coregistration
mesh_spare = mesh;
mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);

figure;
ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
c = datestr(clock); %time and date
camlight; hold on; view([-180,0]);
ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, err(end)));
clear c; OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true; 
CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)

