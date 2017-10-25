%% Set paths to necesary files
dir_name = '/Users/44737483/Documents/alien_data/2660/visual'; cd(dir_name);
confile = '/Users/44737483/Documents/mcq_data/2660/meg/run-alien/2660_EA_ME160_2017_08_28_aliens.con';
mrkfile = '/Users/44737483/Documents/mcq_data/2660/meg/run-alien/2660_EA_ME160_2017_08_28_alien_PRE.mrk';
niifile = '/Users/44737483/Documents/mcq_data/2660/anat/2660.nii';

% I think we should change the parsePolhemus function to accept .elp 
% files directly
fname = ('/Users/44737483/Documents/mcq_data/2660/meg/2660_EA_ME160_2017_08_28'); 
parsePolhemus(fname);

% Again this part could be made more elegant
load shape; % in mm, the headshape and fiducial data from the polhemus 
shape                       = ft_convert_units(shape,'cm');
save shape shape

% Read the grads
grad_con                    = ft_read_sens(confile); %in cm, load grads

% Read the mrk file
mrk                         = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
[R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix

grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
grad_trans.fid              = shape; %add in the head information

%%plot 
figure;ft_plot_sens(grad_trans) %plot channel position : between the 1st and 2nd coils
ft_plot_headshape(grad_trans.fid) %plot headshape
ft_plot_sens(grad_con,'style', 'sk','edgecolor','none','facecolor',[0.9 0.9 0.9]) %plot coil position 

save grad_trans grad_trans

load('shape.mat'); load('grad_trans.mat');

mri_file = niifile; %name of the MRI in nift format
mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
% Changing to cm here messes things up for me - so stay in mm
mri_orig.coordsys = 'neuromag';
mri_orig = ft_convert_units(mri_orig,'cm');
% This flips the MRI to improve coreg
%mri_orig.transform(2,:) = mri_orig.transform(2,:).*-1;
save mri_orig mri_orig


%% Segment the brain and create a single-shell headmodel
cfg           = []; 
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_orig);

% Create singleshell headmodel
cfg = [];
cfg.method='singleshell';
headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, creat headmodel

%% Do manual coreg for initial

cfg                         = [];
cfg.method                  = 'interactive';
cfg.viewmode                = 'ortho';
cfg.coordsys                = 'neuromag';
[mri_realigned]             = ft_volumerealign(cfg, mri_orig); 

mri_realigned = ft_convert_units(mri_realigned,'cm');

%align MRI into polhemus coords using the same method
nas                         = mri_realigned.cfg.fiducial.nas;
lpa                         = mri_realigned.cfg.fiducial.lpa;
rpa                         = mri_realigned.cfg.fiducial.rpa;
transm                      = mri_realigned.transform;

nas                         = ft_warp_apply(transm,nas, 'homogenous');
lpa                         = ft_warp_apply(transm,lpa, 'homogenous');
rpa                         = ft_warp_apply(transm,rpa, 'homogenous');


%%tranform cooridnates into the polhemus coordinate
fids_mir                    = [nas;lpa;rpa];
[R,T,Yf,Err]                = rot3dfit(fids_mir,shape.fid.pnt(1:3,:));%calc rotation transform
mri2head_transm             = [[R;T]'; 0 0 0 1];

mri_realigned =  ft_transform_geometry_PFS(mri2head_transm, mri_realigned);

%% Extract the scalp and create a mesh

cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 3;
cfg.scalpthreshold = 0.05;

scalp  = ft_volumesegment(cfg, mri_realigned);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 50000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'cm');

%%
grad_trans.fid = ft_convert_units(grad_trans.fid,'cm');
grad_trans.fid.pos = grad_trans.fid.pnt;

%% Create a figure to check the scalp mesh and headshape points
figure;ft_plot_mesh(mesh,'facealpha',0.2); hold on;
ft_plot_headshape(grad_trans.fid);

%% Perform ICP algorithm
numiter = 50;
w = ones(size(grad_trans.fid.pos,1),1);
%weights = @(x)assignweights(x,w);
[R, t, err, dummy, info] = icp(mesh.pos', grad_trans.fid.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true, 'WorstRejection', 0.05);

%% Create transformation matrix
trans_matrix = inv([R t;0 0 0 1]);
save trans_matrix trans_matrix

%% Create diagram to show accuracy of ICP co-reg
mesh_spare = mesh;
mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);

figure;ft_plot_mesh(mesh_spare,'facealpha',0.2); hold on;
ft_plot_headshape(grad_trans.fid);

%% Align MRI into polhemus coords using the same method
nas                         = mri_realigned.cfg.fiducial.nas;
lpa                         = mri_realigned.cfg.fiducial.lpa;
rpa                         = mri_realigned.cfg.fiducial.rpa;
transm                      = mri_realigned.transform;

nas                         = ft_warp_apply(transm,nas, 'homogenous');
lpa                         = ft_warp_apply(transm,lpa, 'homogenous');
rpa                         = ft_warp_apply(transm,rpa, 'homogenous');

% tranform cooridnates into the polhemus coordinate
fids_mir                    = [nas;lpa;rpa];
[R,T,Yf,Err]                = rot3dfit(fids_mir,shape.fid.pnt(1:3,:));%calc rotation transform
mri2head_transm             = [[R;T]'; 0 0 0 1];

headmodel_singleshell_realigned =  ft_transform_geometry_PFS(mri2head_transm, headmodel_singleshell_realigned);

%% Use ft_warp_apply and trans_matrix to realign headmodel_singleshell 
headmodel_singleshell_realigned = headmodel_singleshell;
headmodel_singleshell_realigned.bnd.pos = ft_warp_apply(trans_matrix, headmodel_singleshell_realigned.bnd.pos);


%% Plot the result of this
figure('units','normalized','outerposition',[0 0 1 1]);
ft_plot_headshape(grad_trans.fid) %plot headshape
%ft_plot_sens(grad_trans,'coil' ,'true','style', 'sk','edgecolor','none','facecolor',[0.6 0.6 0.8]) %plot coil position 
ft_plot_sens(grad_trans, 'style', 'k*')
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 1; camlight
% legend ('realigned sensor position','headpoints polhemus','realigned coil position','aligned headmodel')
view([90,0]); title('Before Coreg'); hold on;

figure('units','normalized','outerposition',[0 0 1 1]);
ft_plot_headshape(grad_trans.fid) %plot headshape
%ft_plot_sens(grad_trans,'coil' ,'true','style', 'sk','edgecolor','none','facecolor',[0.6 0.6 0.8]) %plot coil position 
ft_plot_sens(grad_trans, 'style', 'k*')
ft_plot_vol(headmodel_singleshell_realigned,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 1; camlight
% legend ('realigned sensor position','headpoints polhemus','realigned coil position','aligned headmodel')
ft_plot_mesh(mesh_spare,'facealpha',0.2); hold on;
view([90,0]); title('After Coreg'); 

save headmodel_singleshell_realigned  headmodel_singleshell_realigned



