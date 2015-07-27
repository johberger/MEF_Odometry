clear all;
dbstop if warning;
addpath('./src/matlab');
addpath('./src/matlab/external/flow_code/');

% SET PARAMETERS FOR THE FILTER
%================================================
n_points = 1000;          % number of points in image domain for which filter is evaluated
q = 2.;                 % weight on the data term (penalty term of abservation  noise)
s1 = 1e-2;              % penalty term of the model noise (rotational components)
s2 = 1e-5;              % penalty term of the model noise (translational components)
alpha = 2;              % decay rate (alpha>0), for alpha = 0 all termporal infomration is taken into account
numSteps = 50;          % number of discretization steps for numerical integration
delta = 1/numSteps;     % integration step size
order = 4;              % order of kinematic model, 1: constant velocity, 2: constant acceleration, 3: constant jerk, etc.
%================================================


% read calibration file
[Camera, baseline] = readKittiCalib(sprintf('data/kitti/00/calib.txt'));
% set path for ground truth motion
path_poses = 'data/kitti/00/poses.txt';

start = 0;      % first frame
stop = 20;      % last frame

% Image dimensions:
dimx1 = 376;    % height
dimx2 = 1241;   % width

x1 = kron((1:dimx2)', ones(dimx1,1));
x2 = kron(ones(dimx2,1), (1:dimx1)');
XY = horzcat(x1,x2); % dimensions of image

% vectors with rotational and translational error, respectively
error_rot = zeros(stop-start,1);
error_trans = zeros(stop-start,1);


% INITIALIZAITON OF THE FILTER
MEF_mex('init', Camera, q, s1, s2, alpha, numSteps, order);

for i=start:stop-1

    % flow and depth computed by the method from Vogel et al.
    flow = readFlowFile(sprintf('data/kitti/00/flow/%06d.flo', i+1));
    depth = h5read(sprintf('data/kitti/00/depth/%06d.h5',i+1),'/depth');
    Depth_vec = depth(:);
    % compute observation from optical flow:
    flowu = flow(:,:,1); flowu = flowu(:);
    flowv = flow(:,:,2); flowv = flowv(:);
    flowUV = horzcat(flowu, flowv);
    
    % extract Sample Points
    % Note, we omit points from the outer regtions
    SampleX = 120 + randsample(dimx1-200,n_points, true);
    SampleY = 100 + randsample(dimx2-2*100,n_points, true);
    Sample = dimx1*(SampleY-1) + SampleX;

    % Generate required data
    XYs = XY(Sample,:);
    XYZs = horzcat(XYs, ones(length(XYs),1));
    % get depth map
    Depths = Depth_vec(Sample);
    % get optical flow
    flowUVs = flowUV(Sample,:);


    % 	ITERATION OF THE FILTER
    E = MEF_mex('iterate', XYs, Depths, flowUVs);
    
    % load ground truth ego motion
    E_gt = readKittiEgoMotionGT(path_poses,i);
    error = errorE(E,E_gt);
    error_rot(i+1,1) = error(2);
    error_trans(i+1,1) = error(1);
    
    % PLOT ERRORS REGARGINH GROUND TRUTH
    figure(1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    % plot left image
    subplot(3,2,1);
    imshow(imread(sprintf('data/kitti/00/image_0/%06d.png',i)));
    title('image of left camera');    
    %plot of right image
    subplot(3,2,2);
    imshow(imread(sprintf('data/kitti/00/image_0/%06d.png',i)));
    title('image of right camera');
    % plot flow
    subplot(3,2,3);
    imshow(flowToColor(flow,100));
    title('Induced flow computed by method from Vogel et al. (2015)');
    % plot depth
    subplot(3,2,4);
    imshow(imagesc2rgb(depth, [0,50]));
    title('Induced flow computed by method from Vogel et al. (2015)');    
    %plot rotational error
    subplot(3,2,5);
    plot((start+1):stop, error_rot,'r-x');
    title('Rotational error in deg');
    xlabel('frames');
    ylabel('rot. err. (deg)');
    % plot translational error
    subplot(3,2,6);
    plot((start+1):stop, error_trans,'r-x');
    title('Translational error in m');
    xlabel('frames');
    ylabel('rot. err. (m)');
    
    drawnow;
end


