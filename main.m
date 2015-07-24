addpath('../matlab/MinimumEnergyFilter/MinimumEnergyFilter_Flow_SE_comparison');

seq = 0;
i=1;
numFix = 4;

start = 0;
ending = 200;

dimx1 = 376;
dimx2 = 1241;

x1 = kron((1:dimx2)', ones(dimx1,1));
x2 = kron(ones(dimx2,1), (1:dimx1)');
XY = horzcat(x1,x2); % dimensions of image

% vectors with rotational and translational error, respectively
error_rot = zeros(ending-start,1);
error_trans = zeros(ending-start,1);
error_rot_orig = zeros(ending-start,1);
error_trans_orig = zeros(ending-start,1);

% read calibration file
CalibPath = sprintf('/home/jberger/Benchmarks/kitti/odometry/sequences/%02d/calib.txt',seq);
[Camera, baseline] = readKittiCalib(CalibPath); % read calib file
Motionpath = ['/home/jberger/Benchmarks/kitti/odometry/poses/', sprintf('%02d',seq), '.txt'];

n_points = 10;
q = 1.;
s1 = 1e-2;
s2 = 1e-5;
alpha = 2;
numSteps = 200;
delta = 1/numSteps;
order = 1;

W = diag([s1,s1,s1,s2,s2,s2]');
E_orig = eye(4);
Q = diag([q/n_points,q/n_points]');
Q1st = eye(6);

MEF_mex('init', Camera, q, s1, s2, alpha, numSteps, order);

for i=start:ending
    % flow and depth from Vogel et al
    flow = readFlowFile(sprintf('/mnt/flow/kitti/PRSM/%02d/flow/%06d.flo',seq, i+1));
    depth = h5read(['/mnt/flow/kitti/PRSM/', sprintf('%02d',seq), '/depth/', sprintf('%06d',i+1) ,'.h5'],'/depth');
    Depth_vec = depth(:);

    % compute observation from optical flow:
    flowu = flow(:,:,1);
    flowu = flowu(:);
    flowv = flow(:,:,2);
    flowv = flowv(:);
    flowUV = horzcat(flowu, flowv);
    
    %extract Sample Points
    SampleX = 120 + randsample(dimx1-200,n_points, true);
    SampleY = 100 + randsample(dimx2-2*100,n_points, true);
    %% Find corresponding pixel positions to SampleX and Sample Y
    Sample = dimx1*(SampleY-1) + SampleX;

    XYs = XY(Sample,:);
    XYZs = horzcat(XYs, ones(length(XYs),1));
    Depths = Depth_vec(Sample);
    flowUVs = flowUV(Sample,:);

    XYZhom = (Camera\XYZs')';

    flowTest = horzcat(flowUVs, zeros(length(flowUVs),1));
    flowHom =  (Camera\flowTest')';

    y = XYZhom + flowHom; % observations
    g = (Depths*ones(1,4)).*horzcat(XYZhom ,1./Depths);
    gy = horzcat(g,y);

    %% call MinimumEnergyFilter
    E = MEF_mex('iterate', XYs, Depths, flowUVs);
    
    %% Comparison
    for j=1:numSteps
      % Implicit midpoint rule;
      grad1st = computeRiemannianGradient(E_orig,gy,Q,1); 
      KE1st = zeros(4,4);
      for p=1:numFix
        KE1st = delta * dynamicsEFlow(Q1st, computeRiemannianGradient(E_orig*expm(0.5*KE1st),gy,Q,1),E_orig,1);
      end
      Hess1st = computeHessian1(E_orig,gy,Q,1);
      E_orig =    E_orig * expm(KE1st);
      Q1st = integrateRiccati(Q1st, Hess1st, alpha, W, grad1st, E_orig, delta, 1);
    end
    
    E_orig
    Q1st
    
    
    E_gt = readKittiEgoMotionGT(Motionpath,i);
    
    error = errorE(E,E_gt);
    error_orig = errorE(E_orig, E_gt);
    error_rot(i+1,1) = error(2);
    error_trans(i+1,1) = error(1);
    error_rot_orig(i+1,1) = error_orig(2);
    error_trans_orig(i+1,1) = error_orig(1);
    
    figure(1);
    semilogy(   (start+1):ending, error_rot,'r-x', ...
                (start+1):ending, error_rot_orig,'b-x');
    title('Rotational error');
    drawnow

    figure(2);
    semilogy((start+1):ending, error_trans,'r-x', ...
        (start+1):ending, error_trans_orig,'b-x');
    title('Translational error');
    drawnow
    
end


