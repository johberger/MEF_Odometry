%readKittiEgoMotionGT

function [ Result ] = readKittiEgoMotionGT(path,i)

%path = '/home/jberger/Benchmarks/kitti/dataset/poses/00.txt'

Motions = textread(path, '%n','delimiter', '\n');
l = length(Motions);
m = l/12;

MotionsMat = reshape(Motions, 12, m)';

M = zeros(4);
M(4,4) = 1.;
Mold = M;
Mnew = M;
Mold(1:3,1:4) = reshape(MotionsMat(i+1,:),4,3)';
Mnew(1:3,1:4) = reshape(MotionsMat(i+2,:),4,3)';
Mdiff = Mnew\Mold;

R = Mdiff(1:3,1:3)';
h = Mold(1:3,1:3)'*(Mnew(1:3,4)- Mold(1:3,4));

Result = [R,h; [0,0,0],1];

end