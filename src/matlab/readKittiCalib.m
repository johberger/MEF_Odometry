function [ C, b ] = readKittiCalib(path)
% read the Kitti calibration files and returns
% C - (internal) camera matrix
% b - baseline between stereo cameras

fileID = fopen(path);
A = textscan(fileID, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
B = cell2mat(A(2:13));

C = reshape(B(1,:),4,3)';
C = C(1:3,1:3);

b = B(2,4);

fclose(fileID);


end