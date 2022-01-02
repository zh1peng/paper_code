% sphere_testing
% Test sphere before extraction
% It caculates percentage of a sphere coverage
% zhipeng 2018/12/4
function results=sphere_testing(vols, mnis,r)
%% important:
% as r specify the mm radiaus, covert xyz to mm
% change T if neccessary!
% mnis input as : nx3
T = ...
    [-3     0     0    81;...
    0     3     0  -115;...
    0     0     3   -53;...
    0     0     0     1];
for voli=1:length(vols)
    disp(['Running testing for vol:' ,num2str(voli)]);
    voldata=spm_read_vols(spm_vol(vols{voli}));
    [X,Y,Z]=ind2sub(size(voldata),1:numel(voldata));
    XYZ=cor2mni([X; Y; Z]',T);
    for mni_i=1:size(mnis,1)
        %     fprintf('mni %s /t',mni_i)
        M =XYZ - repmat(mnis(mni_i,:),numel(voldata),1);
        Q=find(sum(M.^2,2) <= r^2);
        results(mni_i, voli)=sum(voldata(Q))./length(Q);
    end
end
end

function mni = cor2mni(cor, T)
% function mni = cor2mni(cor, T)
% convert matrix coordinate to mni coordinate
%
% cor: an Nx3 matrix
% T: (optional) rotation matrix
% mni is the returned coordinate in mni space
%
% caution: if T is not given, the default T is
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
% last revised: 2005-04-30

if nargin == 1
    T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

cor = round(cor);
mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
mni = mni';
mni(:,4) = [];
end

