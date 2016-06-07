function subfnCalcMagFSL(Image1,Image2,Design,ColumnsOfInterest);
% Name: subfnCalcMagFSL
% Inputs:
%           Image1: full path to beta image corresponding to the primary
%              basis function
%           Image2: full path to beta image corresponding to the secondary
%              basis function
%           Design: design matrix
%           ColumnsOfInterest: The two columns for which the contrast is
%               estimated for. The first column is expected to correspond to
%               the primary basis function.
%
% This program creates the FSL command lines to calculate the magnitude
% image. It checks the design matrix to see if it is normalized.
%
% Written by: Jason Steffener
% js2746@columbia.edu
% date: November 3, 2009
%

[PathName FileName Ext] = fileparts(Image1);

% Check Design matrix to determine if it is normalized
tol = 0.0001;
X = Design(:,ColumnsOfInterest);
SSDesign = X'*X;
%SSDesign = nX'*nX;
if sum(diag(SSDesign) - [1; 1]) < tol
    % Design is normalized
else
    nX(:,1) = X(:,1)./norm(X(:,1));
    nX(:,2) = X(:,2)./norm(X(:,2));
end
if SSDesign(2) > tol % Check off diagonal
    error('Design matrix is not orthogonalized');
end

Weight1 = SSDesign(1,1);
Weight2 = SSDesign(2,2);

line1 = ['fslmaths ' Image1 ' -mul ' Image1 ' -mul ' num2str(Weight1) ' ' fullfile(PathName,'MagTemp1')];
line2 = ['fslmaths ' Image2 ' -mul ' Image2 ' -mul ' num2str(Weight2) ' ' fullfile(PathName,'MagTemp2')]
line3 = ['fslmaths ' fullfile(PathName,'MagTemp1') ' -add ' fullfile(PathName,'MagTemp2') ' -sqrt ' fullfile(PathName, 'MagnitudeFSL')]
line4 = ['rm MagTemp1* MagTemp2*'];
fprintf('\n%s\n%s\n%s\n%s\n\n',line1,line2,line3,line4)