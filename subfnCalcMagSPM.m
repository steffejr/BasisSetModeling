function subfnCalcMagSPM(Image1,Image2,Design, ColumnsOfInterest);
% Name: subfnCalcMagSPM
% Inputs:
%           Image1: full path to beta image corresponding to the primary
%              basis function
%           Image2: full path to beta image corresponding to the secondary
%              basis function
%           Design: design matrix, e.g. SPM.xX.X which is the unfiltered
%               design matrix used stored in the SPM.mat file
%           ColumnsOfInterest: The two columns for which the contrast is
%               estimated for. The first column is expected to correspond to
%               the primary basis function.
%
% This program creates and executes the SPM command lines to calculate the magnitude
% image. It checks the design matrix to see if it is normalized.
%
% Written by: Jason Steffener
% js2746@columbia.edu
% date: November 3, 2009
% Bug Fixes:
% fix output file extensions
% updated: July, 2, 2010
%
% Create a structure of image volumne structures
Vin = spm_vol(Image1);
Vin(2) = spm_vol(Image2);
[PathName FileName Ext] = fileparts(Image1);
% Create the output image structure
%
% The extension has been hardcoded here allowing inputs to be header files,
% or files where the extension is: .img,1 as per the output from the
% spm_select program. --JS 7/2/10
Ext = '.img';
%
VOut = Vin(1);
VOut.fname = fullfile(PathName, ['MagnitudeImage' Ext]);
VOut.descrip = 'Magnitude Image';


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
% Create and execute the spm image calculation. 
Weight1 = SSDesign(1,1);
Weight2 = SSDesign(2,2);
f = ['sqrt(i1.*i1.*' num2str(Weight1) ' + i2.*i2.*' num2str(Weight2) ')']
spm_imcalc(Vin,VOut,f);