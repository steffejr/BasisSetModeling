function [ncR ncL] = subfnCalcCont(Design, ratio, ColumnsOfInterest);
% Name: subfnCalcCont
% Inputs:
%           Design: design matrix
%           ratio: the ratio of interest between the two basis functions
%           ColumnsOfInterest: The two columns for which the contrast is
%               estimated for. The first column is expected to correspond to
%               the primary basis function.
%
% This program calculates the contrast weights required to combine two
% basis functions at the given ratio. The program checks the design matrix
% to determine if it is normalized.
%
% Written by: Jason Steffener
% js2746@columbia.edu
% date: November 3, 2009
%
[N M] = size(Design);
if M < 2
    error('The design should have at least two columns')
end
if length(ColumnsOfInterest) ~= 2
    error('The size of the ColumnsOfInterest vector should be of length 2')
end

% Create contrast vector
c = [1 ratio];
% normalize the contrast
c_norm = c./norm(c);
% Contrast transform matrices
Left =  [0 1; -1 0];
Right = [0 -1; 1 0];
cL = c_norm*Left;
cR = c_norm*Right;
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
% Calculate the contrast weights in the un-normalized design matrix.
% If in fact the design matrix is normalized then the contrast is
% multiplied by an identity matrix.
ncR = cR*SSDesign;
ncL = cL*SSDesign;

