P = spm_select(1,'SPM','Select the SPM.mat file');
load(P)
% Which session?
NSess = length(SPM.Sess);
if NSess > 1
    Sess = input(sprintf('Which session? (1 - %d) ',NSess));
else
    Sess = 1;
end
% How many conditions are there?
NCond = length(SPM.Sess(Sess).U);
% Which condition from this session?
CondSelectionString = sprintf('%s\n','Which condition?');
for i = 1:NCond
    CondSelectionString = sprintf('%s\t%d: %s\n',CondSelectionString,i,SPM.Sess(Sess).U(i).name{1});
end
Cond = input(CondSelectionString);
CondName = SPM.Sess(Sess).U(Cond).name{1};

fprintf(1,'You selected %s\n',SPM.Sess(Sess).U(Cond).name{1});
% How many basis functions are there?
Nbases = length(SPM.Sess(Sess).Fc(Cond).i);
fprintf(1,'There are %d basis functions for this condition\n',Nbases);
% Currently only the first two basis function will be used
% The corresponding columns are
BasesCol = SPM.Sess(Sess).Fc(Cond).i(1:2);
% Find the beta images for theses columns
Image1 = fullfile(SPM.swd,SPM.Vbeta(BasesCol(1)).fname);
Image2 = fullfile(SPM.swd,SPM.Vbeta(BasesCol(2)).fname);

% Get the design matrix
Design = SPM.xX.X;
ColumnsOfInterest = SPM.Sess(Sess).Fc(Cond).i(1:2);
subfnCalcMagSPM(Image1, Image2, Design, ColumnsOfInterest, CondName);

% Create time to peak maps
subfnCreateTimeToPeakMaps(Image1, Image2, CondName, SPM)