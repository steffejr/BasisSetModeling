function subfnCreateTimeToPeakMaps(Image1, Image2, CondName, SPM)
% This is a simple function which creates a time to peak graph based on two
% basis functions. This can easily be modified for use with three basis
% functions. It is based off of SPM5 analyses.
%
% written 7/9/2010
% Jason Steffener
%

% load the headers
    V1 = spm_vol(deblank(Image1));
    V2 = spm_vol(deblank(Image2));
    % load images
    I1 = spm_read_vols(V1);
    I2 = spm_read_vols(V2);

    % load the basis functions
    Hrf = SPM.xBF.bf(:,1:2);
    dt = SPM.xBF.dt;
    Nbf = length(SPM.xBF.bf);
    % create a time vector
    time = [0:dt:(Nbf-1)*dt];
   F = find(~isnan(I1));
   N = length(F);
   OutData = zeros(size(I1));
   for j = 1:N
       % find voxel specific HRF
       temp = [I1(F(j)) I2(F(j))];
       if temp(1) > 0
           tempHRF = Hrf*temp';
           % find max or min
           FT2P = find(max(tempHRF) == tempHRF);
           T2P = time(FT2P);
           OutData(F(j)) = T2P;

       end
   end
   % Write the time to peak map to disk
   Vo = V1;
   Vo.descrip = 'time to peak plot for beta1 > 0';
   Vo.fname = [sprintf('Time2Peak2bf_%s.nii',CondName)];
   spm_write_vol(Vo, OutData);
   
