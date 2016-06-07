% Create the ratio versus latency graph in Henson 2002 and Figure 1 in
% Steffener 2010.
% This plot uses the standard double Gamma hemodynamic response model
% created by default with SPM.
% From this graph one can select any rations of interest.
%
% written by Jason Steffener
% finalized 06/08/2011

dt = 0.05;
[bf p] = spm_hrf(0.05);
p = [6 16 1 1 6 0 32];
%p = [3 10 1 1 6 0 32];
[bf p]         = spm_hrf(0.05, p);
% add time derivative
%---------------------------------------------------------------
dp     = 1;
p(6)   = p(6) + dp;
D      = (bf(:,1) - spm_hrf(dt,p))/dp;
bf     = [bf D(:)];
xBF.bf  =  spm_orth(bf);
BFtime = [0:dt:(length(xBF.bf)-1)*dt];
Hrf = xBF.bf;
SS = Hrf'*Hrf;
xBF.bf(:,1) = Hrf(:,1)./sqrt(SS(1,1));
xBF.bf(:,2) = Hrf(:,2)./sqrt(SS(2,2));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = length(Hrf);
time = [0:dt:(M-1)*dt];
% rescale the Basis functions.
SS = Hrf'*Hrf;
Hrf(:,1) = Hrf(:,1)./sqrt(SS(1,1));
Hrf(:,2) = Hrf(:,2)./sqrt(SS(2,2));

X = xBF.bf;
X(:,1) = X(:,1)./max(X(:,1)).*max(Hrf(:,1));
X(:,2) = X(:,2)./max(X(:,2)).*max(Hrf(:,1));

% create list of potential ratios for the beta values
Ratios = [-4:0.01:-0.5 -0.5:0.001:0.5 0.5:0.01:4];
% preallocate memory
PeakTime = time(find(Hrf(:,1) == max(Hrf(:,1))));
Nr = length(Ratios);
EstT2P = zeros(Nr,1);
nB = zeros(Nr,2);
EstFWHM = zeros(Nr,1);
EstMaxUnder = zeros(Nr,1);
EstTime2Min = zeros(Nr,1);
normTempC = zeros(Nr,1);
PeakValues = zeros(Nr,2);
pHrf = inv(Hrf'*Hrf)*Hrf';
VAF = zeros(Nr,1);
RatioPeak2Min = zeros(Nr,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure of the basis functions
f=figure(1);
clf 
a = axes;
hold on
b = plot(time,Hrf,'k');
set(b(1),'LineWidth',3)
set(b(2),'LineStyle','--','LineWidth',3)

c = plot(BFtime,X,'k');
set(c(1),'LineWidth',1)
set(c(2),'LineStyle','--','LineWidth',1)

grid on
axis([-1 30 -0.13 0.13])
set(a,'YTickLabel',[])
c = xlabel('time (s)')
set(c,'FontSize',24,'FontWeight','b')
c = ylabel('Arbitrary')
set(c,'FontSize',24,'FontWeight','b')
legend('Canonical','Derivative','Location','NorthEast')
legend('Function 1','Function 2','Canonical','Derivative','Location','NorthEast')
set(a,'FontSize',24,'FontWeight','b')
set(f,'Position',[100 200 800 600])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For every ratio calculate the estimated hemodynamic response and then a
% series of parameters. Although only the latency is plotted other
% parameters are calculated.
beta = zeros(2,Nr);
MaxValue= zeros(1,Nr);
MinValue = zeros(1,Nr);
for i = 1:Nr
    % create the HRF from the current ratio
    temp = Hrf(:,1) + Ratios(i).*Hrf(:,2);
    temp = temp./max(temp);
%     plot(time,temp)
%     pause(0.25)
% end
    % Find out how well the basis set fits these estimated HDRs
    beta(:,i) = pHrf*temp;
    fit = Hrf*beta(:,i);
    VAF(i) = var(fit)/var(temp);
    PeakValues(i,:) = [max(Hrf(:,1)) max(temp)];
    % normalize the ratio
    tempC = [1 Ratios(i)];
    nB(i,:) = tempC./norm(tempC);
   % find the time to peak
    EstT2P(i) = time(find(temp == max(temp)));
    MaxValue(i) = temp(find(temp == max(temp)));
    MaxLoc = min(find(temp == max(temp)));
    MinLoc = min(find(temp == min(temp)));
    if MinLoc < MaxLoc
         MinLoc = min(find(temp == min(temp(MaxLoc:end))));
    end
    MinValue(i) = temp(MinLoc);
    % Time to Peak
    Postp2 = (MaxLoc+min(find(temp(MaxLoc:end) < 0.5*temp(MaxLoc))))-1;
    
    Postp1 = Postp2 - 1;
    
    Prep2 = (max(find(temp(1:MaxLoc) < 0.5*temp(MaxLoc))));
    Prep1 = Prep2 + 1;
       
    Slope = (time(Postp2) - time(Postp1))/(temp(Postp2) - temp(Postp1));
    Post = Slope*(0.5*temp(MaxLoc)-temp(Postp1)) + time(Postp1);

    Slope = (time(Prep2) - time(Prep1))/(temp(Prep2) - temp(Prep1));
    Pre = Slope*(0.5*temp(MaxLoc)-temp(Prep1)) + time(Prep1);
     % find the FWHM
    EstFWHM(i) = Post - Pre;
    % Size of Undershoot
    EstMaxUnder(i) = temp(MinLoc);
    % Time of Undershoot
    EstTime2Min(i) = time(MinLoc);
    RatioPeak2Min(i) = abs(MinValue(i));
%      plot(time,temp./max(temp))
%      pause(0.5)
    
end
% find out what the contrast are that correspond to each beta ratio


 
%% %%%%%%%%%%%%%%%%%%%%%%%
% plot the ration versus latency
f = figure(31);
clf
b = axes;
hold on
set(b,'FontSize',24,'FontWeight','b')
a = plot(EstT2P-PeakTime,Ratios,'k')
set(a,'LineWidth',3)
grid
axis([-2.5 2.5 -4 4])
xlabel('Latency Difference, dt (s)');
ylabel({'Derivative:Canonical Ratio', '\beta_2/\beta_1 Ratio'})
set(f,'Position',[100 200 1200 600])
set(b,'YMinorGrid','on')

