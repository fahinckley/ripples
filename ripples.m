%--------------------------------
% Modelling of Aeolian Ripple Formation
%--------------------------------
% Franklin Hinckley
% 31 March 2016
%--------------------------------
%
%--------------------------------

%% Clean up 
clearvars
close all
clc

%% Parameters
D = 1e-3;   % grain diameter [m]

%% Initial bed shape
% Bin widths 
numGr = 10;
binWid = numGr*D;

% Distance array
len = 10;
x = 0.5*binWid : binWid : len;
numBins = length(x);

% Function for bin midpoint heights
binHeight = 2 + 0.01*sin(10*x);

% Save initial shape
binHeight0 = binHeight;

%% Beam
% Parameters
beta = 10; % beam angle [deg]
bS = tand(beta); % beam slope
bHgtMin = min(binHeight); % minimum beam source height [m]
bHgtRange = len*bS; % range of beam source heights [m]

%% Ejected grains
meanEj = 5; % mean number of ejected grains
stdEj = 2; % standard deviation for ejected grains
dist = 4; % transport distance [bins]

%% Main loop
figure
numSteps = 6e4;
for ii = 1:numSteps
    % Beam source height 
    bHgt = bHgtMin + bHgtRange*rand(1);
    
    % Beam height at all x
    beam = bHgt - bS*x;
    
    % Beam intersection with bed
    bInt = find(beam < binHeight,1,'first');
    if isempty(bInt)
        bInt = numBins;
    end
    
    % Ejected grains
    numEj = meanEj + stdEj*randn(1);
    
    % Transport distance
    %dTrans = dist;
    
    % Parabola for grain trajectory
    xPara = bInt : 1 : numBins;
    xPara = [xPara 1:1:bInt-1]; % wrap bins
    hPara = -D * (xPara - (bInt + 2)).^2 + (binHeight(bInt) + 0.01); % y = (x-h)^2 + k
    
    % Find parabola intersections with bed
    bPara = find(hPara < binHeight(xPara),1,'first');
    
    % Landing bins
    %bLand = bInt + dTrans;
    bLand = bInt + bPara;
    
    % Check for wrap
    if bLand > numBins
        bLand = bLand - numBins;
    end
    
    % Error trap
    if isempty(bLand)
        error('No landing bin')
    end
    
    % Update impacted bin
    binHeight(bInt) = binHeight(bInt) - (numEj/numGr)*D;
    
    % Update landing bins
    binHeight(bLand) = binHeight(bLand) + (numEj/numGr)*D;
    
    % Plot
    if mod(ii,100) == 0
        plot(x,binHeight0,'--k')
        hold on
        plot(x,binHeight,'-k')
        plot(x,beam,'--r')
        hold off   
        xlim([0 len])
        ylim([1.9 2.25])
        drawnow
    end
    
end

%% FFTs
% Parameters
Fs = numBins/len;
T = 1/Fs;
L = numBins;
t = (0:L-1)*T;

% FFT of original bins
Y_0 = fft(binHeight0-mean(binHeight0));
P2_0 = abs(Y_0/L);
P1_0 = P2_0(1:L/2+1);
P1_0(2:end-1) = 2*P1_0(2:end-1);
xf_0 = Fs*(0:(L/2))/L;

% FFT of final bins
Y = fft(binHeight-mean(binHeight));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
xf = Fs*(0:(L/2))/L;

%% Final plots
% Height differences
hDiff = binHeight-binHeight0;
figure
plot(x,hDiff)
title('Change in Bed Height')
xlabel('Position [m]')
ylabel('\Delta H [m]') 

% FFTs
figure
hold on
plot(xf,P1_0)
plot(xf,P1)
xlim([0 50])
%ylim([0 0.025])
legend('Original Bed','Final Bed')
title('FFT of Bed Heights')
xlabel('f (Hz)')
ylabel('|P1(f)|') 

%% Checks
% Conservation of mass
con = sum(hDiff);
fprintf('Sum of height differences: %3.2e',con)
fprintf('\n')

