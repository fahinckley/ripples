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
bHgtMin = binHeight(1); % minimum beam source height [m]
bHgtRange = len*bS; % range of beam source heights [m]

%% Ejected grains
meanEj = 5; % mean number of ejected grains
stdEj = 2; % standard deviation for ejected grains
dist = 50; % transport distance

%% Main loop
figure
numSteps = 2e5;
for ii = 1:numSteps
    % Beam source height 
    bHgt = bHgtMin + bHgtRange*rand(1);
    
    % Beam height at all x
    beam = bHgt - bS*x;
    
    % Beam intersection with bed
    bInt = find(beam < binHeight,1,'first');
    
    % Ejected grains
    numEj = meanEj + stdEj*randn(1);
    
    % Transport distance
    dTrans = dist;
    
    % Landing bins
    bLand = bInt + 4;
    
    % Check for wrap
    if bLand > numBins
        bLand = bLand - numBins;
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

%% Final plots
hDiff = binHeight-binHeight0;
figure
plot(x,hDiff)

%% Checks
% Conservation of mass
con = sum(hDiff);
fprintf('Sum of height differences: %3.2e',con)
fprintf('\n')

