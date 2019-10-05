clc; close all; clear variables; format short g

%% Parameters
tho = 2;	% Time Constant of Entropy (Parameter)
w = 0.5;    % Relative Cost of Transition to Similar State
a = .05;    % Significance Level (Alpha)
b = 2;      % Configurations of Variable
b0 = 0;     % ID of Insignificant Bit
D = 7;    % Temporal Length of Variable 2:9
C = 1;      % Entropy Threshold of Filtering 0.1:0.1:2

%% Analysing Range of Configurations
nd = length(D); nc = length(C);
N_Pass = zeros(nc,nd); N_All = zeros(nc,nd);
%% Creating All Possible Scenarios
for di = 1:nd
    d = D(di);
    n0 = b^d;
    X = dec2base(0:n0-1,b,4)-'0';
    
    %% Removing Non-treatment Scenarios
    Keep = ones(n0,1);
    for i = 1:n0
        x = X(i,:);
        if sum(x==1)==0
            Keep(i) = 0;
        end
    end
    X = X(find(Keep),:);
    n = size(X,1);
    
    %% Cataloging Low Entropy Configurations
    X(1,end+1) = 0;
    for j = 1:n
        x = X(j,1:end-1);
        x = x(find((x==b0) == 0,1):end);
        X(j,end) = Entropy(x,tho,w);
    end
    
    %% Filtering Hypotheses Based on Threshold
    X = sortrows(X,d+1);
    for ci = 1:nc
        c = C(ci);
        if c >= X(end,end)
            Cut = n;
        elseif c <= X(1,end)
            Cut = 0;
        else
            Cut = find(X(:,end)>c,1)-1;
        end
        X1 = X(1:Cut,:); X0 = X(Cut+1:end,:); n1 = size(X1,1);      
        N_Pass(ci,di) = n1; N_All(ci,di) = n0;
    end
end

%% Min Achievable P-value From Each Sample Size
P_Power = ones(30,1);
for Samples = 1:30
    Z = sqrt(3*Samples*(Samples+1)/(4*Samples+2));
    P_Power(Samples) = 2*normcdf(-Z);
end

%% Filtering Hypotheses Based on Sample Size
N_Pass1 = N_Pass; N_Pass1(N_Pass1==0) = 1; P_Needed = a./N_Pass1;
M0 = zeros(size(N_Pass));
for i = 1:numel(P_Needed)
    M0(i) = find(P_Power<=P_Needed(i),1);
end

%% Display Results
figure; plot(P_Power,'.-','LineWidth',1,'MarkerSize',14); grid minor;
xlabel('Sample Size'); ylabel('Min Achievable P-value');
if numel(C)>1 && numel(D)>1
figure; surf(D,C,N_Pass+1); Plot=gca; set(Plot,'zscale','log'); grid minor;
title('Scaling of passed Hypotheses'); zlabel('Number of Hypothesis');
xlabel('Temporal Depth'); ylabel('Entropy Threshold');
figure; surf(D,C,N_Pass); grid minor;
title('Number of passed Hypotheses'); zlabel('Number of Hypothesis');
xlabel('Temporal Depth'); ylabel('Entropy Threshold');
figure; surf(D,C,N_Pass./N_All); grid minor;
title('Percentage of passed Hypotheses'); zlabel('Number of Hypothesis');
xlabel('Temporal Depth'); ylabel('Entropy Threshold');
figure; surf(D,C,M0); grid minor;
title('Min Sample Size to Pass FWER'); zlabel('Number of Samples');
xlabel('Temporal Depth'); ylabel('Entropy Threshold');
else
fprintf([repmat('%.4g  ',1,d+1) '\n'],X1');
disp('===========================')
fprintf([repmat('%.4g  ',1,d+1) '\n'],X0');
end
% figure; subplot(1,2,1);
% histogram(X(:,end),'Normalization','cdf','BinWidth',.02);
% grid minor; set(gca,'YTick',0:0.05:1,'XTick',0:.2:2.6);
% title('Cumulative Distribution of Entropy');
% xlabel('Permutation Entropy'); ylabel('Number of Possible Hypotheses');
% subplot(1,2,2);
% histogram(X(:,end),'BinWidth',.1);
% grid minor; set(gca,'XTick',0:.2:2.6);
% title('Distribution of Entropy');
% xlabel('Permutation Entropy'); ylabel('Percentage of Possible Hyp');

%% Function: Modified Permutation Entropy @@ changed for b = 2
function [H] = Entropy(x,tho,w)
V = [];
for i = 0:tho-1
    V = [V; x(i+1:length(x)) zeros(1,i)];
end
V = V(:,1:length(x)-tho+1)';

% [~,~,Temp] = unique(V,'rows','stable');
% f = accumarray(Temp,1); f(f==0) = []; f = f/sum(f);
% H = - f'*log2(f);
f = [sum(ismember(V,[1 1],'rows'));sum(ismember(V,[0 0],'rows'))
    sum(ismember(V,[1 0],'rows'));sum(ismember(V,[0 1],'rows'))];
H = 0; f = f/sum(f);
for i = 1:4
    if f(i)>0
        if i<=2; H = H - w*f(i)*log2(f(i));
        else; H = H - f(i)*log2(f(i)); end
    end
end

end