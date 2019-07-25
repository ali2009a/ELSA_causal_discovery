clc; clear variables; format short g
tho = 2; % Time Constant of Entropy (Parameter)
b = 3;   % Configurations of Variable
d = 3;   % Temporal Length of Variable

%% Cataloging Low Entropy Configurations
n = b^d;
X = dec2base(0:n-1,3,4)-'0'; X(1,end+1) = 0;
for j = 1:n
    x = X(j,1:end-1);
    X(j,end) = Entropy(x,tho);
end
X = sortrows(X,d+1);

%% Removing Non-treatment Scenarios
Keep = ones(n,1);
for i = 1:n
    x = X(i,1:end-1);
    if sum(x==1)==0
        Keep(i) = 0;
    end
end
X = X(find(Keep),:);

%% Display Results
clc; close all
figure; subplot(1,2,1);
histogram(X(:,end),'Normalization','cdf','BinWidth',.02);
grid minor; set(gca,'YTick',0:0.05:1,'XTick',0:.2:2.6);
title('Cumulative Distribution of Entropy');
xlabel('Permutation Entropy'); ylabel('Number of Possible Hypotheses');
subplot(1,2,2);
histogram(X(:,end),'BinWidth',.1);
grid minor; set(gca,'XTick',0:.2:2.6);
title('Distribution of Entropy');
xlabel('Permutation Entropy'); ylabel('Percentage of Possible Hypotheses');

Threshold = find(X(:,end)>1.8,1)-1;
X1 = X(1:Threshold,:); X0 = X(Threshold+1:end,:); 
fprintf([repmat('%.4g  ',1,size(X1,2)) '\n'],X1');
disp('===========================')
fprintf([repmat('%.4g  ',1,size(X0,2)) '\n'],X0');

%% Modified Permutation Entropy
function [H] = Entropy(x,tho)
V = [];
for i = 0:tho-1
   V = [V; x(i+1:length(x)) zeros(1,i)];
end
V = V(:,1:length(x)-tho+1)';
% Keep = ones(size(V,1),1);
% for i = 1:size(V,1)
%     if numel(unique(V(i,:)))==1
%         Keep(i) = 0;
%     end
% end
% V = V(find(Keep),:);
[~,~,Temp] = unique(V,'rows','stable');
f = accumarray(Temp,1); f(f==0) = []; f = f/sum(f);
H = - f'*log2(f);
end