clear;
close all; 

% A = {
%     % 'tk0056-171014-04';
%     % 'tk0056-171015-04';
%     % 'tk0062-171129-04';
%     % 'tk0062-171130-04';
%     'tk0064-180206-04';
%     'tk0064-180207-04';
%     'tk0067-180426-04';
%     'tk0067-180427-04';
%     'tk0068-180530-04';
%     'tk0068-180531-04';
%     'tk0069-180726-04';
%     'tk0069-180727-04';
%     'tk0070-180829-04';
%     'tk0070-180830-04';
%     'tk0072-181025-04';
%     'tk0072-181026-04';
%     'tk0074-181219-04';
%     'tk0074-181220-04';
%     'tk0075-190326-04';
%     'tk0075-190327-04';
%     'tk0076-190424-04';
%     'tk0076-190425-04'
%     };
% 
% addpath('posbytrials')
% savePath = 'testGauss';
% 
% B = strcat(A,'rawData'); 

% for k = 1:size(B,1)

    % session = B{k};
    % load([session '.mat'])

    addpath('preprocess')
    load('tk0070-180829-04')

%% Preprocess head-direction data 
    % angleDouble = zeros(length(headDirection),1);
    % for i = 1:length(headDirection)
    %     if headDirection(i)*2<360
    %         angleDouble(i) = headDirection(i)*2;
    %     else
    %         angleDouble(i) = headDirection(i)*2-360;
    %     end 
    % end 
    % 
    % radDouble = angleDouble*2*pi/360 - pi; 
    % 
    
    %% Define behavior by basis functions (von Mises) 
    kappa = 400;                        % concentration parameter (larger the value, smaller the variance and closer to the mean)  
    n = 75;                         % number of von Mises functions 
    % mu = linspace(0,2,n);
    mu = linspace(0,2*pi-2*pi/n,n); % assumes the mean position on the circle increases at a constant rate  
    K = [];
    for i=1:length(mu)  
        % K(:,i) = exp(-(rpos-mu(i)).^2/k.^2/2)/k/sqrt(2*pi);
        K(:,i) = exp(kappa*cos(rpos*pi-mu(i)));  % rpos is position from the range [0,2], mutiply by pi to change the range to [0,2pi]
        K(:,i) = K(:,i)/exp(kappa);                % normalize (minimize) values to reduce large numbers  
    end
    
    rbs_basis = [];
    pvec = linspace(0,2*pi,256);
    % pvec = linspace(0,2,256);
    for i=1:length(mu)
        % rbs_basis(:,i) = exp(-(pvec-mu(i)).^2/k.^2/2)/k/sqrt(2*pi);
        rbs_basis(:,i) = exp(kappa*cos(pvec-mu(i)));
        rbs_basis(:,i) = rbs_basis(:,i)/exp(kappa);
    end
    
    
    %% Split the data into training and testing 
    numel = size(Xf,2); 
    nfold = 10;
    splits = floor(linspace(0,numel,nfold+1));
    ridx=1:numel;
    for i=1:nfold
        cv.train{i} = ridx([1:splits(i) (splits(i+1)+1):numel]);
        cv.test{i} = ridx([(splits(i)+1):splits(i+1)]);
    end
    
    %% Train data and predict behavioral variables 
    y = zscore([real(Xf.') imag(Xf.')]);
    
    f=[];
    for i = 1:nfold
        g = y(cv.train{i},:)\K(cv.train{i},:);              % compute weights on n number of von Mises basis functions 
        f(cv.test{i},:) = y(cv.test{i},:)*(g*rbs_basis');   
    end
    
    f = zscore(f')';
    idx=1:size(f,1); 
    [tmp,maxpost]=max(f'); 
    
    
    %% Calculate circular distance between actual and predicted data
    errp = circ_dist(rpos(idx),maxpost/length(pvec)*2*pi)*360/pi;
    errm = circ_dist(2*pi-rpos(idx),maxpost/length(pvec)*2*pi)*360/pi;
    
    [err_nd_All,mini] = min([abs(errm) abs(errp)]');
    
        
    %% Bootstat median error in degrees with s.e.
    err_mean = [];
    err_std = [];
    
    bootstat = bootstrp(500,'median',abs(err_nd_All));
    err_mean = mean(bootstat);
    err_std = std(bootstat);
    
    % [fi,xi] = ksdensity(bootstat);
    % plot(xi,fi)

    err_all_trials_mean1(k) = err_mean;
    err_all_trials_std1(k) = err_std;


% end 

% plot the power spectrum


%%
% save(fullfile(savePath,'all_trial_error_headDirection_CA1'),'err_all_trials_mean1','err_all_trials_std1') 

% rather than removing data based on speed, based on variance of the LFP data

% path tortuosity 
% the algo is kind of using angle
% linear kernel basis function! 
% does the head direction and moving direction coincide when the rat is
% traveling? 

% decoding peformance during aimless walking vs. goal-oriented (?) 
% tangent line from the point now to few points in the future 
% take the difference between current head direction and the tangent line
% goal-oriented head direction (?)  

% border as signed distance functions










