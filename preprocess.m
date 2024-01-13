clear;
close all; 

% addpath('anton')
% addpath('toolbox')

% -- session, mapfile, trange --
A = {
% 'tk0056-171014-04', 'Buzsaki256.mat',[0 1200]+20;
% 'tk0056-171015-04', 'Buzsaki256.mat',[0 1200]+15;
% 'tk0062-171129-04', 'Buzsaki256.mat',[0 1200]+10;
% 'tk0062-171130-04', 'Buzsaki256.mat',[0 1200]+30;
% 'tk0064-180206-04','Buzsaki256.mat',[0 1200]+15;
% 'tk0064-180207-04','Buzsaki256.mat',[0 1200]+15;
% 'tk0067-180426-04', 'Buzsaki256.mat',[0 1200]+10;
% 'tk0067-180427-04','Buzsaki256.mat',[0 1200]+20;
% 'tk0068-180530-04','A8x32Edge.mat', [0 1200]+20;
% 'tk0068-180531-04','A8x32Edge.mat', [0 1200]+15;
% 'tk0069-180726-04','A8x32_5mm_35_300_160.mat',[0 1200]+10;
% 'tk0069-180727-04','A8x32_5mm_35_300_160.mat',[0 1200]+40;
'tk0070-180829-04','A8x32_5mm_35_300_160.mat',[0 1200]+12; % best performing rat 
% 'tk0070-180830-04','A8x32_5mm_35_300_160.mat',[0 1200]+10;
% 'tk0072-181025-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0072-181026-04','A8x32Edge.mat',[0 1200]+12;
% 'tk0074-181219-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0074-181220-04','A8x32Edge.mat',[0 1200]+20;
% 'tk0075-190326-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0075-190327-04','A8x32Edge.mat',[0 1200]+15;
% 'tk0076-190424-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0076-190425-04','A8x32Edge.mat',[0 1200]+15
};

basePath = 'data';      % Change to your data foler
savePath = 'preprocess_001'; 


%%
for k = 1:size(A,1)
    
    session = A{k,1};  % linear track data
    mapfile = A{k,2};  % channelmap ('Buzsaki256.mat', 'A8x32Edge.mat', or 'A8x32_5mm_35_300_160.mat')
    trange = A{k,3};   % time range (sec)

    
    %% Index session 
    idx = strfind(session,'-');
    rat = session(1:idx(1)-1);
    day = session(idx(1)+1:idx(2)-1);
    sess= session(idx(2)+1:end);
    
    % sampling rate of camera
    fsv = 39.0625; % hz
    dt = 1/fsv;


    %% Partition ROI
    % T = readtable('histology_test.xlsx','Sheet',rat,'Range','A1:H32');
    % 
    % TT = T.Variables;
    % 
    % lfpSUB = TT == 1;
    % lfpSUBmol = TT == 2; 
    % lfpCA1 = TT == 3;
    % lfpCA3 = TT == 5; 
    % 
    % % indices for ROI
    % lfpSUB = reshape(fliplr(lfpSUB),[],1); 
    % lfpSUBmol = reshape(fliplr(lfpSUBmol),[],1);
    % lfpCA1 = reshape(fliplr(lfpCA1),[],1); 
    % lfpCA3 = reshape(fliplr(lfpCA3),[],1);


    %% Linearize rpos
    % whl file path
    whlPath = fullfile(basePath,rat,[day '-' sess]);
    
    % position scaling
    % 0.5217 (230pixels = 120cm) for open field, T-maze, Zigzag maze
    % 0.4130 (569pixels = 235cm) for linear track
    pixel2cm = 0.4130;
    
    % load pos (cm)
    % position is sampled at 39.0625 Hz (25.6 ms/sample)
    [t,x1,y1,x2,y2] = loadWhl(whlPath,pixel2cm);
    
    % flip y data
    y1 = -y1;
    y2 = -y2;

    % trim position data
    if trange(1)<t(1) || t(end)<trange(2)
        error 'Error: trange setting is out of range!'
    end
    idx = trange(1)<=t & t<trange(2);
    x1(~idx)=[];
    y1(~idx)=[];
    x2(~idx)=[];
    y2(~idx)=[];
    t(~idx)=[];
    
    % use the better-tracking LED
    if sum(x1==-1)>sum(x2==-1)
        posx = x2;
        posy = y2;
    else
        posx = x1;
        posy = y1;
    end
    post = t;

    % input for fix_position 
    pos = [x1 y1 x2 y2];
    [~,rpos,~,trial_num] = fix_position(pos); 


    %% Spike counts (39.0625 Hz = 25.6 ms time bin)
    % spike file path
    spkPath = fullfile(basePath,rat,[day '_sorted']);

    % load spike data 
    [ts,clu]=loadSpkHist(spkPath,[],session);

    % trim time range
    idx = trange(1)<ts & ts<trange(2);
    ts  = ts(idx)-trange(1);
    clu = clu(idx);

    % spike counts
    CLU = clu==1:max(clu);  % <-cell selection
    TS  = ts(any(CLU,2));
    CLU = CLU(any(CLU,2),:);

    spk = zeros(length(pos(:,1)),max(clu));
    for i=1:size(spk,1)
        idx = post(i)-dt/2<=TS & TS<post(i)+dt/2;
        spk(i,:) = sum(CLU(idx,:),1);
    end


    %% Local-field potential (1250-Hz)
    % lfp file path
    lfpPath = fullfile(basePath,rat,[day '-' sess]);

    % load lfp
    [lfp,t,~] = loadLfp(lfpPath,[],mapfile);
    lfp = lfp(1:256,:); 
    idx = trange(1)<=t & t<=trange(2);
    lfp = lfp(:,idx);
    t   = t(idx)-trange(1);


    %% Decimate by a factor of 32 (39.0625 Hz) 
    indlfp = 32;
    X = []; 
    for ii = 1:size(lfp,1)
        X(ii,:) = decimate(lfp(ii,:),indlfp,'fir');
    end


    %% Plot raw and decimated signal 
    lx = length(lfp); 
    figure(1); clf
    plot(0:lx-1,lfp(1,:))
    hold on 
    % stem(0:indlfp:lx-1,X(1,:),'ro','filled','MarkerSize',4)
    plot(0:indlfp:lx-1,X(1,:),'-r')
    hold off
    xlim([2000 2500])
    legend('Raw','Decimated')
    xlabel('Sample Number')
    ylabel('Signal [mV]')
    title('Raw and decimated signal [Ch #1]')

    print(figure(1),fullfile(savePath,[session '-3']),'-djpeg')


    %% Demodulation 
    % Hilbert transform 
    Fc = 8.4;             % center freq (Hz)
    Fs = 1250/indlfp;     % sampling rate (Hz)
    Xf = morFilter(X,Fc,Fs);

    % Demodulate signal 
    [coeff,score,latent,tsquared,explained] = pca(Xf.');
    v = score(:,1);
    Xf_demod = Xf .* exp(-1i * angle(v.'));


    %% Summmary plot 
    % plot
    figure(2); clf   
    idx = (2000:2500);
    taxis = idx/fsv;
    pad = repmat(-(1:size(X(:,idx),1))',1,size(X(:,idx),2))*100;

    subplot(3,1,1)
    imagesc(taxis,[],X(:,idx))
    title('raw LFP')
    ylabel('Channel')

    subplot(3,1,2)
    imagesc(taxis,[],abs(Xf_demod(:,idx)))
    colormap hsv
    ylabel('Channel')
    title('Amplitude of demodulated LFP')

    subplot(3,1,3)
    imagesc(taxis,[],angle(Xf_demod(:,idx)))
    colormap hsv
    xlabel('Time (s)')
    ylabel('Channel')
    title('Phase of demodulated LFP')
    drawnow

    print(figure(2),fullfile(savePath,[session '-4']),'-djpeg')


    %% Save results 
    % save(fullfile(savePath,[session '-1']),'pos','post','rpos','trial_num','TS','CLU','spk','t','X','Xf','Xf_demod', ...
    % 'coeff','score','latent','tsquared','explained','lfpCA1','lfpCA3','lfpSUB','lfpSUBmol') 
    

    %%
    save(fullfile(savePath,[session '-1']),'trial_num','-append') 

   
end




%%
% function roiName = getRoiName(roi)
% 
% switch roi
%     case -1; roiName = 'ND';
%     case 0; roiName = 'Other';
%     case 1; roiName = 'SUB';
%     case 2; roiName = 'SUBmol';
%     case 3; roiName = 'CA1';
%     case 4; roiName = 'DG';
%     case 5; roiName = 'CA3';
%     otherwise; roiName = 'ND';
% end


