clear;
close all;

rats = {
    'tk0056-171014-04-1';
    'tk0056-171015-04-1';
    'tk0062-171129-04-1';
    'tk0062-171130-04-1';
    'tk0064-180206-04-1';
    'tk0064-180207-04-1';
    'tk0067-180426-04-1';
    'tk0067-180427-04-1';
    'tk0068-180530-04-1';
    'tk0068-180531-04-1';
    'tk0069-180726-04-1';
    'tk0069-180727-04-1';
    'tk0070-180829-04-1';
    'tk0070-180830-04-1';
    'tk0072-181025-04-1';
    'tk0072-181026-04-1';
    'tk0074-181219-04-1';
    'tk0074-181220-04-1';
    'tk0075-190326-04-1';
    'tk0075-190327-04-1';
    'tk0076-190424-04-1';
    'tk0076-190425-04-1'
    };

addpath('preprocess_001')
savePath = 'results_001';


%% 
err_all_trials_mean = cell(0);
err_all_trials_std = cell(0); 

for k = 1:size(rats,1)

    session = rats{k};
    load([session '.mat'])


    %% Thin channel by ROI 
    CA1_channel_num = sum(lfpCA1(:)==1);
    SUB_channel_num = sum(lfpSUB(:)==1);
    min_channel_num = min(CA1_channel_num,SUB_channel_num);
    if min_channel_num < SUB_channel_num
        num2remove = SUB_channel_num - min_channel_num; 
        inds = find(lfpSUB == 1); 
        SUB_inds = lfpSUB; 
        for i = 1:num2remove
           SUB_inds(inds(i)) = 0;
        end
        CA1_Xf = Xf(logical(lfpCA1),1:46875);
        SUB_Xf = Xf(logical(SUB_inds),1:46875);
    elseif min_channel_num < CA1_channel_num 
        num2remove = CA1_channel_num - min_channel_num; 
        CA1_inds = lfpCA1; 
        inds = find(lfpCA1 == 1); 
        for i = 0:num2remove-1
            CA1_inds(inds(end-i)) = 0; 
        end
        CA1_Xf = Xf(logical(CA1_inds),1:46875);
        SUB_Xf = Xf(logical(lfpSUB),1:46875);
    elseif CA1_channel_num == SUB_channel_num
        CA1_Xf = Xf(logical(lfpCA1),1:46875);
        SUB_Xf = Xf(logical(lfpSUB),1:46875);
    end 
    

    %% Demodulate by ROI 
    demodAct = Xf_demod(:,1:46875).'; 

    [~,score1] = pca(CA1_Xf.');
    v1 = score1(:,1); 
    demodCA1 = CA1_Xf .* exp(-1i * angle(v1.')); 
    demodCA1 = demodCA1.'; 

    [~,score2] = pca(SUB_Xf.');
    v2 = score2(:,1); 
    demodSUB = SUB_Xf .* exp(-1i * angle(v2.'));  
    demodSUB = demodSUB.'; 

    
    %% 
    vel = [0; diff(pos(:,1))];
    spd = smooth(sqrt(vel.^2));
    sidx = spd>(mean(spd));
    rpos = rpos(sidx,:);
    pos = pos(sidx,:); 

    % Signals to use for decoding
    yAll=cell(0);
    yAll{1}=zscore([real(demodCA1) imag(demodCA1)]);
    yAll{2}=zscore([real(demodSUB) imag(demodSUB)]); 
    yAll{3}=zscore([real(demodAct) imag(demodAct)]);
    mnames={'CA1(A_d)','SUB(A_d)','ALL(A_d)'};
    
    clearvars errAll errpAll errmAll err_nd_All err_nd_All pos_hat
    
    for model_num = 1:length(yAll)
        y = yAll{model_num}(sidx,:);

        % Generate basis...
        [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*pi,400);
        pvec = linspace(0,2*pi,256);
        [tmp,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);
      
        % Cross-validated predictions
        m.cv = getCVidx(size(y,1),10);
        f=[];
        
        for i=1:m.cv.nfoldcv
            % Fit
            g = y(m.cv.tr{i},:)\Xrbf(m.cv.tr{i},:);
            f(m.cv.ts{i},:) = y(m.cv.ts{i},:)*(g*dbasis');
        end
        
        % Plotting and error calcs...
        % xl=[4000 8000]; % range to show predictions
        f = zscore(f')';
        idx=1:size(f,1); 
        [tmp,maxpost]=max(f'); 
        
        % decoderSummary
        trackLen = 235; 

        % Circular error...
        err = circ_dist(rpos(idx)*pi,maxpost'/length(pvec)*2*pi)*trackLen/pi;
        errp = circ_dist(rpos(idx)*pi,maxpost'/length(pvec)*2*pi)*trackLen/pi;
        errm = circ_dist(2*pi-rpos(idx)*pi,maxpost'/length(pvec)*2*pi)*trackLen/pi;
        
        errAll(:,model_num) = err;
        errpAll(:,model_num) = errp;
        errmAll(:,model_num) = errm;
        % non-directional error
        [err_nd_All(:,model_num),mini] = (min([abs(errmAll(:,model_num)) abs(errpAll(:,model_num))]'));
        pos_hat(:,model_num) = maxpost'/length(pvec)*2*pi;
    end


    %% Median error in cm with s.e.
    err_mean = [];
    err_std = [];
    
    for model_num=1:length(yAll)
        bootstat = bootstrp(500,'median',abs(err_nd_All(:,model_num)));
        err_mean(model_num) = mean(bootstat);
        err_std(model_num) = std(bootstat);
    end
  

    %% Index prediction error 
    err_all_trials_mean{k} = err_mean;
    err_all_trials_std{k} = err_std;


    %% Plot results 
    f = figure(1); clf
    f.WindowState = 'maximized';
    edges = linspace(0,2*pi,256);
    e{1}=edges;
    e{2}=edges;
    ex = edges+mean(diff(edges))/2;
    [mx,my] = meshgrid(ex,ex);

    pvec =[1 2 3]; 
    for i=1:length(pvec)
        subplot(1,3,i)
        [n,c]=hist3([rpos*pi pos_hat(:,pvec(i))],'Edges',e);
        nnorm = bsxfun(@rdivide,n,sum(n));
        idx = nnorm>10e-3;
        scatter(mx(idx)/pi*trackLen-trackLen,my(idx)/pi*trackLen-trackLen,nnorm(idx)*100,'filled','sk')
        axis image
        xlim([-1 1]*trackLen)
        ylim([-1 1]*trackLen)
        set(gca,'TickDir','out')
        set(gca,'XTick',[-trackLen 0 trackLen])
        set(gca,'YTick',[-trackLen 0 trackLen])
        title(mnames{pvec(i)})
    end

    % print(figure(1),fullfile(savePath,[session '-dec']),'-djpeg')


    %% Save results 
    % save(fullfile(savePath,[session '-removed']),'rpos','err_mean','err_std','err_nd_All','pos_hat','n','nnorm')


end 


%%
x = 1:20;
for i = x
    err_mean_CA1(i) = err_all_trials_mean{i}(1);
    err_mean_SUB(i) = err_all_trials_mean{i}(2);
    err_mean_ALL(i) = err_all_trials_mean{i}(3);
    err_std_CA1(i) = err_all_trials_std{i}(1);
    err_std_SUB(i) = err_all_trials_std{i}(2);
    err_std_ALL(i) = err_all_trials_std{i}(3);
end


%%
figure(2); clf
boxplot([err_mean_CA1' err_mean_SUB' err_mean_ALL'],{'CA1@LFP','SUB@LFP','ALL@LFP'})
ylabel('Predicition Error [cm]')
ylim([0 100])
grid on
hold on
plot(mean([err_mean_CA1' err_mean_SUB' err_mean_ALL']), 'dk')
hold off 
title('Prediction Error for CA1, SUB, and ALL')
dim = [0.6 0.5 0.3 0.3];
a = annotation('textbox',dim,'String','p-value: 0.2358','FitBoxToText','on');
a.Color = 'red';
a.FontSize = 10;
a.HorizontalAlignment = 'center';

% not signficant, p-value: 0.2358 
% [tbl,chi,pval] = crosstab(err_mean_CA1,err_mean_SUB);

% print(figure(2),fullfile(savePath,[session '-box-removed']),'-djpeg')


%%
boxplot([err_mean_CA1; err_mean_SUB])


%%
figure(3); clf;
b = bar(x,err_mean_CA1,'k'); 
% b.CData(14,:) = [1 0 0];
hold on 
er = errorbar(x,err_mean_CA1,err_std_CA1,err_std_CA1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth = 1;
hold off
ylim([0 100])
xlabel('Rat')
ylabel('Median Error [cm]')
title('CA1 prediction error')

% print(figure(3),fullfile(savePath,[session '-CA1-bar-removed']),'-djpeg')


%%
figure(4); clf;
b = bar(x,err_mean_SUB,'k'); 
% b.CData(14,:) = [1 0 0];
hold on 
er = errorbar(x,err_mean_SUB,err_std_SUB,err_std_SUB);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth = 1;
hold off
ylim([0 100])
xlabel('Rat')
ylabel('Median Error [cm]')
title('SUB prediction error')

% print(figure(4),fullfile(savePath,[session '-SUB-bar-removed']),'-djpeg')


%%
figure(5); clf;
b = bar(x,err_mean_ALL,'k'); 
% b.CData(14,:) = [1 0 0];
hold on 
er = errorbar(x,err_mean_ALL,err_std_ALL,err_std_ALL);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth = 1;
hold off
ylim([0 100])
xlabel('Rat')
ylabel('Median Error [cm]')
title('ALL prediction error')

% print(figure(5),fullfile(savePath,[session '-bar-removed']),'-djpeg')


%% 
% what percentage is in the right position given that the rat moves in the
% correct direction (right/left)


% save(fullfile(savePath,'summary-stats-removed'),'err_mean_ALL','err_mean_CA1','err_mean_SUB')



