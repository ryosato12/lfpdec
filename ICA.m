clear;
close all; 

% -- session --
rats = {
    'tk0056-171014-04';
    % 'tk0056-171015-04';
    % 'tk0062-171129-04';
    % 'tk0062-171130-04';
    % 'tk0064-180206-04';
    % 'tk0064-180207-04';
    % 'tk0067-180426-04';
    % 'tk0067-180427-04';
    % 'tk0068-180530-04';
    % 'tk0068-180531-04';
    % 'tk0069-180726-04';
    % 'tk0069-180727-04';
    % 'tk0070-180829-04';
    % 'tk0070-180830-04';
    % 'tk0072-181025-04';
    % 'tk0072-181026-04';
    % 'tk0074-181219-04';
    % 'tk0074-181220-04';
    % 'tk0075-190326-04';
    % 'tk0075-190327-04';
    % 'tk0076-190424-04';
    % 'tk0076-190425-04'
    };

rats = strcat(rats,'rawData-7-18'); 
addpath('posbytrials');      % Change to your data foler
savePath = 'posbytrials'; 


%%

for k = 1:size(rats,1)

    session = rats{k};
    load([session '.mat'])    

    %% 
    % Use only data where spd>5% max
    vel = [0; diff(pos(:,1))];
    spd = smooth(sqrt(vel.^2));
    sidx = spd>(max(spd)*0.05);

    Xf = Xf(:,1:46875);
    [A,W] = ACMNsym(Xf(:,sidx),'mle_circ'); 
    Act = W*Xf; % W or conj of W 
    
    [~,score] = pca(Xf.');
    v = score(:,1); 
    demodAct = Act .* exp(-1i * angle(v.'));


    %% 
    bins = 100; 
    pos_binned = pos(:,1) - min(pos(:,1)) + eps;
    pos_binned = ceil(pos_binned/max(pos_binned)*bins);

    trial_num = w;


    %% 
    for i = 1:size(Act,1) % use unwrapped position 
        act_by_trial_and_pos{i} = accumarray([trial_num(sidx).' pos_binned(sidx)],demodAct(i,sidx).',[],@mean);
    end


    %%
    % for i = 1:size(ActSUB,1)
    %     im = complexIm(act_by_trial_and_posSUB{i});
    %     imshow(im);
    %     magnificationFactor = 5; 
    %     im = imresize(im,magnificationFactor,"nearest");
    %     imwrite(im,['posbytrials/',[session, '-SUB-Ch', num2str(i)],'.jpg']); 
    % end

    % for i = 1:size(ActSUB,1)
    %     im = complexIm(act_by_trial_and_pos{i});
    %     imshow(im);
    %     magnificationFactor = 5; 
    %     im = imresize(im,magnificationFactor,"nearest");
    %     imwrite(im,['posbytrials/',[session, '-ALL-Ch', num2str(i)],'.jpg']); 
    % end


    %% save figures
    f = figure(1);
    % f.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        % im = complexIm(postact_by_ch{i});
        % im = complexIm(velact_by_ch{i});
        % im = complexIm(accact_by_ch{i});
        % im = complexIm(hdact_by_ch{i});
        im = complexIm(act_by_trial_and_pos{i});
        imagesc(im)
        ylabel('Trial #')
        title(['Ch #',num2str(i)])
    end

    g = figure(2);
    % g.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        % im = complexIm(postact_by_ch{i+64});
        % im = complexIm(velact_by_ch{i+64});
        % im = complexIm(accact_by_ch{i+64});
        % im = complexIm(hdact_by_ch{i+64});
        im = complexIm(act_by_trial_and_pos{i+64});
        imagesc(im)
        ylabel('Trial #')
        title(['Ch #',num2str(i+64)])

    end

    h = figure(3);
    % h.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        % im = complexIm(postact_by_ch{i+128});
        % im = complexIm(velact_by_ch{i+128});
        % im = complexIm(accact_by_ch{i+128});
        % im = complexIm(hdact_by_ch{i+128});
        im = complexIm(act_by_trial_and_pos{i+128});
        imagesc(im)
        ylabel('Trial #')
        title(['Ch #',num2str(i+128)])
    end

    k = figure(4);
    % k.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        % im = complexIm(postact_by_ch{i+192});
        % im = complexIm(velact_by_ch{i+192});
        % im = complexIm(accact_by_ch{i+192});
        % im = complexIm(hdact_by_ch{i+192});
        im = complexIm(act_by_trial_and_pos{i+192});
        imagesc(im)
        ylabel('Trial #')
        title(['Ch #',num2str(i+192)])
    end

    a = figure(5);
    % f.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i)       
        im = complexIm(reshape(A(:,i),[],8));
        imagesc(im)
        ylabel('Trial #')
        title(['Ch #',num2str(i)])
    end
    
    b = figure(6);
    % g.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        im = complexIm(reshape(A(:,i+64),[],8));    
        imagesc(im)
        title(['Ch #',num2str(i+64)])
    
    end
    
    c = figure(7);
    % h.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        im = complexIm(reshape(A(:,i+128),[],8));    
        imagesc(im)
        title(['Ch #',num2str(i+128)])
    end
    
    d = figure(8);
    % k.Position = [100 100 4961 3508];
    for i = 1:64
        subplot(8,8,i) 
        im = complexIm(reshape(A(:,i+192),[],8));    
        imagesc(im)
        title(['Ch #',num2str(i+192)])
    end


    % print(figure(1),fullfile(savePath,[session '1-64-ICs']),'-dtiff')
    % print(figure(2),fullfile(savePath,[session '65-128-ICs']),'-dtiff')
    % print(figure(3),fullfile(savePath,[session '129-192-ICs']),'-dtiff')
    % print(figure(4),fullfile(savePath,[session '193-256-ICs']),'-dtiff')

    % print(figure(5),fullfile(savePath,[session '1-64-EMs']),'-dtiff')
    % print(figure(6),fullfile(savePath,[session '65-128-EMs']),'-dtiff')
    % print(figure(7),fullfile(savePath,[session '129-192-EMs']),'-dtiff')
    % print(figure(8),fullfile(savePath,[session '193-256-EMs']),'-dtiff')

    

end













