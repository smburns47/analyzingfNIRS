%-------------------------
%PHASE SYNCHRONY THINGS
%-------------------------
%depends on nirs toolbox

%instantaneous phase synchrony for all possible channel pairs - code from: Pedersen et al., 2017
dataprefix = 'DYAD'; %change this to whatever your data's prefix is
currdirpath = strcat('demo_data',filesep,'demo_dyads_preproc');
currdir = dir(strcat(currdirpath,filesep,dataprefix,'*'));
numdyads = length(currdir);
numchannels = 35; %change this to whatever number of channels you have in each person
samprate = 3.91; %change this to whatever sampling rate you had 

%need a narrow frequency band for phase synchrony to work
dyads_phasesynch_fast = nan(numchannels*2,numchannels*2,numdyads); %slow-4 freq band
dyads_phasesynch_slow = nan(numchannels*2,numchannels*2,numdyads); %slow-5 freq band

num_good_channels = nan(1,numdyads);

for a=1:numdyads
    dyad=currdir(a).name;
    scan = dir(strcat(currdirpath,filesep,dyad,filesep,'*convo')); %rename wildcard to your scan name 
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj1_preprocessed.mat'),'-mat');
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj2_preprocessed.mat'),'-mat');
    subj1 = z_oxy1;
    subj2 = z_oxy2;
  
    goodinds1 = find(~isnan(subj1(1,:)));
    goodinds2 = find(~isnan(subj2(1,:)));
    goodindscheck1 = ~isnan(subj1(1,:));
    goodindscheck2 = ~isnan(subj2(1,:));
    goodindscheck = [goodindscheck1 goodindscheck2];
    subj1brief = subj1(:,goodinds1);
    subj2brief = subj2(:,goodinds2);
    
    subj1filtered = nan(size(subj1brief));
    subj2filtered = nan(size(subj2brief));
    for x=1:size(subj1brief,2)
        subj1filtered(:,x) = hmrBandpassFilt(subj1brief(:,x), 3.91, 0.03, 0.07);
    end
    for x=1:size(subj2brief,2)
        subj2filtered(:,x) = hmrBandpassFilt(subj2brief(:,x), 3.91, 0.03, 0.07);
    end
    subjboth = [subj1filtered subj2filtered];
    [subjbothinnov,f] = nirs.math.innovations(subjboth,20);
    num_good_channels(a)=size(subjboth,2);
    PSfast = angle(hilbert(subjbothinnov')); % instantaneous phases of channelxtime data
    PSfast_matrix = zeros(size(PSfast,1),size(PSfast,1),size(PSfast,2));
    for time_point = 1:size(PSfast, 2)
        %channel by channel phase synchrony at each time point - CxCxT
        PSfast_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSfast(:,time_point)', PSfast(:,time_point))));
    end
    for x=1:(size(subj1,2)*2)
        if ~goodindscheck(x)
            PSfast_matrix = [PSfast_matrix(:,1:(x-1),:) nan(size(PSfast_matrix,1),1,size(PSfast,2)) PSfast_matrix(:,x:end,:)];
        end
    end
    for x=1:(size(subj1,2)*2)
        if ~goodindscheck(x)
            PSfast_matrix = [PSfast_matrix(1:(x-1),:,:); nan(1,size(PSfast_matrix,2),size(PSfast,2)); PSfast_matrix(x:end,:,:)];
        end
    end
    dyads_phasesynch_fast(:,:,a) = nanmean(PSfast_matrix,3);

    subj1filtered = nan(size(subj1brief));
    subj2filtered = nan(size(subj2brief));
    for x=1:size(subj1brief,2)
        subj1filtered(:,x) = hmrBandpassFilt(subj1brief(:,x), 3.91, 0.01, 0.03);
    end
    for x=1:size(subj2brief,2)
        subj2filtered(:,x) = hmrBandpassFilt(subj2brief(:,x), 3.91, 0.01, 0.03);
    end
    subjboth = [subj1filtered subj2filtered];
    [subjbothinnov,f] = nirs.math.innovations(subjboth,20);
    num_good_channels(a)=size(subjboth,2);
    PSslow = angle(hilbert(subjbothinnov')); % instantaneous phases of channelxtime data
    PSslow_matrix = zeros(size(PSslow,1),size(PSslow,1),size(PSslow,2));
    for time_point = 1:size(PSslow, 2)
        %channel by channel phase synchrony at each time point - CxCxT
        PSslow_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSslow(:,time_point)', PSslow(:,time_point))));
    end
    for x=1:(size(subj1,2)*2)
        if ~goodindscheck(x)
            PSslow_matrix = [PSslow_matrix(:,1:(x-1),:) nan(size(PSslow_matrix,1),1,size(PSslow,2)) PSslow_matrix(:,x:end,:)];
        end
    end
    for x=1:(size(subj1,2)*2)
        if ~goodindscheck(x)
            PSslow_matrix = [PSslow_matrix(1:(x-1),:,:); nan(1,size(PSslow_matrix,2),size(PSslow,2)); PSslow_matrix(x:end,:,:)];
        end
    end
    
    dyads_phasesynch_slow(:,:,a) = nanmean(PSslow_matrix,3);
end




%phase synch null dist: pairing random subjects who weren't in same dyad together.
iterations = 10; %number of bootstrapped null samples to take; small number here for speedier demo but should be 100+
nulldist_fast = nan(70,70,iterations);
nulldist_slow = nan(70,70,iterations);

randsubj1 = randi([1 numdyads],1,iterations);
randsubj2 = randi([1 numdyads],1,iterations);
evenodd1 = randi([0 1],1,iterations);
evenodd2 = randi([0 1],1,iterations);

fprintf('\n\t Bootstrapping Phase Synchrony Null Distribution ...\n')
reverseStr = '';
Elapsedtime = tic;
for a=1:iterations
    if randsubj1(a)~=randsubj2(a)
        dyad1=currdir(randsubj1(a)).name;
        dyad2=currdir(randsubj2(a)).name;
        scan1 = dir(strcat(currdirpath,filesep,dyad1,filesep,'*convo')); %rename wildcard to your scan name 
        scan2 = dir(strcat(currdirpath,filesep,dyad2,filesep,'*convo')); %rename wildcard to your scan name 
        if evenodd1(a)
            load(strcat(currdirpath,filesep,dyad1,filesep,scan1(1).name,filesep,scan1(1).name,'_subj1_preprocessed.mat'),'-mat');
            subj1 = z_oxy1;
        else
            load(strcat(currdirpath,filesep,dyad1,filesep,scan1(1).name,filesep,scan1(1).name,'_subj2_preprocessed.mat'),'-mat');
            subj1 = z_oxy2;
        end
        if evenodd2(a)
            load(strcat(currdirpath,filesep,dyad2,filesep,scan2(1).name,filesep,scan2(1).name,'_subj1_preprocessed.mat'),'-mat');
            subj2 = z_oxy1;
        else
            load(strcat(currdirpath,filesep,dyad2,filesep,scan2(1).name,filesep,scan2(1).name,'_subj2_preprocessed.mat'),'-mat');
            subj2 = z_oxy2;
        end
        length1 = length(subj1);
        length2 = length(subj2);
        smallerlength = min(length1,length2);
        subj1 = subj1(1:smallerlength,:);
        subj2 = subj2(1:smallerlength,:);
    
        goodinds1 = find(~isnan(subj1(1,:)));
        goodinds2 = find(~isnan(subj2(1,:)));
        goodindscheck1 = ~isnan(subj1(1,:));
        goodindscheck2 = ~isnan(subj2(1,:));
        goodindscheck = [goodindscheck1 goodindscheck2];
        subj1brief = subj1(:,goodinds1);
        subj2brief = subj2(:,goodinds2);
    
        subj1fast = nan(size(subj1brief));
        subj2fast = nan(size(subj2brief));
        subj1slow = nan(size(subj1brief));
        subj2slow = nan(size(subj2brief));
    
        for x=1:size(subj1brief,2)
            subj1fast(:,x) = hmrBandpassFilt(subj1brief(:,x), 3.91, 0.03, 0.07);
            subj1slow(:,x) = hmrBandpassFilt(subj1brief(:,x), 3.91, 0.01, 0.027);
        end
        for x=1:size(subj2brief,2)
            subj2fast(:,x) = hmrBandpassFilt(subj2brief(:,x), 3.91, 0.03, 0.07);
            subj2slow(:,x) = hmrBandpassFilt(subj2brief(:,x), 3.91, 0.01, 0.027);
        end
    
        oxybothfast = [subj1fast subj2fast];
        oxybothslow = [subj1slow subj2slow];
    
        PSfast = angle(hilbert(oxybothfast')); % instantaneous phases of channelxtime data
        PSfast_matrix = zeros(size(PSfast,1),size(PSfast,1),size(PSfast,2));
        PSslow = angle(hilbert(oxybothslow')); 
        PSslow_matrix = zeros(size(PSslow,1),size(PSslow,1),size(PSslow,2));
    
        for time_point = 1:size(PSfast, 2)
        %channel by channel phase synchrony at each time point - CxCxT
            PSfast_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSfast(:,time_point)', PSfast(:,time_point))));
            PSslow_matrix(:,:,time_point) = 1-abs(sin(bsxfun(@minus, PSslow(:,time_point)', PSslow(:,time_point))));
        end
    
        for x=1:70
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(:,1:(x-1),:) nan(size(PSfast_matrix,1),1,size(PSfast,2)) PSfast_matrix(:,x:end,:)];
                PSslow_matrix = [PSslow_matrix(:,1:(x-1),:) nan(size(PSslow_matrix,1),1,size(PSslow,2)) PSslow_matrix(:,x:end,:)];
            end
        end
        for x=1:70
            if ~goodindscheck(x)
                PSfast_matrix = [PSfast_matrix(1:(x-1),:,:); nan(1,size(PSfast_matrix,2),size(PSfast,2)); PSfast_matrix(x:end,:,:)];
                PSslow_matrix = [PSslow_matrix(1:(x-1),:,:); nan(1,size(PSslow_matrix,2),size(PSslow,2)); PSslow_matrix(x:end,:,:)];
            end
        end
    
        nulldist_fast(:,:,a) = nanmean(PSfast_matrix,3);
        nulldist_slow(:,:,a) = nanmean(PSslow_matrix,3);
    else
        nulldist_fast(:,:,a) = nan(70,70);
        nulldist_slow(:,:,a) = nan(70,70);
    end
    msg = sprintf('\n\t iteration %d/%d ...',a,iterations);
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);

nulldist_fast_inter = nulldist_fast(1:35,36:end,:);
nulldist_slow_inter = nulldist_slow(1:35,36:end,:);
nulldist_fast_inter = nulldist_fast_inter(~isnan(nulldist_fast_inter));%makes a single vector of inter connections, ignores channel label
nulldist_slow_inter = nulldist_slow_inter(~isnan(nulldist_slow_inter));
fisherz = @(r)(log(1+r)-log(1-r))/2;
nulldist_fast_z = arrayfun(fisherz, nulldist_fast_inter);
nulldist_slow_z = arrayfun(fisherz, nulldist_slow_inter);
dyad_fast_z = arrayfun(fisherz, dyads_phasesynch_fast(1:35,36:end,:));
dyad_slow_z = arrayfun(fisherz, dyads_phasesynch_slow(1:35,36:end,:));
slst_fast = sort(nulldist_fast_z);
slst_slow = sort(nulldist_slow_z);
slength_fast = length(slst_fast);
slength_slow = length(slst_slow);
phasecutoff_fast = slst_fast(end-round(slength_fast*5/100));
phasecutoff_slow = slst_slow(end-round(slength_slow*5/100));

%network estimation based on similarity matrices and null dist made above
binarynetwork_fast = zeros(numchannels,numchannels,numdyads);
binarynetwork_slow = zeros(numchannels,numchannels,numdyads);
for z=1:numdyads
binarynetwork_fast(:,:,z) = dyad_fast_z(:,:,z)>phasecutoff_fast;
binarynetwork_slow(:,:,z) = dyad_slow_z(:,:,z)>phasecutoff_slow;
end

%Find p-values of all channel pairings compared to the null dist
pvals_fast = ones(numchannels,numchannels);
pvals_slow = ones(numchannels,numchannels);
for x=1:numchannels
    for y=1:numchannels
        fastval = nanmean(dyad_fast_z(x,y,:));
        slowval = nanmean(dyad_slow_z(x,y,:));
        fast_pval = slst_fast>fastval;
        slow_pval = slst_slow>slowval;
        pvals_fast(x,y)=sum(fast_pval)/slength_fast;
        pvals_slow(x,y)=sum(slow_pval)/slength_slow;
    end
end
