%-------------------------------
%ROBUST REGRESSION CORRELATION
%-------------------------------
%depends on nirs toolbox

%overall correlations - all possible channel pairs
dataprefix = 'DYAD'; %change this to whatever your data's prefix is
currdirpath = strcat('demo_data',filesep,'demo_dyads_preproc');
currdir = dir(strcat(currdirpath,filesep,dataprefix,'*'));
numdyads = length(currdir);
numchannels = 35; %change this to whatever number of channels you have in each person
samprate = 3.91; %change this to whatever sampling rate you had 
dyad_corrs = nan(numchannels,numchannels,numdyads);

fprintf('\n\t Calculating prewhitened, preweighted correlations for all channel pairs ...\n')
reverseStr = '';
Elapsedtime = tic;
for i=1:numdyads
    dyad=currdir(i).name;
    scan = dir(strcat(currdirpath,filesep,dyad,filesep,'*convo')); %rename wildcard to your scan name 
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj1_preprocessed.mat'),'-mat');
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj2_preprocessed.mat'),'-mat');
    sz = size(z_oxy1);
    if sz(1)>=940 %if scan lasted longer than 4 min (or however long you care about)
        subj1 = z_oxy1(1:940,:); %just looking at oxy for now - do deoxy as well if you want
        subj2 = z_oxy2(1:940,:);
    
    for k=1:35
        for j=1:35
              if ~isnan(subj1(1,k)) && ~isnan(subj2(1,j))
                data = cat(2,subj1(round(samprate*10):end-round(samprate*10),k),subj2(round(samprate*10):end-round(samprate*10),j)); %trim off first and last ten seconds
                corr = nirs.sFC.ar_corr(data);
                dyad_corrs(k,j,i) = corr(1,2);
              else
                dyad_corrs(k,j,i) = NaN;
              end 
        end
    end
    end
msg = sprintf('\n\t participant number %d/%d ...',i,numdyads);
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);



%corr null dist: pairing random subjects who weren't in same dyad together.
iterations = 10; %number of bootstrapped null samples to take; small number here for speedier demo but should be 100+
nulldist_corrs = nan(numchannels,numchannels,iterations);

randsubj1 = randi([1 numdyads],1,iterations);
randsubj2 = randi([1 numdyads],1,iterations);
evenodd1 = randi([0 1],1,iterations);
evenodd2 = randi([0 1],1,iterations);

fprintf('\n\t Bootstrapping Correlation Null Distribution ...\n')
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
            
        for k=1:35
            for j=1:35
              if ~isnan(subj1(1,k)) && ~isnan(subj2(1,j))
                data = cat(2,subj1(round(samprate*10):end-round(samprate*10),k),subj2(round(samprate*10):end-round(samprate*10),j)); %trim off first and last ten seconds
                corr = nirs.sFC.ar_corr(data);
                nulldist_corrs(k,j,a) = corr(1,2);
              else
                nulldist_corrs(k,j,a) = NaN;
              end 
            end
        end
    end         
    msg = sprintf('\n\t iteration %d/%d ...',a,iterations);
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);

nulldist_corrs = nulldist_corrs(~isnan(nulldist_corrs));%makes a single vector, ignores channel label
fisherz = @(r)(log(1+r)-log(1-r))/2;
nulldist_corrs_z = arrayfun(fisherz, nulldist_corrs);
dyad_corrs_z = arrayfun(fisherz, dyad_corrs);
dyad_corrs_z = reshape(dyad_corrs_z,[numchannels,numchannels,numdyads]);
slst=sort(nulldist_corrs_z);
slength=length(slst);
corrcutoff=slst(end-round(slength*5/100));

%network estimation based on similarity matrices and null dist made above
binarynetwork = zeros(numchannels,numchannels,numdyads);
for z=1:numdyads
binarynetwork(:,:,z) = dyad_corrs_z(:,:,z)>corrcutoff;
end

%Find p-values of all channel pairings compared to the null dist
channel_pvals = ones(35,35);
for x=1:35
    for y=1:35
        testval = nanmean(dyad_corrs_z(x,y,:));
        pval=slst>testval;
        channel_pvals(x,y)=sum(pval)/slength;
    end
end