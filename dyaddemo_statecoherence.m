%------------------------
% STATE COHERENCE
%------------------------

dataprefix = 'DYAD'; %change this to whatever your data's prefix is
currdirpath = strcat('demo_data',filesep,'demo_dyads_preproc');
currdir = dir(strcat(currdirpath,filesep,dataprefix,'*'));
numdyads = length(currdir);
numchannels = 35; %change this to whatever number of channels you have in each person
samprate = 3.91; %change to your sampling rate
offdiag = round(60*samprate); %consider every off-diagonal up to 60sec of time difference

cohesionlist = nan(numdyads,offdiag); %consider every off-diagonal up to 1-minute of time difference
dyadlengthlist = nan(numdyads,1);
numchannelslist = nan(numdyads,1);

fprintf('\n\t Calculating Experimental State Cohesion ...\n')
reverseStr = '';
Elapsedtime = tic;
for a=1:numdyads
    dyad=currdir(a).name;
    scan = dir(strcat(currdirpath,filesep,dyad,filesep,'*convo')); %rename wildcard to your scan name 
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj1_preprocessed.mat'),'-mat');
    load(strcat(currdirpath,filesep,dyad,filesep,scan(1).name,filesep,scan(1).name,'_subj2_preprocessed.mat'),'-mat');
    %comment in if wanting just first/last 4 minutes
    %add just channel range of dmn or executive network 
    %oxy1 = oxy1(1:940,:);
    %oxy2 = oxy2(1:940,:);
    %oxy1 = oxy1(end-940:end,:);
    %oxy2 = oxy2(end-940:end,:);
    subj1=z_oxy1;
    subj2=z_oxy2;
   
    goodinds1 = find(~isnan(subj1(1,:)));
    goodinds2 = find(~isnan(subj2(1,:)));
    goodindscheck1 = ~isnan(subj1(1,:));
    goodindscheck2 = ~isnan(subj2(1,:));
    goodindscheck = goodindscheck1 & goodindscheck2;
    subj1brief = subj1(:,goodindscheck);
    subj2brief = subj2(:,goodindscheck);
    
    dyadlengthlist(a) = length(subj1brief);
    numchannelslist(a) = sum(goodindscheck);
    
    simdot1 = zeros(length(subj1brief), length(subj1brief));
    simdot2 = zeros(length(subj2brief), length(subj2brief));
    for i = 1:length(subj1brief)
        data = [subj1brief(i,:); subj1brief];
        ms = mean(data,2);
        datam = data - repmat(ms, 1, size(data,2));
        datass = sqrt(sum(datam'.^2)');
        temp = dot(datam(2:end,:)', repmat(datam(1,:), size(datam(2:end,:),1),1)');
        rs = temp' ./ (datass(2:end) * datass(1));
        simdot1(:,i) = rs;
        
        data = [subj2brief(i,:); subj2brief];
        ms = mean(data,2);
        datam = data - repmat(ms, 1, size(data,2));
        datass = sqrt(sum(datam'.^2)');
        temp = dot(datam(2:end,:)', repmat(datam(1,:), size(datam(2:end,:),1),1)');
        rs = temp' ./ (datass(2:end) * datass(1));
        simdot2(:,i) = rs;
    end
    
    for t=1:offdiag
        subj1=nan(length(subj1brief)-t,1);
        subj2=nan(length(subj1brief)-t,1);
        for j=1:length(subj1)
            subj1(j,1)=simdot1(j,j+t);
            subj2(j,1)=simdot2(j,j+t);
        end
        cohesion = corrcoef(subj1,subj2);
        cohesionlist(a,t) = cohesion(1,2);
        
    end
    msg = sprintf('\n\t dyad number %d/%d ...',a,numdyads);
    fprintf([reverseStr,msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);


%bootstrap the null state dependency distribution (pairing random,
%non-iteracting partners)
iterations = 50; %number of bootstrapped null samples to take; small number here for speedier demo but should be 100+
null_cohesionlist = nan(iterations,offdiag);
null_dyadlengthlist = nan(iterations);
null_numchannelslist = nan(iterations);
randsubj1 = randi([1 numdyads],1,iterations);
randsubj2 = randi([1 numdyads],1,iterations);
evenodd1 = randi([0 1],1,iterations);
evenodd2 = randi([0 1],1,iterations);

fprintf('\n\t Bootstrapping Experimental State Dependency Null Distribution ...\n')
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
        goodinds1 = find(~isnan(subj1(1,:)));
        goodinds2 = find(~isnan(subj2(1,:)));
        goodindscheck1 = ~isnan(subj1(1,:));
        goodindscheck2 = ~isnan(subj2(1,:));
        goodindscheck = goodindscheck1 & goodindscheck2;
        subj1brief = subj1(1:smallerlength,goodindscheck);
        subj2brief = subj2(1:smallerlength,goodindscheck);
            
        %comment in if wanting just first/last 4 minutes
        %add just channel range of dmn or executive network
        %only works if all dyads in data are >4min
        %subj1brief = subj1brief(1:940,:);
        %subj2brief = subj2brief(1:940,:);
        %subj1brief = subj1brief(end-940:end,:);
        %subj2brief = subj2brief(end-940:end,:);
    
        null_dyadlengthlist(a) = length(subj1brief);
        null_numchannelslist(a) = sum(goodindscheck);
    
        simdot1 = zeros(length(subj1brief), length(subj1brief));
        simdot2 = zeros(length(subj2brief), length(subj2brief));
        for i = 1:length(subj1brief)
            data = [subj1brief(i,:); subj1brief];
            ms = mean(data,2);
            datam = data - repmat(ms, 1, size(data,2));
            datass = sqrt(sum(datam'.^2)');
            temp = dot(datam(2:end,:)', repmat(datam(1,:), size(datam(2:end,:),1),1)');
            rs = temp' ./ (datass(2:end) * datass(1));
            simdot1(:,i) = rs;
        
            data = [subj2brief(i,:); subj2brief];
            ms = mean(data,2);
            datam = data - repmat(ms, 1, size(data,2));
            datass = sqrt(sum(datam'.^2)');
            temp = dot(datam(2:end,:)', repmat(datam(1,:), size(datam(2:end,:),1),1)');
            rs = temp' ./ (datass(2:end) * datass(1));
            simdot2(:,i) = rs;
        end
    
        for t=1:offdiag
            subj1=nan(length(subj1brief)-t,1);
            subj2=nan(length(subj1brief)-t,1);
            for j=1:length(subj1)
                subj1(j,1)=simdot1(j,j+t);
                subj2(j,1)=simdot2(j,j+t);
            end
            cohesion = corrcoef(subj1,subj2);
            null_cohesionlist(a,t) = cohesion(1,2);
        end
    end
msg = sprintf('\n\t dyad/iteration %d/%d ...',a,iterations);
fprintf([reverseStr,msg]);
reverseStr = repmat(sprintf('\b'),1,length(msg));
end
Elapsedtime = toc(Elapsedtime);
fprintf('\n\t Elapsed time: %g seconds ...\n', Elapsedtime);

cohesioncutoff = nan(1,t);
t_pvals = ones(1,t);
fisherz = @(r)(log(1+r)-log(1-r))/2;
for t=1:offdiag
    nullt = null_cohesionlist(~isnan(null_cohesionlist(:,t)),t);
    nullt_z = arrayfun(fisherz, nullt);
    cohesiont = arrayfun(fisherz, cohesionlist(:,t)); 
    slst=sort(nullt);
    slength=length(slst);
    cohesioncutoff(1,t)=slst(end-round(slength*5/100));
    testval = nanmean(cohesiont);
    pval = slst>testval;
    t_pvals(1,t) = sum(pval)/slength;
end