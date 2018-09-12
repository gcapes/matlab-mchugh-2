function [fittedParams,fitRes,sToUse,...
    paramsToUse,weightsThisVoxToUse]=fitMM_voxelwise(ims,imsNorm,...
    scanParams,initialParamsMat,voxelwiseInitialisation,...
    riceNoiseStd,lb,ub,weights,diffModel,ls_mle_str)
%% Function for fitting microstructural model on a voxel-wise basis
% Inputs: ims - image signal intensities
%         imsNorm - normalised signal intensities (normalised to G=0)
%         scanParams - matrix of scan parameters
%         initalParams - matrix of initial parameters for fitting
%         riceNoiseStd - Rician noise values
%         lb,ub - lower and upper fit bounds
%         diffModel - string specifying model to fit

%% Display model and fit type
disp(strcat('Fitting',{' '},diffModel,{' '},'using',{' '},ls_mle_str))

%% Noise - Gudbjartsson H, Patz S. Magn Reson Med 1995; 34:910–914.
noiseMean=riceNoiseStd.*sqrt(pi/2);

%% Plot or not?
plotYesNo='n'; %'n'/'y'

%% Check that number of signals matches number of scan parameters
if size(imsNorm,3)~=size(scanParams,1)
	error(['!!! Error - number of signals and number of scan '...
	        'parameters do not match']);
end

%% Option for fitting ls or mle
switch ls_mle_str
case 'mle'
	ls_mle='mle'
    disp('*********************************')
    disp('******TESTING MLE FITTING******')
    disp('*********************************')
otherwise 
    ls_mle='ls';
end

%% Check dimensions of initialisation maps if being used 
switch voxelwiseInitialisation
    case 'voxelInit'
        if all([size(ims(:,:,1)) numel(lb)]==size(initialParamsMat))~=1
            error('!!! error with dimensions of voxelwise initialisation maps')
        end
    case 'standardInit'
        %nothing to do
end
   
%% Initialise output maps
switch diffModel
    case {'DiDe','DiDe_with_rise_time'}
        rMap=zeros(size(ims,1),size(ims,2));
        diMap=zeros(size(ims,1),size(ims,2));
        deMap=zeros(size(ims,1),size(ims,2));
        fMap=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
    case 'DiDe_unnormalised'
        rMap=zeros(size(ims,1),size(ims,2));
        diMap=zeros(size(ims,1),size(ims,2));
        deMap=zeros(size(ims,1),size(ims,2));
        fMap=zeros(size(ims,1),size(ims,2));
        s0Map=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
    case 'DiDe_notort'
        rMap=zeros(size(ims,1),size(ims,2));
        diMap=zeros(size(ims,1),size(ims,2));
        deMap=zeros(size(ims,1),size(ims,2));
        fMap=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));    
    case 'DiDe_fiso'
        rMap=zeros(size(ims,1),size(ims,2));
        diMap=zeros(size(ims,1),size(ims,2));
        deMap=zeros(size(ims,1),size(ims,2));
        fiMap=zeros(size(ims,1),size(ims,2));
        feMap=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        extFlgMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
    case 'DiDe_fit_fiso'
        rMap=zeros(size(ims,1),size(ims,2));
        diMap=zeros(size(ims,1),size(ims,2));
        deMap=zeros(size(ims,1),size(ims,2));
        fiMap=zeros(size(ims,1),size(ims,2));
        feMap=zeros(size(ims,1),size(ims,2));
        fisoMap=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        extFlgMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
    case 'D_fiso'
        rMap=zeros(size(ims,1),size(ims,2));
        dMap=zeros(size(ims,1),size(ims,2));
        fiMap=zeros(size(ims,1),size(ims,2));
        feMap=zeros(size(ims,1),size(ims,2));
        s0Map=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        extFlgMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));  
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
    case 'D'
        rMap=zeros(size(ims,1),size(ims,2));
        dMap=zeros(size(ims,1),size(ims,2));
        fMap=zeros(size(ims,1),size(ims,2));
        s0Map=zeros(size(ims,1),size(ims,2));
        rsqMap=zeros(size(ims,1),size(ims,2));
        strtValMap=zeros(size(ims,1),size(ims,2));
        numSig=zeros(size(ims,1),size(ims,2));
        fitRes=zeros(size(ims,1),size(ims,2));
        sToUse=zeros(size(ims,1),size(ims,2));
        paramsToUse=zeros(size(ims,1),size(ims,2));
        weightsThisVoxToUse=zeros(size(ims,1),size(ims,2));
	otherwise
		error('!!! Unrecognised diffModel')
end

switch voxelwiseInitialisation
case 'standardInit'
	initialParams=initialParamsMat;
end

%% Get signals to fit, and pass to fitMM_ls_mle
for i=1:size(ims,1)
    for j=1:size(ims,2)
        s=squeeze(ims(i,j,:));
        % If any signals are 0, don't perform fitting
        if all(s(:))==1
        sNorm=squeeze(imsNorm(i,j,:));       
        weightsThisVox=squeeze(weights(i,j,:));
        % Perform fitting if we're not in the noise
        % 3xricianNoiseStd mostly excludes the noise
        if s(1)>3*riceNoiseStd(1)%2*riceNoiseStd/sqrt(2/pi);
            switch voxelwiseInitialisation
                case 'voxelInit'
                    % Pass voxel-specific initialisation values
                    initialParams=squeeze(...
                        initialParamsMat(i,j,:));
                case 'standardInit'
                    %nothing to do
            end
            % Exclude signals where intensity is lower than 2*noise ROI
            % mean (Portnoy et al. MRM 2013)
            % get corresponding normalised signals
            switch ls_mle
                case 'ls'
                    if numel(noiseMean)>1
                        error('!!! LS fitting and >1 noiseMean value !')
                    end
                    sToUse=sNorm(s>2*noiseMean(1));
                    paramsToUse=scanParams(s>2*noiseMean(1),:);
                    weightsThisVoxToUse=weightsThisVox(s>2*noiseMean(1));
                case 'mle'
                    sToUse=sNorm(s>2*0);
                    paramsToUse=scanParams(s>2*0,:);
                    weightsThisVoxToUse=weightsThisVox(s>2*0);
                    switch diffModel
                        case 'DiDe_unnormalised'
                            riceNoiseStdNorm=riceNoiseStd;
                        otherwise
                            % Normalise rician noise (see signalNormalisation.m for tests)
                            % Use b0 signals as scaling
                            riceNoiseStdNorm=zeros(numel(sToUse),1);
                            [ind,~]=find(sNorm==1);
                            %riceNoiseStdNorm=riceNoiseStd./s(ind(1));
                            %{%
                            riceNoiseStdNorm(ind)=(riceNoiseStd./s(ind));
                            id = find(riceNoiseStdNorm);
                            riceNoiseStdNorm(id(2:end)) = ...
                                diff(riceNoiseStdNorm(id));
                            riceNoiseStdNorm = ...
                                cumsum(riceNoiseStdNorm);
                            riceNoiseStdNorm = ...
                                riceNoiseStdNorm.*(sqrt(1+sToUse.^2));
                            %}
                    end
            end
            
            switch diffModel
                case {'DiDe','DiDe_with_rise_time'}
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                diffModel,'ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                diffModel,'mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
                    % Collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    diMap(i,j)=paramEstimatesLS(2);
                    deMap(i,j)=paramEstimatesLS(3);
                    fMap(i,j)=paramEstimatesLS(4);
                    rsqMap(i,j)=paramEstimatesLS(5);
                    strtValMap(i,j)=paramEstimatesLS(6);
                    numSig(i,j)=numel(sToUse);
                case 'DiDe_unnormalised'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_unnormalised','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_unnormalised','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
                    % Collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    diMap(i,j)=paramEstimatesLS(2);
                    deMap(i,j)=paramEstimatesLS(3);
                    fMap(i,j)=paramEstimatesLS(4);
                    s0Map(i,j)=paramEstimatesLS(5);
                    rsqMap(i,j)=paramEstimatesLS(6);
                    strtValMap(i,j)=paramEstimatesLS(7);
                    numSig(i,j)=numel(sToUse);
                case 'DiDe_notort'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_notort','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_notort','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
                    % Collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    diMap(i,j)=paramEstimatesLS(2);
                    deMap(i,j)=paramEstimatesLS(3);
                    fMap(i,j)=paramEstimatesLS(4);
                    rsqMap(i,j)=paramEstimatesLS(5);
                    strtValMap(i,j)=paramEstimatesLS(6);
                    numSig(i,j)=numel(sToUse);
                case 'DiDe_fiso'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_fiso','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_fiso','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
                    % collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    diMap(i,j)=paramEstimatesLS(2);
                    deMap(i,j)=paramEstimatesLS(3);
                    fiMap(i,j)=paramEstimatesLS(4);
                    feMap(i,j)=paramEstimatesLS(5);
                    rsqMap(i,j)=paramEstimatesLS(6);
                    strtValMap(i,j)=paramEstimatesLS(7);
                    extFlgMap(i,j)=paramEstimatesLS(8);
                    numSig(i,j)=numel(sToUse);
                case 'DiDe_fit_fiso'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_fit_fiso','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'DiDe_fit_fiso','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
					% collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    diMap(i,j)=paramEstimatesLS(2);
                    deMap(i,j)=paramEstimatesLS(3);
                    fiMap(i,j)=paramEstimatesLS(4);
                    feMap(i,j)=paramEstimatesLS(5);
                    fisoMap(i,j)=paramEstimatesLS(6);
                    rsqMap(i,j)=paramEstimatesLS(7);
                    strtValMap(i,j)=paramEstimatesLS(8);
                    extFlgMap(i,j)=paramEstimatesLS(9);
                    numSig(i,j)=numel(sToUse);
                case 'D_fiso'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'D_fiso','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'D_fiso','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
					% collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    dMap(i,j)=paramEstimatesLS(2);
                    fiMap(i,j)=paramEstimatesLS(3);
                    feMap(i,j)=paramEstimatesLS(4);
                    s0Map(i,j)=paramEstimatesLS(5);
                    rsqMap(i,j)=paramEstimatesLS(6);
                    strtValMap(i,j)=paramEstimatesLS(7);
                    extFlgMap(i,j)=paramEstimatesLS(8);
                    numSig(i,j)=numel(sToUse);
                case 'D'
                    switch ls_mle
                        case 'ls'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'D','ls',[]);
                        case 'mle'
                            [paramEstimatesLS, fitRes]=mm_fitting.fitMM_ls_mle(...
                                sToUse,paramsToUse,initialParams,lb,ub,plotYesNo,...
                                'D','mle',riceNoiseStdNorm);
                        otherwise
                            error('!!! ls_mle must be ls or mle')
                    end
                    
                    % collect each parameter for map
                    rMap(i,j)=paramEstimatesLS(1);
                    dMap(i,j)=paramEstimatesLS(2);
                    fMap(i,j)=paramEstimatesLS(3);
                    s0Map(i,j)=paramEstimatesLS(4);
                    rsqMap(i,j)=paramEstimatesLS(5);
                    strtValMap(i,j)=paramEstimatesLS(6);
                    numSig(i,j)=numel(sToUse);
            end
        else
            switch diffModel
                case {'DiDe','DiDe_with_rise_time'}
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'DiDe_unnormalised'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;    
                case 'DiDe_notort'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;    
                case 'DiDe_fiso'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'DiDe_fit_fiso'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    fisoMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;    
                case 'D_fiso'
                    rMap(i,j)=0;
                    dMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'D'
                    rMap(i,j)=0;
                    dMap(i,j)=0;
                    fMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
            end
        end
        else
            switch diffModel
                case {'DiDe','DiDe_with_rise_time'}
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'DiDe_unnormalised'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;    
                case 'DiDe_notort'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'DiDe_fiso'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'DiDe_fit_fiso'
                    rMap(i,j)=0;
                    diMap(i,j)=0;
                    deMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    fisoMap(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'D_fiso'
                    rMap(i,j)=0;
                    dMap(i,j)=0;
                    fiMap(i,j)=0;
                    feMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    extFlgMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
                case 'D'
                    rMap(i,j)=0;
                    dMap(i,j)=0;
                    fMap(i,j)=0;
                    s0Map(i,j)=0;
                    rsqMap(i,j)=0;
                    strtValMap(i,j)=0;
                    numSig(i,j)=0;
                    fitRes(i,j)=NaN;
                    sToUse(i,j)=NaN;
                    paramsToUse(i,j)=NaN;
            end
        end  
    end
end

switch diffModel
    case {'DiDe','DiDe_with_rise_time'}
        fittedParams=cat(3,rMap,diMap,deMap,fMap,...
            rsqMap,strtValMap,numSig);
    case 'DiDe_unnormalised'
        fittedParams=cat(3,rMap,diMap,deMap,fMap,s0Map,...
            rsqMap,strtValMap,numSig);    
    case 'DiDe_notort'
        fittedParams=cat(3,rMap,diMap,deMap,fMap,...
            rsqMap,strtValMap,numSig);    
    case 'DiDe_fiso'
        fittedParams=cat(3,rMap,diMap,deMap,fiMap,feMap,...
            rsqMap,strtValMap,extFlgMap,numSig);
    case 'DiDe_fit_fiso'
        fittedParams=cat(3,rMap,diMap,deMap,fiMap,feMap,fisoMap,...
            rsqMap,strtValMap,extFlgMap,numSig);    
    case 'D_fiso'
        fittedParams=cat(3,rMap,dMap,fiMap,feMap,s0Map,...
            rsqMap,strtValMap,extFlgMap,numSig);    
    case 'D'
        fittedParams=cat(3,rMap,dMap,fMap,s0Map,...
            rsqMap,strtValMap,numSig);
end