%% Function for processing microstructural modelling data
function exitcode = ...
        process_dwi_mm_multiFixDiff(dirAndBoundInd)
    
    %% Args come in as strings. Convert to numbers:
    dirAndBoundInd=str2num(dirAndBoundInd);
    
    %% Process invivo data or simulated data?
    invivo_or_sims='invivo';%'invivo','sims'
    
    %% Acquisition type
    acqType='DEL10_DEL40';
    
    %% Which which type of fitting
    fitType='voxelwise'; %'voxel' or 'ROIav' or 'allVoxels' or 'allVoxelsMedian'
    
    %% Specify model to fit
    diffModel='DiDe_with_rise_time';%'D' or 'DiDe' or 'DiDe_with_rise_time'
    
    %% Include or exclude high G, high DELTA point?
    highG_highDEL='include';%'exclude' or 'include'
    
    %% Scale variables?
    scaling='unscaled';%'scaled' or 'unscaled'
    
    %% Standard fit initialisation or voxel-wise based on previously-run fit?
    initType='standardInit';%'standardInit','voxelInit';
    
    % if using standardInit, choose to run through fixed diffs, or fit all or fit adc
    standardInitBnd='fitall';%'fitall','fixdiff','fitadc'
    
    % if using voxelInit, provide path to and filename of maps
    pathToInitMaps=...
        ['microstructuralModelling_NewInitialisation_Seed*100valsDiDe_'...
        'unscaled/derivedMaps_fromMultiFixDiff_using16fits/'];
    initMapFile='cfFixedDiffFits.mat';
    
    %% When permuting residuals, permute only within a given DEL, or allow
    %  residuals to move between DELs
    permuteResOption='withinDEL';%'withinDEL' or 'betweenDEL'
    
    %% Use original data or median-filtered data?
    %% 30/03/17 - can now also pass data averages from different repetitions
    dataToUse='unfiltered';%'3x3_apply2';%'unfiltered',
    %'2x2_apply1','2x2_apply2','2x2_apply3',
    %'3x3_apply1','3x3_apply2','3x3_apply3',
    %'4x4_apply1','4x4_apply2','4x4_apply3',
    
    %% Specify masks, if using
    maskFileDir='';%'diff_DEL10_del5_G_0_113_207_293/';
    maskFile={''};
    
    %% Specific information for each scan
    switch invivo_or_sims
        case 'invivo'
            % Specify data directories to use
            % 11/09/18 - just use one test dataset
            pathToDataList={'./20180622_112000_180127_1_1/Analysis_180127_1_1/'};
        case 'sims'
            %Simulation directories
            pathToDataList={...
                'R5_Di0.5_De0.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De0.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De1.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De1.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De1_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De1_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De2_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.5_De2_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.751_De1.2019_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.751_De1.2019_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.8205_De0.667_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.8205_De0.667_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.8349_De1.6302_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di0.8349_De1.6302_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.0808_De1.4818_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.0808_De1.4818_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.3873_De1.1778_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.3873_De1.1778_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.4992_De1.3733_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.4992_De1.3733_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De0.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De0.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De1.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De1.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De1_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De1_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De2_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.5_De2_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.7091_De1.7407_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.7091_De1.7407_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.843_De1.3526_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1.843_De1.3526_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De0.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De0.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De1.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De1.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De1_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De1_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De2_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di1_De2_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De0.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De0.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De1.5_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De1.5_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De1_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De1_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De2_fi30_t2125_snr25.3423_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                'R5_Di2_De2_fi30_t2125_snr253.4233_paramScheme_preclin_DEL9Pt862_DEL40_numAverages_1/';
                };
        otherwise
            error('!!! Unrecognised invivo_or_sims')
    end
    
    %% Generate cell linking all subjects with all DiDe combinations (if fitting with multiFixedDiff)
    % or just a cell of subjects if fitting all
    % Index into cell with dirAndBoundInd
    switch initType
        case 'standardInit'
            switch standardInitBnd
                case 'fixdiff'
                    fixDiList=linspace(0.1,3.0,7).*1e-9;
                    fixDeList=linspace(0.1,3.0,7).*1e-9;
                    fixDiffMat=combvec(fixDiList,fixDeList)';
                    
                    % Group together so that dirAndBoundInd number iterates over all
                    % combinations of data directories and fixed diffusivities
                    groupDirAndBounds=[];
                    c=0;
                    for i=1:numel(pathToDataList)
                        for j=1:size(fixDiffMat,1)
                            c=c+1;
                            groupDirAndBounds{c,1}=pathToDataList(i);
                            groupDirAndBounds{c,2}=fixDiffMat(j,:);
                        end
                    end
                    
                    % Uncomment below if we need to rerun a subset of subjects/diffusivities
                    % indToInc=[6 28 37 57 59 60 72 73 86 87 130 131 136 261 262 264 265 267 273 280 281 282];
                    % groupDirAndBounds=...
                    %	groupDirAndBounds(ismember(1:size(groupDirAndBounds,1),indToInc),:);
                    
                    X=groupDirAndBounds(dirAndBoundInd,:);
                    pathToData=X{1};
                    fixDi=X{2}(1);
                    fixDe=X{2}(2);
                case {'fitall','fitadc'}
                    groupDirAndBounds=[];
                    c=0;
                    for i=1:numel(pathToDataList)
                        c=c+1;
                        groupDirAndBounds{c,1}=pathToDataList(i);
                    end
                    X=groupDirAndBounds(dirAndBoundInd,:);
                    pathToData=X{1};
            end
            
        case 'voxelInit'
            groupDirAndBounds=[];
            c=0;
            for i=1:numel(pathToDataList)
                c=c+1;
                groupDirAndBounds{c,1}=pathToDataList(i);
            end
            X=groupDirAndBounds(dirAndBoundInd,:);
            pathToData=X{1};
    end
    
    % Check that number of days matches number of masks
    if numel(pathToData)~=numel(maskFile)
        error('!!! Number of paths does not match number of mask files')
    end
    
    %% Loop over scans-----------------------------------------------------
    for scanInd=1:numel(pathToData)
        clearvars -except scanInd maskFile maskFileDir pathToData fixDi fixDe ...
            standardInitBnd initType dataToUse permuteResOption initMapFile ...
            pathToInitMaps highG_highDEL diffModel fitType acqType ...
            invivo_or_sims scaling
        
        pathToPhntmMask=fullfile(pathToData{scanInd},maskFileDir,maskFile{scanInd});
        
        %% Load data------------------------------------------------------
        %  - either unfiltered (i.e. original) dataset, or filtered dataset with
        %    specified neighbourhood and number of applications
        switch dataToUse
            case 'unfiltered'
                dataFileName='data.mat';
                noiseFileName='ricianNoiseStd_diff.mat';
            otherwise
                dataFileName=strcat('data_',dataToUse,'.mat');
                noiseFileName=strcat('ricianNoiseStd_diff_',dataToUse,'.mat');
        end
        
        switch invivo_or_sims
            case 'invivo'
                del10=loadRename(strcat(pathToData{scanInd},...
                    'diff_DEL10_del5_G_0_113_207_293/',dataFileName));
                del40=loadRename(strcat(pathToData{scanInd},...
                    'diff_DEL40_del5_G_0_113_207_293/',dataFileName));
                noiseLevel=loadRename(strcat(pathToData{scanInd},...
                    'diff_DEL10_del5_G_0_113_207_293/',noiseFileName));
            case 'sims'
                del10=loadRename(strcat(pathToData{scanInd},...
                    'diffData_scan1/',dataFileName));
                del40=loadRename(strcat(pathToData{scanInd},...
                    'diffData_scan2/',dataFileName));
                noiseLevel=loadRename(strcat(pathToData{scanInd},...
                    'diffData_scan1/',noiseFileName));
            otherwise
                error('!!! Unrecognised invivo_or_sims')
        end
        
        saveDir=pathToData{scanInd};
        
        %% Scan params
        protonGamma=2*pi*42.57746778e6;
        switch acqType
            case 'DEL10_DEL40'
                gradList=[0 113.461484 207.151381 292.956293]./1e3;
                DEL=[9.862 40]./1e3;
                del=4.65/1e3;
                scanParams=combvec(gradList,DEL,del,protonGamma)';
                switch diffModel
                    case 'DiDe_with_rise_time'
                        scanParams(:,5)=0.245e-3;
                    otherwise
                        % Nothing to do
                end
            otherwise
                error('!!! Unrecognised acqType')
        end
        pathToData{scanInd}
        
        %% Sort data for fitting---------------------------------
        %  - Format - x,y,gradient,slice,direction,DELTA
        switch acqType
            case 'DEL10_DEL40'
                dataMatrixLoad=mm_fitting.sortModellingStructs(del10,del40);
        end
        
        switch invivo_or_sims
            case 'invivo'
                % Roughly mask out background
                m=squeeze(dataMatrixLoad(:,:,1,:,1,1))>15;
                for i=1:size(dataMatrixLoad,3)
                    for j=1:size(dataMatrixLoad,5)
                        for k=1:size(dataMatrixLoad,6)
                            dataMatrix(:,:,i,:,j,k)=...
                                squeeze(dataMatrixLoad(:,:,i,:,j,k)).*m;
                        end
                    end
                end
            case 'sims'
                % No masking for simulated data, so just assign dataMatrix
                dataMatrix=dataMatrixLoad;
        end
        % 11/09/18 - just test on a few voxels
        disp('Testing specific voxels')
        dataMatrix(:,:,:,[1:10 12:end],:,:)=0;
        dataMatrix([1:39 42:end],:,:,:,:)=0;
        dataMatrix(:,[1:29 32:end],:,:,:)=0;
        
        %% Remove highest-G, highest-DEL data from fitting, if needed
        switch highG_highDEL
            case 'exclude'
                dataMatrix(:,:,end,:,:,end)=...
                    zeros(size(dataMatrix,1),size(dataMatrix,2),size(dataMatrix,4))+eps;
            case 'include'
                % nothing to do
            otherwise
                error('!!! Unrecognised highG_highDEL')
        end
        
        %% Initialisation
        switch diffModel
            case 'D'
                error('Model with single diffusivitity is not implemented here')
            case {'DiDe','DiDe_notort','DiDe_with_rise_time'}
                rLwrLim=0.1e-6;
                rUpprLim=25e-6;
                diLwrLim=0.1e-9;
                diUpprLim=3e-9;
                deLwrLim=0.1e-9;
                deUpprLim=3e-9;
                fLwrLim=0.01;
                fUpprLim=1;
                
                % Initialisation matrix not generated here ...
                %{
        	initialParamsSeedOrMatrix=134564;
        	radRange=[rLwrLim rUpprLim];
        	diffIRange=[diLwrLim diUpprLim];
        	diffERange=[deLwrLim deUpprLim];
        	fRange=[fLwrLim fUpprLim];%[0 1];%
        	RandStream.setGlobalStream ...
            	(RandStream('mt19937ar','seed',initialParamsSeedOrMatrix));
        	for i=1:100
            initialParams(1,i)=(radRange(1) + (radRange(2)-radRange(1)).*rand(1,1));
            initialParams(2,i)=(diffIRange(1) + (diffIRange(2)-diffIRange(1)).*rand(1,1));
            initialParams(3,i)=(diffERange(1) + (diffERange(2)-diffERange(1)).*rand(1,1));
            initialParams(4,i)=(fRange(1) + (fRange(2)-fRange(1)).*rand(1,1));
        	end
                %}
                
                % ... load pre-generated matrix instead
                initialParamSeed=...
                    loadRename('./initialParams.mat');
                
            case 'DiDe_fiso'
                error('DiDe_fiso not implemented - check!')
                %{
        	rLwrLim=0.1e-6;
        	rUpprLim=25e-6;
        	diLwrLim=0.1e-9;
        	diUpprLim=3e-9;
        	deLwrLim=0.1e-9;
        	deUpprLim=3e-9;
        	fiLwrLim=0;
        	fiUpprLim=1;
        	feLwrLim=0;
        	feUpprLim=1;
        
        	initialParamsSeedOrMatrix=134564;
        	radRange=[rLwrLim rUpprLim];
        	diffIRange=[diLwrLim diUpprLim];
        	diffERange=[deLwrLim deUpprLim];
        	fiRange=[fiLwrLim fiUpprLim];
        	feRange=[feLwrLim feUpprLim];%[0 1];%
        	RandStream.setGlobalStream ...
            	(RandStream('mt19937ar','seed',initialParamsSeedOrMatrix));
        	for i=1:100
            initialParams(1,i)=(radRange(1) + (radRange(2)-radRange(1)).*rand(1,1));
            initialParams(2,i)=(diffIRange(1) + (diffIRange(2)-diffIRange(1)).*rand(1,1));
            initialParams(3,i)=(diffERange(1) + (diffERange(2)-diffERange(1)).*rand(1,1));
            initialParams(4,i)=(fiRange(1) + (fiRange(2)-fiRange(1)).*rand(1,1));
            
            switch feInitialisation
                case 'fe_useLimits'
                    initialParams(5,i)=(feRange(1) + (feRange(2)-feRange(1)).*rand(1,1));
                case 'fe_1minusFi'
                    initialParams(5,i)=1-initialParams(4,i);
                case 'fe_betweenLwrLimAnd1minusFi'
                    initialParams(5,i)=(feRange(1) + (1-initialParams(4,i)-feRange(1)).*rand(1,1));
            end
        	end
        initialParamSeed=initialParams;
                %}
            case 'DiDe_fit_fiso'
                error('DiDe_fit_fiso not implemented - check!')
                %{
        rLwrLim=0.1e-6;
        rUpprLim=25e-6;
        diLwrLim=0.1e-9;
        diUpprLim=3e-9;
        deLwrLim=0.1e-9;
        deUpprLim=3e-9;
        fiLwrLim=0;
        fiUpprLim=1;
        feLwrLim=0;
        feUpprLim=1;
        fisoLwrLim=0;
        fisoUpprLim=1;
        
        initialParamsSeedOrMatrix=134564;
        radRange=[rLwrLim rUpprLim];
        diffIRange=[diLwrLim diUpprLim];
        diffERange=[deLwrLim deUpprLim];
        fiRange=[fiLwrLim fiUpprLim];
        feRange=[feLwrLim feUpprLim];%[0 1];%
        RandStream.setGlobalStream ...
            (RandStream('mt19937ar','seed',initialParamsSeedOrMatrix));
        for i=1:100
            initialParams(1,i)=(radRange(1) + (radRange(2)-radRange(1)).*rand(1,1));
            initialParams(2,i)=(diffIRange(1) + (diffIRange(2)-diffIRange(1)).*rand(1,1));
            initialParams(3,i)=(diffERange(1) + (diffERange(2)-diffERange(1)).*rand(1,1));
            initialParams(4,i)=(fiRange(1) + (fiRange(2)-fiRange(1)).*rand(1,1));
            % Pick fe between its lower lim and 1-fi's initial value
            initialParams(5,i)=(feRange(1) + (1-initialParams(4,i)-feRange(1)).*rand(1,1));
            % Initialise fiso to 1-fi-fe
            initialParams(6,i)=1-initialParams(4,i)-initialParams(5,i);
        end
        initialParamSeed=initialParams;
                %}
            otherwise
                error('!!! Unrecognised diffModel')
        end
        %{%
        
        %% If using voxelwise initialisation, overwrite above assignment
        switch initType
            case 'voxelInit'
                initialParams=[];
                initialParamSeed=[];
                initfile=rdir(fullfile(pathToData{scanInd},pathToInitMaps,initMapFile));
                %check there's only one
                if numel(initfile)~=1
                    error('!!! Either none or more than 1 initfile found')
                else
                    initialParamSeed=loadRename(initfile.name);
                end
            case 'standardInit'
                %nothing to do
        end
        
        %% Create bounds for fit
        switch diffModel
            case 'D'
                error('!!! Model with single diffusivitity is not implemented here')
            case {'DiDe','DiDe_notort','DiDe_with_rise_time'}
                switch initType
                    case 'standardInit'
                        switch standardInitBnd
                            case 'fixdiff'
                                % Get fixed Di and De from matrix above
                                lwrBnd=[rLwrLim fixDi fixDe fLwrLim];
                                upprBnd=[rUpprLim fixDi fixDe fUpprLim];
                            case 'fitall'
                                lwrBnd=[rLwrLim diLwrLim deLwrLim fLwrLim];
                                upprBnd=[rUpprLim diUpprLim deUpprLim fUpprLim];
                            case 'fitadc'
                                lwrBnd=[rLwrLim diLwrLim deLwrLim 0];
                                upprBnd=[rUpprLim diUpprLim deUpprLim 0];
                        end
                    case 'voxelInit'
                        lwrBnd=[rLwrLim diLwrLim deLwrLim fLwrLim];%
                        upprBnd=[rUpprLim diUpprLim deUpprLim fUpprLim];
                end
                
            case 'DiDe_fiso'
                error('DiDe_fiso not implemented - check!')
            case 'DiDe_fit_fiso'
                error('DiDe_fit_fiso not implemented - check!')
            otherwise
                error('!!! Unrecognised diffModel')
        end
        
        %% Save options
        switch diffModel
            case {'DiDe','DiDe_notort','DiDe_with_rise_time'}
                if strcmp(fitType,'ROIav')==1 || strcmp(fitType,'allVoxelsMean')==1 ...
                        || strcmp(fitType,'allVoxelsMedian')==1
                    saveStr=strcat('data_',dataToUse,'_',fitType,'_',...
                        maskFile{scanInd},'_','LB_',...
                        num2str(lwrBnd(1)),'_',num2str(lwrBnd(2)),'_',...
                        num2str(lwrBnd(3)),'_',num2str(lwrBnd(4)),'_UB_',...
                        num2str(upprBnd(1)),'_',num2str(upprBnd(2)),'_',...
                        num2str(upprBnd(3)),'_',num2str(upprBnd(4)));
                elseif strcmp(fitType,'voxelwise')==1
                    saveStr=strcat('data_',dataToUse,'_',fitType,'_','LB_',...
                        num2str(lwrBnd(1)),'_',num2str(lwrBnd(2)),'_',...
                        num2str(lwrBnd(3)),'_',num2str(lwrBnd(4)),'_UB_',...
                        num2str(upprBnd(1)),'_',num2str(upprBnd(2)),'_',...
                        num2str(upprBnd(3)),'_',num2str(upprBnd(4)));
                end
            case 'DiDe_fiso'
                error('DiDe_fiso not implemented - check!')
            case 'DiDe_fit_fiso'
                error('DiDe_fit_fiso not implemented - check!')
        end
        
        disp(saveStr);
        
        if strcmp(diffModel,'DiDe_fiso')
            fisoOpt=strcat('_',feInitialisation);
        else
            fisoOpt=[];
        end
        
        switch acqType
            case 'DEL10_DEL40'
                switch initType
                    case 'standardInit'
                        subDir=strcat(saveDir,'mm_seed',...
                            '_',num2str(size(initialParamSeed,2)),'vals',...
                            diffModel,fisoOpt,'_',scaling,'/',saveStr,...
                            initType,'_/');
                    case 'voxelInit'
                        subDir=strcat(saveDir,'mm_',...
                            diffModel,fisoOpt,'_',scaling,'/',saveStr,...
                            initType,'_/');
                end
        end
        
        if ~exist(subDir,'dir')
            mkdir(subDir)
        else
            disp('Directory already exists')
        end
        
        %% For consistency with previous code, replicate data along fifth
        %  dimension, so that it appears there are 3 directions
        %  This simply means matrix dimensions are consistent (so
        %  squeeze() doesn't collape more than one dimension), and doesn't
        %  affect the fitting
        dataMatrixRep=cat(5,dataMatrix,dataMatrix,dataMatrix);
        size(dataMatrixRep)
        
        %% Deal with scaling if necessary
        switch scaling
            case 'scaled' % scale tissue properties, bounds and scan parameters
                switch diffModel
                    case {'DiDe','DiDe_with_rise_time'}
                        modelParamScalings=[1e6 1e9 1e9 1];
                    case 'DiDe_notort'
                        modelParamScalings=[1e6 1e9 1e9 1];
                end
                switch initType
                    case 'voxelInit'
                        % loop over parameter maps and scale each model param
                        for p=1:numel(modelParamScalings)
                            initialParamSeedScaled(:,:,p,:)=...
                                squeeze(initialParamSeed(:,:,p,:)).*modelParamScalings(p);
                        end
                    case 'standardInit'
                        initialParamSeedScaled=bsxfun(@times,initialParamSeed,modelParamScalings');
                end
                lwrBndScaled=lwrBnd.*modelParamScalings;
                upprBndScaled=upprBnd.*modelParamScalings;
                switch diffModel
                    case 'DiDe_with_rise_time'
                        scanParamsScaled=bsxfun(@times,scanParams,[1e3 1e3 1e3 1e-12 1e3]);
                    otherwise
                        scanParamsScaled=bsxfun(@times,scanParams,[1e3 1e3 1e3 1e-12]);
                end
            case 'unscaled'
                % Reassign, so correct values are used on each loop
                initialParamSeedScaled=initialParamSeed;
                lwrBndScaled=lwrBnd;
                upprBndScaled=upprBnd;
                scanParamsScaled=scanParams;
        end
        
        %% Fitting
        switch fitType
            case {'allVoxelsMean','ROIav','allVoxelsMedian'}
                phntmMaskTmp=loadRename(pathToPhntmMask);
                phntmMask=phntmMaskTmp;
                [paramMaps,fitRes,sigs,scanParamsUsed]=...
                    mm_fitting.prepareForFitting(fitType,...
                    dataMatrixRep,scanParamsScaled,initialParamSeedScaled,initType,...
                    noiseLevel,lwrBndScaled,upprBndScaled,'ls',phntmMask,diffModel);
            case 'voxelwise'
                [paramMaps,fitRes,sigs,scanParamsUsed]=...
                    mm_fitting.prepareForFitting(fitType,...
                    dataMatrixRep,scanParamsScaled,initialParamSeedScaled,initType,...
                    noiseLevel,lwrBndScaled,upprBndScaled,'ls',[],diffModel);
            otherwise
                error('!!! Unrecognised fitType - should not get here!')
        end
        
        %% Permute residuals for bootstrapping
        switch acqType
            case 'DEL10_DEL40'
                if strcmp(fitType,'allVoxelsMean') || strcmp(fitType,'allVoxelsMedian')
                    fitResPerm=permuteResiduals(fitRes,4000,50);
                end
        end
        
        %% Save output
        %{%
        if ~exist(strcat(subDir,'paramMaps_',saveStr,'.mat'),'file')
            save(strcat(subDir,'paramMaps_',saveStr,'.mat'),'paramMaps')
            save(strcat(subDir,'fitResiduals.mat'),'fitRes')
            save(strcat(subDir,'fitSignals.mat'),'sigs')
            if strcmp(fitType,'allVoxels') || ...
                    strcmp(fitType,'allVoxelsMedian')
                save(strcat(subDir,'fitResidualsPerm.mat'),'fitResPerm')
            end
            save(strcat(subDir,'fitScanParams.mat'),'scanParamsUsed')
            save(strcat(subDir,'fitLowerBounds.mat'),'lwrBnd')
            save(strcat(subDir,'fitUpperBounds.mat'),'upprBnd')
        else
            disp(['paramMaps variable already exists - '...
                'not saving any output'])
        end
        %}%
        
        %% display fitted parameters if ROI-fitting
        switch fitType
            case {'ROIav','allVoxels','allVoxelsMedian'}
                switch scaling
                    case 'scaled'
                        squeeze(paramMaps)
                    case 'unscaled'
                        switch diffModel
                            case {'DiDe','DiDe_notort','DiDe_with_rise_time'}
                                squeeze(paramMaps).*[1e6 1e9 1e9 1 1 1 1]'
                            case 'DiDe_fiso'
                                squeeze(paramMaps).*[1e6 1e9 1e9 1 1 1 1 1]'
                            case 'DiDe_fit_fiso'
                                squeeze(paramMaps).*[1e6 1e9 1e9 1 1 1 1 1 1 1]'
                        end
                end
            otherwise
        end
    end
    
    %% Output
    exitcode=1;
