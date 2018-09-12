function [paramMaps,fitResiduals,signalsUsed,scanParamsUsed,weightsUsed]=...
    prepareForFitting(fitType,dataMatrix,scanParams,initialParamSeed,...
    voxelwiseInitialisation,noiseLevel,lwrBnd,upprBnd,ls_mle,mask,diffModel)
%% Function which organises data and performs normalisation before fitting
% Inputs: fitType - string specifying 'voxelwise' or 'allVoxelsMean'
%         dataMatrix - 6D matrix with diffusion images
%         scanParams - matrix specifying scan parameters
%         initialParamSeed - matrix of initial values to use in fitting
%         voxelwiseInitialisation - string specifying initialisation
%         noiseLevel - Rician noise from background
%         lwrBnd,upprBnd - lower and upper bounds for fit constraints
%         ls_mle - 'ls' or 'mle' fitting
%         mask - binary mask specifying voxels to use when averaging
%         diffModel - string specifying model to fit, 'D' or 'De'
% Outputs: paramMaps - matrix of fitted parameter maps
%          fitResiduals - for allVoxelsAveraged, vector of fit residuals
%          signalsUsed - for allVoxelsAveraged, vector of signals used in fitting
%          scanParamsUsed - for allVoxelsAveraged, matrix of scan parameters used in fitting

% 11/09/18 - check correct package is used
disp('!!! 11/09/18 !!!')

switch fitType
    case 'voxelwise' % Fit voxelwise signals
    
    	% Pre-allocate output
    	switch diffModel
    		case {'D','DiDe','DiDe_notort','DiDe_with_rise_time'}
        		paramMaps=zeros(size(dataMatrix,1),size(dataMatrix,2),...
            		7,size(dataMatrix,4));
            case 'DiDe_unnormalised'
                paramMaps=zeros(size(dataMatrix,1),size(dataMatrix,2),...
            		8,size(dataMatrix,4));
        	case 'DiDe_fiso'
            	paramMaps=zeros(size(dataMatrix,1),size(dataMatrix,2),...
            		9,size(dataMatrix,4));
			case 'DiDe_fit_fiso'
            	paramMaps=zeros(size(dataMatrix,1),size(dataMatrix,2),...
            		10,size(dataMatrix,4));
        	case 'D_fiso'
            	paramMaps=zeros(size(dataMatrix,1),size(dataMatrix,2),...
            		9,size(dataMatrix,4));
        	otherwise
        		error('!!! Unrecognised diffModel')
    	end
        %% Loop over slices and fit voxel wise
        for sliceInd=1:size(dataMatrix,4)
            
            % Get this slice and corresponding noise
            thisSliceAllDir=squeeze(dataMatrix(:,:,:,sliceInd,:,:));
            sliceNoiseLevel=noiseLevel(sliceInd);
            
            % Average across directions
            thisSliceMean=squeeze(mean(thisSliceAllDir,4));
            
            % Normalise to G0?
            switch diffModel
                case 'DiDe_unnormalised'
                    disp('!!! NOT NORMALISING SIGNALS - CHECK !!!')
                    for DELind=1:size(dataMatrix,6)
                        for gradInd=1:size(dataMatrix,3)
                            thisSliceNorm(:,:,gradInd,DELind)=...
                                thisSliceMean(:,:,gradInd,DELind);
                            % weights
                            weightsTmp=(sliceNoiseLevel^2).*(thisSliceNorm(:,:,gradInd,DELind).^2).*( (1./(thisSliceMean(:,:,gradInd,DELind).^2)) + (1./(thisSliceMean(:,:,1,DELind).^2)));
                            weights(:,:,gradInd,DELind)=1./weightsTmp;
                        end
                    end
                otherwise
                    for DELind=1:size(dataMatrix,6)
                        for gradInd=1:size(dataMatrix,3)
                            thisSliceNorm(:,:,gradInd,DELind)=...
                                thisSliceMean(:,:,gradInd,DELind)./...
                                thisSliceMean(:,:,1,DELind);
                            % weights
                            weightsTmp=(sliceNoiseLevel^2).*(thisSliceNorm(:,:,gradInd,DELind).^2).*( (1./(thisSliceMean(:,:,gradInd,DELind).^2)) + (1./(thisSliceMean(:,:,1,DELind).^2)));
                            weights(:,:,gradInd,DELind)=1./weightsTmp;
                        end
                    end
            end
            
            % Are we using voxelwise initialisation?
            switch voxelwiseInitialisation
                case 'voxelInit'
                	% Indices for initialisation
            		switch diffModel
            			case {'DiDe','DiDe_with_rise_time'}
                			initMapInd=1:4;
                        case 'DiDe_unnormalised'
                            initMapInd=1:5;
            			case 'D'
                			initMapInd=1:3;
            			case 'DiDe_notort'
                			initMapInd=1:4;
            		end
                    % fitting
                    paramMaps(:,:,:,sliceInd)=fitSubFn(thisSliceNorm,thisSliceMean,...
                        scanParams,initialParamSeed(:,:,initMapInd,sliceInd),...
                        voxelwiseInitialisation,sliceNoiseLevel,...
                        lwrBnd,upprBnd,diffModel,weights,ls_mle);
                case 'standardInit'
                    % fitting
                    paramMaps(:,:,:,sliceInd)=fitSubFn(thisSliceNorm,thisSliceMean,...
                        scanParams,initialParamSeed,...
                        voxelwiseInitialisation,sliceNoiseLevel,...
                        lwrBnd,upprBnd,diffModel,weights,ls_mle);
            end
        end
        
    case 'ROIav'
        %% Fit slice-wise averaged signals
        
        % First check that the mask was supplied as an argument
        if isempty(mask)
            error('Error - mask file must be supplied for allVoxelsMean fitting')
        end
        
        % Loop over slices
        for sliceInd=1:size(dataMatrix,4)
            
            % Get this slice and corresponding noise
            thisSliceAllDir=squeeze(dataMatrix(:,:,:,sliceInd,:,:));
            sliceNoiseLevel=noiseLevel(sliceInd);
            
            % Apply mask
            thisSliceMask=mask(:,:,sliceInd);
            maskRep=repmat(thisSliceMask,[1 1 ...
                size(thisSliceAllDir,3) size(thisSliceAllDir,4) ...
                size(thisSliceAllDir,5)]);
            applyMask=thisSliceAllDir.*maskRep;
            
            % Average over ROI and directions, giving a mean signal for each
            % G and DEL; take mean of voxels where signal>0 so background isn't
            % included
            for gradInd=1:size(applyMask,3)
                for DELind=1:size(applyMask,5)
                    thisSliceAll=applyMask(:,:,gradInd,:,DELind);
                    thisSliceAllVec=thisSliceAll(:);
                    thisSliceMeanTemp(gradInd,DELind)=...
                        mean(thisSliceAllVec(thisSliceAllVec>0));
                end
            end
            
            % Reshape so it can be passed to same code as voxelwise fitting,
            % which assumes we're working with images
            thisSliceMean=reshape(thisSliceMeanTemp,1,1,size(applyMask,3),...
                size(applyMask,5));
            
            % Normalise to G0
            for DELind=1:size(applyMask,5)
                for gradInd=1:size(applyMask,3)
                    thisSliceNorm(:,:,gradInd,DELind)=...
                        thisSliceMean(:,:,gradInd,DELind)./...
                        thisSliceMean(:,:,1,DELind);
                end
            end
            
            % Fitting
            [paramMaps(:,:,:,sliceInd),residuals{sliceInd},...
                signals{sliceInd},scanParamsUsed{sliceInd}]...
                =fitSubFn(thisSliceNorm,thisSliceMean,...
                scanParams,initialParamSeed,sliceNoiseLevel,...
                lwrBnd,upprBnd,diffModel,ls_mle);
        end
        
    case {'allVoxelsMean','allVoxelsMedian'}
        %% Fit average of all signals (i.e. all voxels in all slices)
        
        % First check that the mask was supplied as an argument
        if isempty(mask)
            error('Error - mask file must be supplied for allVoxelsMean fitting')
        end
        
        % Create empty array to hold all voxel signals
        collectAll=[];
        
        % Average noise level over all slices
        avNoiseLevel=mean(noiseLevel(:));
        
        % Loop over slices
        for sliceInd=1:size(dataMatrix,4)
            
            % Index for collectAll update
            i1=size(collectAll,1)+1;
            
            % Get this slice 
            thisSliceAllDir=squeeze(dataMatrix(:,:,:,sliceInd,:,:));
                                    
            % Apply mask
            thisSliceMask=mask(:,:,sliceInd);
            maskRep=repmat(thisSliceMask,[1 1 ...
                size(thisSliceAllDir,3) size(thisSliceAllDir,4) ...
                size(thisSliceAllDir,5)]);
            applyMask=thisSliceAllDir.*maskRep;
            
            % Collect signals for all directions for this slice, as a
            % function of G and DEL
            for gradInd=1:size(applyMask,3)
                for DELind=1:size(applyMask,5)
                    thisSliceAll=applyMask(:,:,gradInd,:,DELind);
                    i2=i1-1+size(thisSliceAll(:),1);
                    collectAll(i1:i2,gradInd,DELind)=thisSliceAll(:);
                end
            end
        end
        
        % Check that we've collected the expected total number of voxels
        % NOTE: this includes phantom voxels and background!
        numExpectedVoxels=size(dataMatrix,1)*size(dataMatrix,2)*...
            size(dataMatrix,4)*size(dataMatrix,5);
        numCollectedVoxels=size(collectAll,1);
        if numExpectedVoxels==numCollectedVoxels
            % if we're here we're ok
        else
            numExpectedVoxels
            numCollectedVoxels
            error('!!! Error - incorrect total number of voxels collected')
        end
        
        % Average over all voxels
        for gradInd=1:size(applyMask,3)
            for DELind=1:size(applyMask,5)
                thisSigVec=squeeze(collectAll(:,gradInd,DELind));
                thisSigVecMaskedVox=thisSigVec(thisSigVec>0);
                % check we're averaging over the expected number of masked
                % voxels
                if size(thisSigVecMaskedVox,1)==...
                       size(find(mask(:)==1),1)*size(dataMatrix,5)
                    if strcmp(fitType,'allVoxelsMean')
                        meanTemp(gradInd,DELind)=...
                            mean(thisSigVecMaskedVox);
                    elseif strcmp(fitType,'allVoxelsMedian')
                        meanTemp(gradInd,DELind)=...
                            median(thisSigVecMaskedVox);
                    else
                        error('!!! Should not get here!')
                    end
                else
                    size(thisSigVecMaskedVox,1)
                    size(find(mask(:)==1),1)
                    error('!!! Error - incorrect number of masked voxels')
                end
            end
        end
        % Reshape so it can be passed to same code as voxelwise fitting,
        % which assumes we're working with images
        thisMean=reshape(meanTemp,1,1,size(applyMask,3),...
            size(applyMask,5));
        
        % Normalise to G0
        for DELind=1:size(applyMask,5)
            for gradInd=1:size(applyMask,3)
                thisNorm(:,:,gradInd,DELind)=...
                    thisMean(:,:,gradInd,DELind)./...
                    thisMean(:,:,1,DELind);
                weightsTmp=(avNoiseLevel^2).*(thisNorm(:,:,gradInd,DELind).^2).*( (1./(thisMean(:,:,gradInd,DELind).^2)) + (1./(thisMean(:,:,1,DELind).^2)));
                weights(:,:,gradInd,DELind)=1./weightsTmp;
            end
        end
        %{
        testOut.collectAll=collectAll;
        testOut.thisNorm=thisNorm;
        testOut.thisMean=thisMean;
        testOut.avNoiseLevel=avNoiseLevel;    
        testOut.lb=lwrBnd;
        testOut.ub=upprBnd;
        %}
     
        % Fitting
        [paramMaps,residuals,signals,scanParamsUsed,weightsUsed]=fitSubFn(thisNorm,...
            thisMean,scanParams,initialParamSeed,avNoiseLevel,...
            lwrBnd,upprBnd,diffModel,weights,ls_mle);
        
    otherwise
        error('!!! Unrecognised fitType')
end

%% Assign output
if strcmp(fitType,'voxelwise')==1
    fitResiduals=[];
    signalsUsed=[];
    scanParamsUsed=[];
    weightsUsed=[];
else
    fitResiduals=residuals;
    signalsUsed=signals;
end

end

%% Subfunction which reshapes signal matrices and
%  calls fitMM_voxelwise
function [fitSubFn,fitRes,sigUsed,scanParamsUsed,weightsUsed]=...
    fitSubFn(thisSliceNorm,thisSliceMean,...
    scanParams,initialParamSeed,voxelwiseInitialisation,sliceNoiseLevel,...
    lwrBnd,upprBnd,diffModel,weightsForFitting,ls_mle)

% Reshape to 3D matrix - format: (:,:,firstG & firstDEL, secondG &
% firstDEL, thirdG & firstDEL...)
signalsToFitNorm=reshape(thisSliceNorm,size(thisSliceNorm,1),...
    size(thisSliceNorm,2),...
    size(thisSliceNorm,3)*size(thisSliceNorm,4));
signalsToFit=reshape(thisSliceMean,size(thisSliceMean,1),...
    size(thisSliceMean,2),...
    size(thisSliceMean,3)*size(thisSliceMean,4));

switch ls_mle
    case 'ls'
        [fitSubFn,fitRes,sigUsed,scanParamsUsed,weightsUsed]=...
            mm_fitting.fitMM_voxelwise(signalsToFit,signalsToFitNorm,...
            scanParams,initialParamSeed,voxelwiseInitialisation,sliceNoiseLevel,...
            lwrBnd,upprBnd,weightsForFitting,diffModel,'ls');
    case 'mle'
        [fitSubFn,fitRes,sigUsed,scanParamsUsed,weightsUsed]=...
            mm_fitting.fitMM_voxelwise(signalsToFit,signalsToFitNorm,...
            scanParams,initialParamSeed,voxelwiseInitialisation,sliceNoiseLevel,...
            lwrBnd,upprBnd,weightsForFitting,diffModel,'mle');
    otherwise 
        error('!!! ls_mle must be ls or mle')
end
end
