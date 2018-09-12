function sortModellingStructs=sortModellingStructs(varargin)
%% Function for sorting separate diffusion data structures so they can be
%  used collectively for microstructural modelling
%  Inputs: variable number of diffusion data structures, with the following fields:
%          diffIms, nSlices, del, DELTA
%  Output: sortModellingStructs - 6D matrix containing diffusion images in the format:
%          x,y,gradient,slice,direction,DEL

%% NB. CURRENTLY ASSUMES ONE DIRECTION, AND THAT EACH ARGUMENT IS A
%      STRUCTURE FOR A DIFFERENT DEL
disp(['!***Currently assumes one direction, and that each argument is a'...
      ' structure for a different DELTA***!'])

%% Check structures have expected fields
for i=1:nargin
	fn=fieldnames(varargin{i});
	if all(ismember({'diffIms','nSlices','DEL','del'},fn))==0
		disp(i)
		error('!!! Structure is missing expected fields')
	end
end

%% Loop over inputs and sort 
sortModellingStructs=[];
for i=1:nargin
    x=[];
    thisStruct=varargin{i};
    ims=thisStruct.diffIms;
    numSlices=thisStruct.nSlices;
    b0=ims(:,:,1:numSlices);
    non_b0=ims(:,:,numSlices+1:end);
    for j=1:numSlices
        x(:,:,:,j)=cat(3,b0(:,:,j),...
            non_b0(:,:,j:numSlices:end));
    end
    sortModellingStructs=cat(6,sortModellingStructs,x);
end