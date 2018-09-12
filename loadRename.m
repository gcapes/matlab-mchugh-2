function loadRename = loadRename(file)
%% function which loads a matlab file, and renames the variable based on 
%  the assignment with which the function is called, e.g.:
%  x=loadRename('./thisFile') - loads the variable saved in 'thisFile', and
%  renames it 'x' in the workspace
%  NOTE: ONLY WORKS FOR FILES CONTAINING ONE VARIABLE

%% Load the variable into a structure
initialStruct=load(file);

%% Get variable name from structure fieldname, and check there is only one 
%  variable
fn=fieldnames(initialStruct);
if numel(fn)~=1
    disp(['!!! Warning - loadRename designed for files in which one ' ...
        'variable is saved - WILL LOAD FIRST FILE'])
else
    %% Assign output, so that the variable is now named X, where X is the 
	%  assignment which was used when loadRename was called
	loadRename=initialStruct.(fn{1});
end