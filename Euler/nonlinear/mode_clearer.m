function modified_C = mode_clearer(C,a)
%
% Clears out the modes corresponding to the indices in a
%
% for clearing out unresolved modes,
%
% a = a_tilde = N+1:M
%
% for clearing out dealiasing modes,
%
% a = a_tilde2 = 2*N+1:M
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  C  =  the array whose entries we wish to clear out
%
%  a  =  the indices we wish to eliminate
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  modified_C  =  the C array with the indicies associated with a set to
%                 zero

% copy the array
modified_C = C;

% eliminate all modes with those indices (positive and negative)
modified_C(a,:,:,:,1) = 0;
modified_C(a,:,:,:,2) = 0;
modified_C(a,:,:,:,3) = 0;
modified_C(a,:,:,:,4) = 0;

modified_C(:,a,:,:,1) = 0;
modified_C(:,a,:,:,2) = 0;
modified_C(:,a,:,:,3) = 0;
modified_C(:,a,:,:,4) = 0;

modified_C(:,:,a,:,1) = 0;
modified_C(:,:,a,:,2) = 0;
modified_C(:,:,a,:,3) = 0;
modified_C(:,:,a,:,4) = 0;
