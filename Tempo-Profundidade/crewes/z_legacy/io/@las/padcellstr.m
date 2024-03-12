function cs = padcellstr(cs)
%
% cs = function padcellstr(cs)
%
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

if ~iscellstr(cs)
    warning('crewes:las:padcellstr','Input must be a cell array of strings')
    cs = [];
    return
end

% Get maximum and minimum lengths of the strings
n = cellfun(@length,cs);
maxlength = max(n);
minlength = min(n);

% Which columns have different minimums and maximums
idx = (minlength ~= maxlength);

% Replicate maxlength for columns that need padding
[m n] = size(cs(:,idx));
len = num2cell(repmat(-maxlength(idx),m,1));

% Pad columns (right justified)
cs(:,idx) = cellfun(@(X,Y) sprintf('%*s',-Y,X),...
            cs(:,idx),len,...
            'UniformOutput',false);

end %end function