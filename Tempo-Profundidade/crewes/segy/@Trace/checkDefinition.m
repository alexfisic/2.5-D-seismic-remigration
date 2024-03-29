function checkDefinition(obj, bd)
%
% function checkDefinition(obj, bd)
%
% CHECKDEFINITION checks a trace header definition cell array for:
% 1. size (five columns)
% 2. Column 1 (field names)
%     a. cell contents are char
%     b. duplicate field names
%     c. spaces within the the field names
% 3. Column 2 (data types)
%     a. cell contents are char and are one of:
%     b. 'uint8','int8','uint16','int16','uint32','int32','uint64',
%        'int64','ibm32','ieee32','ieee64'
% 4. Column 3 (start byte)
%     a. cell contents are numeric and
%     b. cross-reference column 2 with columns 3 and verify that the total
%        number of defined bytes matche the total number of bytes
%        expected in a SEG-Y binary file header
% 5. Column 4 (co-ordinate scalar field name to apply)
%     a. cell contents are char
%     b. [duplicates are OK]
%     c. spaces within the the field names
% 6. Column 5 (description)
%     a. cell contents are char
%
% Authors: Kevin Hall, 2009, 2017
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

%make sure bd is a cell array
if ~iscell(bd)
    error('@Trace/checkDefinition: Definition must be a cell array')
end

[nfields,ncols] = size(bd);

%make sure we have five columns
if ~isequal(ncols,5)
    error('@Trace/checkDefinition: Definition must have 5 columns')
end

%make sure 1st, 2nd or 4th columns are char
if any(cellfun(@(X) ~ischar(X), bd(:,1))) || ...
        any(cellfun(@(X) ~ischar(X), bd(:,2))) || ...
        any(cellfun(@(X) ~isnumeric(X), bd(:,3))) || ...
        any(cellfun(@(X) ~ischar(X), bd(:,4))) || ...
        any(cellfun(@(X) ~ischar(X), bd(:,5)))
    error('@Trace/checkDefinition: Definition cols 1,2,4,5 must be ''char'', col 3 must be ''numeric''')
end

%check to see if column1 contains duplicate field names
nfieldnames = length(unique(bd(:,1)));
if ~isequal(nfieldnames,nfields)
    error('@Trace/checkDefinition: Definition col 1 contains duplicate field names ')
end

%check to see if column1 field names have spaces in them
c = cellfun(@(X) sum(isspace(X)), bd(:,1),'UniformOutput',false);
if sum(cell2mat(c))
   error('@Trace/checkDefinition: Definition col 1 contains field names with spaces in them')
end

%check to see if column4 field names have spaces in them
c = cellfun(@(X) sum(isspace(X)), bd(:,1),'UniformOutput',false);
if sum(cell2mat(c))
   error('@Trace/checkDefinition: Definition col 4 contains field names with spaces in them')
end

%create a numbered row vector so we can output useful error messages
rownums = 1:nfields;

%convert single to ieee32 and double to ieee64
bd(:,2) = strrep(bd(:,2), 'single', 'ieee32');
bd(:,2) = strrep(bd(:,2), 'double', 'ieee64');

%check to see if column2 contains valid datatypes
vals = {'uint8','int8','uint16','int16','uint32','int32',...
        'uint64','int64','ibm32','ieee32','ieee64'};
vt = false(nfields,1);
for k = 1:length(vals)
    vt = vt | strcmp(bd(:,2),vals(k));
end

if sum(vt) ~= nfields
    br = num2str(rownums(~vt),'%d, ');
    error(['@Trace/checkDefinition: Definition col 2  has invalid datatype(s) at row(s) '...
        br(1:end-1) ])
end

%Check total number of bytes in definition
% convert data type to number of bytes (eg. 'int32' -> 4 bytes)
bytes = cellfun(@str2double,regexp(bd(:,2),'\d*','match'))/8.0;

if mod(sum(bytes),obj.HDRSIZE)
    if isempty(obj.GUI) || isa(obj.GUI,'handle')
        mm_warndlg(['Total number of bytes represented by Definition col 2 is not a multiple of ' ...
            num2str(obj.HDRSIZE)],'Warning!',gcbf);
    elseif objGUI
        warning(['Total number of bytes represented by Definition col 2 is not a multiple of ' ...
            num2str(obj.HDRSIZE)])
    end
end    

cb=cumsum(bytes); %cumulative sum list of data sizes from col2 to compare to start bytes in col3
sbyte = cell2mat(bd(:,3)); %convert col3 to doubles
vt = zeros(nfields,1);
vt(2:end) = cb(1:end-1)-sbyte(2:end); %subtract start bytes from cumulative sum
vt = logical(vt);

if sum(vt)
    br = num2str(rownums(vt),'%d, ');
    error(['@Trace/checkDefinition: Definition col 3 has inconsistent startbyte at row(s) '...
        br(1:end-1) ])
end

end %end function