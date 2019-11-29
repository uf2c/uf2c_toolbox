function selection = select_files(folder,pref,mid,suf,ext,type,path)
%
% SELECTION = SELECT_FILES(FOLDER,PREF,MID,SUF,EXT,TYPE,PATH) selects files
% (if TYPE='files') or FOLDERS (if TYPE='folders') with prefix PREF, middle
% part MID, suffix SUF and/or extension EXT in folder FOLDER. PREF, MID,
% SUF and EXT can be empty strings. Option EXT is ignored if
% TYPE='folders'.
%
% SELECTION is a cell array containing the full path if PATH='path' or only
% the file or folder names if PATH='nopath'.
%

% Guilherme C. Beltramini - 2012-Feb-03, 10:41 am


% Get files with chosen pattern
%--------------------------------------------------------------------------
curr_dir = pwd;
cd(folder)
switch lower(type)
    case 'folders'
        tmp = dir(sprintf('%s*%s*%s',pref,mid,suf));
    case 'files'
        tmp = dir(sprintf('%s*%s*%s.%s',pref,mid,suf,ext));
end
cd(curr_dir)


% Find the directories in the selection
%--------------------------------------------------------------------------
N           = size(tmp,1);
direct      = cell(N,1);
[direct{:}] = deal(tmp.isdir);
direct      = cell2mat(direct);


% Get chosen type
%--------------------------------------------------------------------------
switch lower(type)
    case 'files'
        tmp = tmp(~direct);
    case 'folders'
        tmp = tmp(direct);
    otherwise
        error('Unknown option for TYPE')
end


N         = size(tmp,1);
selection = cell(N,1);
switch lower(path)
    
    % Full path
    %----------------------------------------------------------------------
    case 'path'
        for f=1:N
            selection{f} = fullfile(folder,tmp(f).name);
        end
    
    % Only file name
    %----------------------------------------------------------------------
    case 'nopath'
        [selection{:}] = deal(tmp.name);
        
    otherwise
        error('Unknown option for PATH')
end