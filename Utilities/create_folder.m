function create_folder(newdir)
%
% CREATE_FOLDER(NEWDIR) creates a folder called NEWDIR. If there is already
% a non-empty folder named NEWDIR, it will be renamed to NEWDIR-backup (or
% NEWDIR-backup02, NEWDIR-backup03, ...). If NEWDIR is empty, nothing is
% done.
%

% Guilherme C. Beltramini - 2012-Jan-29, 04:10 pm


if length(ls(newdir))==2    % Directory exists and is empty
    
    % Do nothing
    %----------------------------------------------------------------------
    fprintf('Folder %s already exists and is empty\n',newdir)
    fprintf('Nothing was done\n')
    return
    
elseif exist(newdir,'dir')  % Directory exists and is not empty
    
    % Rename existing directory to a non-existing name
    %----------------------------------------------------------------------
    if exist([newdir '-backup'],'dir')
        aux = 2;
        while exist([newdir sprintf('-backup%.2d',aux)],'dir')
            aux = aux + 1;
        end
        s = movefile(newdir,[newdir sprintf('-backup%.2d',aux)]);
    else
        s = movefile(newdir,[newdir '-backup']);
    end
    if ~s
        error('Check folder permission')
    end
    
end

% Create new directory
%--------------------------------------------------------------------------
s = mkdir(newdir);
pause(0.5) % give some time for folder to appear