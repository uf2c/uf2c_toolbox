function varargout = RenameFiles(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% NCA - Neuroimaging Computational Analisys %%%%%%%%%%%%%%%%%%%%
%%%%%% Brunno Machado de Campos - University of Campinas - 2012

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RenameFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @RenameFiles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function RenameFiles_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = RenameFiles_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addfiles_Callback(hObject, eventdata, handles)
global nn file path
[file,path] = uigetfile({'*.*','ALL FILES!'},'Add the Files','MultiSelect','on');
file = cellstr(file);
nn = size(file,2);

set(handles.addfilestxt,'String', sprintf('%d File(s) added!',nn))

function addprefixtxt_Callback(hObject, eventdata, handles)
if isequal(get(handles.addprefixtxt,'Value'),1)
    set(handles.addprefix,'Enable','on')
else
    set(handles.addprefix,'Enable','off')
    set(handles.addprefix,'String','')
end

function addprefix_Callback(hObject, eventdata, handles)

function espaco_Callback(hObject, eventdata, handles)

function addsufixtxt_Callback(hObject, eventdata, handles)
if isequal(get(handles.addsufixtxt,'Value'),1)
    set(handles.addsufix,'Enable','on')
else
    set(handles.addsufix,'Enable','off')
    set(handles.addsufix,'String','')
end

function addsufix_Callback(hObject, eventdata, handles)

function jumpertxt_Callback(hObject, eventdata, handles)

if isequal(get(handles.jumpertxt,'Value'),1)
    set(handles.jumper,'Enable','on')
else
    set(handles.jumper,'Enable','off')
    set(handles.jumper,'String','')
end

function jumper_Callback(hObject, eventdata, handles)

function stoptxt_Callback(hObject, eventdata, handles)
if isequal(get(handles.stoptxt,'Value'),1)
    set(handles.stop,'Enable','on')
else
    set(handles.stop,'Enable','off')
    set(handles.stop,'String','')
end

function stop_Callback(hObject, eventdata, handles)

function checkinitials_Callback(hObject, eventdata, handles)
if isequal(get(handles.checkinitials,'Value'),1)
    set(handles.espacotxt,'Enable','on')
    set(handles.espaco,'Enable','on')
    set(handles.COL,'Enable','on')
    set(handles.CLN,'Enable','on')
else
    set(handles.espacotxt,'Enable','off')
    set(handles.espaco,'Enable','off')
    set(handles.COL,'Enable','off')
    set(handles.CLN,'Enable','off')
end

function runboton_Callback(hObject, eventdata, handles)
global nn file path
h = waitbar(0,sprintf('Please wait...%d files selected',nn));
recovery = cell(2,nn);

horas = clock;
namefinal = sprintf('%s%d-%d-%d__%d-%d-%.f_RecoveryNames',...
   path,horas(3),horas(2),horas(1),horas(4),horas(5),horas(6));
   
fide = fopen(sprintf('%s.txt',namefinal),'w+');

order = numel(num2str(nn));


Y=1;

for i = 1:nn
    waitbar(i/nn,h)
    file2Compl = fullfile(path,file{i});
    recovery{i,1} = file2Compl;
    
    [xd,xn,xext] = fileparts(file{i});
    
    newname = xn;
    clear addprefix addsufix
    
    if get(handles.addprefixtxt,'Value')
        addprefix = get(handles.addprefix,'String');
        if get(handles.crescenteprefix,'Value')
            addprefix = [eval(['sprintf(''%.' num2str(order) 'd'',i)']) '-' addprefix];
        end
    else
        if get(handles.crescenteprefix,'Value')
            addprefix = [eval(['sprintf(''%.' num2str(order) 'd'',i)']) '-'];
        else
            addprefix = '';
        end
    end
    
    if get(handles.addsufixtxt,'Value')
        addsufix = get(handles.addsufix,'String');
        if get(handles.crescentesufix,'Value')
            addsufix = [addsufix '-' eval(['sprintf(''%.' num2str(order) 'd'',i)'])];
        end
    else
        if get(handles.crescentesufix,'Value')
            addsufix = ['-' eval(['sprintf(''%.' num2str(order) 'd'',i)'])];
        else
            addsufix = '';
        end
    end

    if get(handles.jumpertxt,'Value') &&  ~isempty(get(handles.jumper,'String'))%Verifica se deseja-se pular alguns caracteres iniciais
        newname = newname(str2double(get(handles.jumper,'String'))+1:end); 
    end
    
    if get(handles.delfend,'Value') &&  ~isempty(get(handles.delendN,'String'))%Verifica se deseja-se pular alguns caracteres iniciais
        newname = newname(1:end-str2double(get(handles.delendN,'String'))); 
    end
    
    if get(handles.stoptxt,'Value') && ~isempty(get(handles.stop,'String'))
        Stopstr = strfind(newname,get(handles.stop,'String'));
        if ~isempty(Stopstr)
            newname = newname(1:Stopstr-1);
        else
            fprintf('The text ''%s'' no longer exist at this point.\r\n',get(handles.stop,'String'))
        end
    end

    if get(handles.starttxt,'Value') && ~isempty(get(handles.start,'String'))
        Startstr = strfind(newname,get(handles.start,'String'));
        if ~isempty(Startstr)
            newname = newname(Startstr+size(get(handles.start,'String'),2):end);
        else
            fprintf('The text ''%s'' no longer exist at this point.\r\n',get(handles.start,'String'))
        end
    end
    
    if ~isempty(get(handles.findreplace,'String')) %remove string ou substitui
        if ~isempty(get(handles.replacetxt,'String'))   %se este não estiver vazio então haverá substituição. o string será substituido
            newname = regexprep(newname,get(handles.findreplace,'String'),get(handles.replacetxt,'String'));
        else
            newname = regexprep(newname,get(handles.findreplace,'String'),''); %se o campo replace estiver vazio, então o string será somente removido.
        end
    end
    
    if get(handles.checkinitials,'Value') && ~isempty(get(handles.espaco,'String'))%Entra na opção onde se quer manter o nome (não se adiciona separadores), mas pode se adicionar sufixo e/ou prefixo jumper...stop....
        
        if get(handles.COL,'Value')
            textInc = 'alpha';
        else
            textInc = 'alphanum';
        end
        
        Espaco = strfind(newname,get(handles.espaco,'String')); %define as posições dos separadores descontanto os que estão após o texto de pausa (se adicionado) e exatamente antes dele.
        if ~isempty(Espaco)
            newnameTMP = newname(1);
            for ir = 1:size(Espaco,2)
                if Espaco(ir)+1<=numel(newname)
                    newnameTMP = [newnameTMP newname(Espaco(ir)+1)];
                end
            end
            TF = isstrprop(newnameTMP,textInc);
            newname = newnameTMP(find(TF)); 
        else
            fprintf('The name separator defined (''%s'') no longer exist.\r\n',get(handles.espaco,'String'))
        end
    end
    
    aVer = [path,addprefix,newname,addsufix,xext]; %simula o nome final

    if exist(aVer,'file')  % Caso haja duas pessoas com iniciais iguais ele atribui à segunda o nome completo + as iniciais.
        newname = [path,addprefix,newname,addsufix,'--',xn,xext]; %atribui o primeiro nome
    end

    filenamefinal = [path,addprefix,newname,addsufix,xext];   
    movefile(file2Compl,filenamefinal,'f')

    recovery{i,2} = filenamefinal;

    fprintf(fide,'%d - %s%s \r\n    %s\r\n\r\n',i,path,xn,filenamefinal);   % salva as alteraçoes em um txt.
end

close(h)

save(sprintf('%s.mat',namefinal),'recovery') %salava as alteraçoes em um .mat que poderá ser usado para restaurar os nomes antigos
fclose(fide);

function recoveryopen_Callback(hObject, eventdata, handles)
wtf = warndlg('The files must to be in the same folder that the program saved them!','Attention!');

waitfor(wtf);
[file2,path2] = uigetfile({'*.mat','Matlab File'},'Add the "RecoveryNames" file created by the program','MultiSelect','off');
load(fullfile(path2,file2));
nn2 = size(recovery,1);
h = waitbar(0,sprintf('Please wait...%d files selected',nn2));
for i = 1:nn2
   if isempty(recovery{i,1})
       break
   end
   waitbar(i/nn2,h)
   movefile(recovery{i,2},recovery{i,1},'f');
end
close(h)

function findreplace_Callback(hObject, eventdata, handles)


function starttxt_Callback(hObject, eventdata, handles)
if get(handles.starttxt,'Value')
    set(handles.start,'Enable','on')
else
    set(handles.start,'Enable','off')
    set(handles.start,'String','')
end

function start_Callback(hObject, eventdata, handles)

function crescenteprefix_Callback(hObject, eventdata, handles)

function replacetxt_Callback(hObject, eventdata, handles)

function crescentesufix_Callback(hObject, eventdata, handles)

function delfend_Callback(hObject, eventdata, handles)
if get(handles.delfend,'Value')
    set(handles.delendN,'Enable','on')
else
    set(handles.delendN,'Enable','off')
    set(handles.delendN,'String','')
end
    
function delendN_Callback(hObject, eventdata, handles)

function COL_Callback(hObject, eventdata, handles)
if get(handles.COL,'Value')
    set(handles.CLN,'Value',0)
else
    set(handles.CLN,'Value',1)
end

function CLN_Callback(hObject, eventdata, handles)
if get(handles.CLN,'Value')
    set(handles.COL,'Value',0)
else
    set(handles.COL,'Value',1)
end

function delendN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
