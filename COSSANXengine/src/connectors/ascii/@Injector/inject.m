function inject(Xi,Tinput)
%%
%  Inject the input values into INPUT ASCII FILE
%
%   Arguments:
%   ==========
%   Xi                  Injector object
%   Pinput          Input Object/Structure
%   varargin       optional parameters
%   varargout     Log results
%
%   Usage:  Clog  = inject(Xi,Xinput)
%
%   see also: connector, extractor, injector
%

%%. Processing Inputs

% Process all the optional arguments and assign them the corresponding
% default value if not passed as argument

assert(isa(Tinput,'struct'),'OpenCossan:injector:inject', ...
    'At least the input values must be passed as a structure');
assert(length(Tinput)==1,'OpenCossan:injector:inject', ...
    'Only one realization of the input parameters is allowed');

%% Check open files

VfIDs = fopen('all');
for ifid=1:length(VfIDs)
    [SfilenameDB, SpermissionDB] = fopen(VfIDs(ifid));
    OpenCossan.cossanDisp(['[OpenCossan.Injector.inject] Open filename: ' SfilenameDB ' permission: ' SpermissionDB ],3 )
end

if isempty(VfIDs)
    OpenCossan.cossanDisp('[OpenCossan.Injector.inject] No currently open files: ;)',4 )
end



%% 2.  Update input files

OpenCossan.cossanDisp(['[OpenCossan.Injector.inject] Current directory: ' pwd ],4 )
OpenCossan.cossanDisp('[OpenCossan.Injector.inject] Current directory using system commands: ' ,4 )
if OpenCossan.getVerbosityLevel>3
    if isunix
        system('pwd');
    else
        system('cd');
    end
end

%% Replace values in the FE input file
SfullName=fullfile(Xi.Sworkingdirectory,Xi.Srelativepath,Xi.Sfile);

OpenCossan.cossanDisp(['[OpenCossan.Injector.inject] Replacing value in the file: ' SfullName ],2)

OpenCossan.cossanDisp('[OpenCossan.Injector.inject] File content before injecting values',4)

if OpenCossan.getVerbosityLevel>3
    type(SfullName)
end

if isa(Xi,'TableInjector')    
    Xi.doInject(Tinput);
else
    
    [Nfid,Serror] = fopen(SfullName,'r+'); % open ASCII file
    
    OpenCossan.cossanDisp(['[OpenCossan.Injector.inject] Open file (relative path): ' ...
        SfullName],3)
    if Nfid~=-1
        OpenCossan.cossanDisp(['[OpenCossan.Injector.inject] Fid: ' num2str(Nfid)],4)
    end
    
    if ~isempty(Serror)
        OpenCossan.cossanDisp('[OpenCossan.Injector.inject] Content of the current directory: ',1)
        error('OpenCossan:Injector:inject',...
            ['Error opening file ' SfullName ' \n Error: ' Serror ]);
    end
    
    if ~isempty(Xi.Xidentifier)
        Xi.Xidentifier.replaceValues('Nfid',Nfid,'Tinput',Tinput);
    else
        warning('OpenCossan:Injector:inject',['No identifiers defined in input file.\n'...
            'Consider put your input file in CadditionalFiles of Connector instead']);
    end
    %% Close open files
    fclose(Nfid);
end


OpenCossan.cossanDisp('[OpenCossan.Injector.inject] File contents after injecting values',4)
if OpenCossan.getVerbosityLevel>3
    type(SfullName)
end

OpenCossan.cossanDisp('[OpenCossan.Injector.inject] Show folder contents',4)
if OpenCossan.getVerbosityLevel>3
    if ispc
        system('dir /A')
    else
        system('dir -all')
    end
end
OpenCossan.cossanDisp('[OpenCossan.Injector.inject] Inject completed',4)

