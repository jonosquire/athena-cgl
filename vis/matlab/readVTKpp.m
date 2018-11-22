
function V = readVTKpp(vtkfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: V = readVTK(vtkfile)
%
%   V:       Structure in which the data are stored
%   vtkfile: The filename
%   notes:   Only reads binary STRUCTURED_POINTS
%
% Erik Vidholm 2006
% Geoffroy Lesur 2009
%       Extended to include several fields in the same file (as Snoopy
%       does)
%       The output is now a structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% open file (OBS! big endian format)
fid = fopen(vtkfile,'r','b');

if( fid == -1 )
  disp('File not found');
end

fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET RECTILINEAR_GRID

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').';
sz=sz-1;
sz(sz==0)=1;

% Generate coordinates
s = fgetl(fid); % X_COORDINATES nx float
nx = sscanf(s, '%*s%d%*s');
if nx-1==sz(1)
    V.xe = fread(fid,nx,'*single');
    V.x = (V.xe(1:end-1)+V.xe(2:end))/2;
else
    error('Dimensions and x-size dont match')
end
fgetl(fid);
% Generate coordinates
s = fgetl(fid); % Y_COORDINATES nx float
ny = sscanf(s, '%*s%d%*s');
if ny-1==sz(2) || ny==1
    V.ye = fread(fid,ny,'*single');
    V.y = (V.ye(1:end-1)+V.ye(2:end))/2;
else
    error('Dimensions and y-size dont match')
end
fgetl(fid);
% Generate coordinates
s = fgetl(fid); % X_COORDINATES nz float
nz = sscanf(s, '%*s%d%*s');
if nz-1==sz(3) || nz==1
    V.ze = fread(fid,nz,'*single');
    V.z = (V.ze(1:end-1)+V.ze(2:end))/2;
else
    error('Dimensions and z-size dont match')
end
fgetl(fid);

V.nx = sz(1);
V.ny = sz(2);
V.nz = sz(3);

s = fgetl(fid); % CELL_DATA NXNYNZ
n2read = sscanf(s, '%*s%d');

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
varname = sscanf(s, '%*s%s%*s');

fgetl(fid); % LOOKUP_TABLE default

% read data
Q = fread(fid,n2read,'*single');
disp(['Adding ' varname])
V = setfield(V,varname,reshape(Q,[1 sz ]));

%Let's now fight with the FIELD region for the other components.
linesBetweenData = 1;
for a=0:linesBetweenData;s=fgetl(fid); end
szr=sz; % size to read
while s ~= -1

    vtype = sscanf(s, '%s%*s%*s');
    varname = sscanf(s, '%*s%s%*s');
    disp(['Adding ' varname])

    if strcmp(vtype,'VECTORS') 
        ncmpt=3;
        szr = sz;
    elseif strcmp(vtype,'SCALARS') 
        ncmpt=1;
        szr = sz;
        fgetl(fid); % LOOKUP_TABLE default -- only with pressure for some reason?
    else
        warning([varname ' NOT A VECTOR OR SCALAR'])
    end

    % read data
    Q = fread(fid,ncmpt*n2read,'*single');
    V = setfield(V,varname,reshape(Q,[ncmpt szr ]));

    for a=0:linesBetweenData;s=fgetl(fid); end
end


% Put velocity and magnetic field into form that matches hdf5 output
if isfield(V,'vel')
    V.vel1 = V.vel(1,:,:,:);
    V.vel2 = V.vel(1,:,:,:);
    V.vel3 = V.vel(1,:,:,:);
end
if isfield(V,'Bcc')
    V.Bcc1 = V.Bcc(1,:,:,:);
    V.Bcc2 = V.Bcc(1,:,:,:);
    V.Bcc3 = V.Bcc(1,:,:,:);
end
    
fclose(fid);