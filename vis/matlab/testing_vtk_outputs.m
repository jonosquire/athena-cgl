% Reads vtk files and calculates energy spectra
% Requries join_vtk.sh has been run if MPI used for original
function testing_vtk_outputs()
folder = '~/Research/athena/turb-tests/decay'; % Folder with outputs
file = 'Turb'; % Name of output
output_id = 2; % Output id (set in input file) 
MHD = 0; % logical for MHD or not

snapshot_num = 20; 

filename = @(n) [folder '/' file '.out' num2str(output_id) '.'  sprintf('%05d',n) '.athdf'];
    
% Number of elements in k grid for spectrum. If numKgrid=0, it uses
% "standard" grid (2pi/L,4pi/L,6pi/L,...)

% Form grid of K using first time step
D = readHDF5(filename(snapshot_num));
dx = D.x(2)-D.x(1);
dy = D.y(2)-D.y(1);
dz = D.z(2)-D.z(1);
Ls = [max(D.x)+dx max(D.y)+dy max(D.z)+dz];
Ns = [length(D.x) length(D.y) length(D.z)];
% Grid in k space, for each dimension, conforming to FT standard [0 1 2 3 ... -N/2 -N/2+1 ... -1]
for kk=1:3
    K{kk} = 2i*pi/Ls(kk)*[0:(Ns(kk)/2-1) -Ns(kk)/2 -Ns(kk)/2+1:-1].';
end
[KX, KY, KZ] = ndgrid(K{1},K{2},K{3}); % 3D grid in K
Kmag = sqrt(abs(KX).^2 + abs(KY).^2 + abs(KZ).^2); % |K|

V = double(D.vel2);
gradV = ifftn(KX.*fftn(double(D.vel1)) + KZ.*fftn(double(D.vel2)) + KY.*fftn(double(D.vel3)));
divV = ifftn(KX.*fftn(double(D.vel1)) + KY.*fftn(double(D.vel2)) + KZ.*fftn(double(D.vel3)));


VolumePlot(real(divV),double(D.x),double(D.y),double(D.z))

end