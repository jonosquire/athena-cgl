% Reads vtk files and calculates energy spectra
% Requries join_vtk.sh has been run if MPI used for original
function spectrum_vtk()
folder = '~/Research/athena/turb-tests/decay'; % Folder with outputs
file = 'Turb'; % Name of output
output_id = 2; % Output id (set in input file) 
MHD = 0; % logical for MHD or not

snapshot_num = 20; 

filename = @(n) [folder '/' file '.out' num2str(output_id) '.'  sprintf('%05d',n) '.athdf'];
    
% Number of elements in k grid for spectrum. If numKgrid=0, it uses
% "standard" grid (2pi/L,4pi/L,6pi/L,...)
numKgrid = 0; 

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
Kmag = sqrt(abs(KX).^2*0 + abs(KY).^2 + abs(KZ).^2); % |K|
% Bins for k
if numKgrid == 0
    kgrid = (0:2*pi/Ls(2):max(imag(K{2}))).'+1e-4; % Use ky for spectrum binning
else
    % Equal spacing in log space for grid.
    kgrid = logspace(log10(min(imag(K{1}))/2),max(min(imag(K{1}))) ,numKgrid ).'; 
end

% To hold the spectrum
S.Nk = length(kgrid)-1;
S.kgrid = (kgrid(1:end-1) +  kgrid(2:end))/2;

% Count the number of modes in each bin to normalize later -- this gives a
% smoother result, since we want the average energy in each bin.
oneG = ones(size(KX));
S.nbin = spect1D(oneG,oneG,Kmag,kgrid)*numel(oneG)^2;
S.nnorm = S.nbin./S.kgrid.^2; % k^2 accounts for the fact that in 3D, number of modes in shell increases with k^2
S.nnorm = S.nnorm/mean(S.nnorm); % Normalization by number of modes

S.EK = 0;
for var = {'vel1','vel2','vel3'}
    ft = fftn(D.(var{1}));
    S.(var{1}) = spect1D(ft,ft,Kmag,kgrid);
    S.EK = S.EK + S.(var{1}); % Total spectrum is the sum of each component
end
S.EM = 0;
if MHD
    for var = {'Bcc1','Bcc2','Bcc3'}
        ft = fftn(D.(var{1}));
        S.(var{1}) = spect1D(ft,ft,Kmag,kgrid);
        S.EM = S.EM + S.(var{1});
    end
end


%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%
loglog(S.kgrid, S.vel1, S.kgrid, S.kgrid.^(-5/3),'k:' )
ylabel('$E_K$','interpreter','latex')
xlabel('$k$','interpreter','latex')
if MHD
    hold on 
    loglog(S.kgrid, S.EM)
end


end


function out = spect1D(v1,v2,K,kgrid)
% Function to find the spectrum <v1 v2>, 
% K is the kgrid associated with v1 and v2
% kgrid is the grid for spectral shell binning

nk = length(kgrid)-1;
out = zeros(nk,1);
NT2 = numel(K)^2;
for kk = 1:nk
    out(kk) = sum( real(v1(K<kgrid(kk+1) & K>kgrid(kk)).*conj(v2(K<kgrid(kk+1) & K>kgrid(kk)))) )/NT2;
end

end

