function analyzeCGLWaves()

% Analyze the frequency of CGL waves from athena
posx_maxmin=[0.25 0.75];
posy = 0.;posz=0.;
expected_freq = 3.3166247903554003;
k = 2*pi/1;
vars = {'vel1'};%,'vel2','vel3'
folder = '/Users/jsquire/Research/athena/wave-tests/cgl';
nums = 0:50;

filename = @(n) [folder '/LinWave.block0.out2.' sprintf('%05d',n) '.vtk'];

V = readVTKpp(filename(0));
% Find indicies
iy = find(V.y>posy,1);
iz = find(V.z>posz,1);
if isempty(iz);iz=1;end % For 2D
if isempty(iy);iy=1;end % For 1D
ix = [find(V.x>posx_maxmin(1),1); find(V.x>posx_maxmin(2),1)];
f = [];t=[];
for nnn=nums
    V = readVTKpp(filename(nnn));
    fl=[];
    for v=vars
        if isempty(V.y); field = V.(v{1}).';else; field = V.(v{1});end
        fl = [fl;field(ix,iy,iz)];
    end
    f = [f fl];
    t = [t V.t];
    drawnow
end
% f=f./repmat(f(:,1),[1,size(f,2)]);
expected = repmat(f(:,1),[1,length(t)]).*repmat(cos(k*expected_freq*t),[2,1]);
hold on
plot(t,f,t,expected,'--')
legend(vars)


end
