function test_athena_CGL_waves


% Function to test the evolution of linear waves in Athena, and how much
% they damp compared to the expected amount. Also check convergence with
% resolution.
% Also adapt to measure Braginskii damping 

% Useful bits
P.tests = '~/Research/athena/wave-tests/';
P.home = '~/Research/athena/athena/vis/matlab/';
P.vtkout = '~/Research/athena/wave-tests/cgl/';
P.save_name = '~/Research/athena/wave-tests/waveTestData/threeDoblique_IdealCGLandISO.mat';
runFullSimulations = 1;

% Parameters.
P.pprp = 8;
P.pprl = 5;
P.Bx = 1.0;
P.By = sqrt(2);
P.Bz = 0.5;
P.rho = 1;
P.vx = 0.0;

if runFullSimulations 
% load(P.save_name,'P','T'); % So you don't immediately overwrite by mistake
% Variables to loop over
rm = @(f,n) repmat(f,[1,n]);

% % From oneD_IdealCGLandISO
% P.typesV = [rm({'hlle'},4*8) rm({'hlld'},4*8) rm({'iso'},3*8)];
% P.nxV = [rm([8,16,32,64,128,256,512,1024],4),...
%     rm([8,16,32,64,128,256,512,1024],4),...
%     rm([8,16,32,64,128,256,512,1024],3)];
% P.wavesV = [rm(5,8) rm(6,8) rm(7,8) rm(4,8) ... % Alf, slow, fast, weird cgl wave
%     rm(5,8) rm(6,8) rm(7,8) rm(4,8) ...
%     rm(4,8) rm(3,8) rm(5,8)];
% P.varV = [rm({'mom2'},3*8) rm({'dens'},1*8) ...
%     rm({'mom2'},3*8) rm({'dens'},1*8) ...
%     rm({'mom2'},3*8)];

% From threeDoblique_IdealCGLandISO
P.typesV = [rm({'hlle'},3*4) rm({'hlld'},3*4) rm({'iso'},3*4)];
P.nxV = [rm([8,16,32,64],3),...
    rm([8,16,32,64],3),...
    rm([8,16,32,64],3)];
P.wavesV = [rm(5,4) rm(6,4) rm(7,4) ... % Alf, slow, fast
    rm(5,4) rm(6,4) rm(7,4)  ...
    rm(4,4) rm(3,4) rm(5,4)];

% % From viscousTest_ISO
% P.typesV = [rm({'iso'},3*8) ];
% P.wavesV = [rm(5,8) rm(3,8) rm(5,8)]; % Alf, slow, fast
% P.visV = [rm([0.001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3],3)];

P.nsims = length( P.typesV );
for simn = 32%1:P.nsims


    % Simulation parameters
    P.waveflag = P.wavesV(simn);
    P.nu_iso = 0.;%P.visV(simn);
    P.dtout = 0.1;
    P.tf = 2;
    P.n = [P.nxV(simn) P.nxV(simn) P.nxV(simn)];
    P.l = [1 1 1 ];
    P.oblique = 1;
    P.type = P.typesV{simn};
    % Testing parameters
    P.yzind = [P.nxV(simn)/2 P.nxV(simn)/2]; % y and z ind of the x line to compute error on
    P.var = 'mom2'; % Variable to compare, number is position in matrix

    % Run the simulation
    writeInputFile(P);
    T{simn}.sim_time = runAthena(P); % T stores the test outputs


    % Collect the results, and compute errors 
    filename = @(n) [P.vtkout '/LinWave.block0.out2.' sprintf('%05d',n) '.vtk'];
    nums = 0:round(P.tf/P.dtout);
    
    % Choose points to compare and extract grid
    V = readVTKpp(filename(0));
    [X,Y,Z] = ndgrid(V.x,V.y,V.z);
    [kwave, Xk]=workOutK(P,X,Y,Z);
    % Work out expected wave
    if strcmp(P.type,'iso')
        [R, vals, L] = formIsoMHDMatrix(kwave,sqrt(P.pprl),P.Bx,P.By,P.Bz,P.rho,P.vx,P.nu_iso);
    else
        [R, vals, L] = formCGLMatrix(P.pprp,P.pprl,P.Bx,P.By,P.Bz,P.rho,P.vx);
    end
    freq = vals(P.waveflag+1);
    evec = R(:,P.waveflag+1);
    omega = freq*kwave;
    % Expected form of wave
    athsol = makeSol2Compare(V,P);
    x = double(V.x);
    [phase0,amp0] = findPhase(x,P.l(1),athsol); % Fit sin to find phase at t=0
    wave_expect = @(t) amp0*sin(2*pi/P.l(1)*V.x + phase0 - real(omega)*t)*exp(imag(omega)*t);
    T{simn}.errs = zeros(1,length(nums));
    T{simn}.t = zeros(1,length(nums));
    T{simn}.amps = zeros(9,length(nums));
    T{simn}.energy = zeros(4,length(nums));
    for nnn=nums
        V = readVTKpp(filename(nnn));
        athsol =  makeSol2Compare(V,P);
        T{simn}.errs(nnn+1) = sqrt( sum( (athsol - wave_expect(V.t)).^2) )./ ...
            sqrt( sum( (0.5*athsol+ 0.5*wave_expect(V.t)).^2) );
        T{simn}.t(nnn+1) = V.t;
        T{simn}.amps(:,nnn+1) = perturbationAmplitudes(V,P);
        T{simn}.energy(:,nnn+1) = perturbationEnergy(V,P);
        plot(V.x,athsol,V.x,wave_expect(V.t))
        title(['t=' num2str(V.t)])
        drawnow;ginput(1);
    end
    T{simn}.decay_rate = mesureDecayRate(T{simn}.t,T{simn}.amps,P);
    T{simn}.decay_all = T{simn}.decay_rate(~isnan(T{simn}.decay_rate));
    T{simn}.decay_all = mean(T{simn}.decay_all);
    T{simn}.expected_decay = imag(omega);
    T{simn}.expected_freq = real(omega);
%     save(P.save_name,'P','T')
end
    
else 
    load(P.save_name,'P','T')
end


% Plotting for oneD_IdealCGLandISO

% Plot error at a givevn time as a function of resolution, for each of the
% different waves
nt = length(T{1}.t);
nr = 8;
error_hlld = zeros(4,nr,nt);
timing_hlld = zeros(4,nr);
error_hlle = zeros(4,nr,nt);
timing_hlle = zeros(4,nr);
error_iso = zeros(3,nr,nt);
timing_iso = zeros(3,nr);
hlld_wn = [0 0 0 0 4 1 2 3]; % hlld_wave_num(wave)=index for plotting
hlle_wn = [0 0 0 0 4 1 2 3]; % hlld_wave_num(wave)=index for plotting
iso_wn = [0 0 0 2 1 3]; % hlld_wave_num(wave)=index for plotting
nxind = @(nx) log2(nx)-2;
for simn = 1:P.nsims
    if strcmp(P.typesV{simn},'hlld')
        error_hlld(hlld_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = reshape(T{simn}.errs,[1,1,nt]);
        timing_hlld(hlld_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = T{simn}.sim_time;
    elseif strcmp(P.typesV{simn},'hlle')
        error_hlle(hlle_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = reshape(T{simn}.errs,[1,1,nt]);
        timing_hlle(hlle_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = T{simn}.sim_time;
    elseif strcmp(P.typesV{simn},'iso')
        error_iso(iso_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = reshape(T{simn}.errs,[1,1,nt]);
        timing_iso(iso_wn(P.wavesV(simn)+1),nxind(P.nxV(simn)),:) = T{simn}.sim_time;
    end
end

resV = [8,16,32,64,128,256,512,1024];
wavestyle = {'-','--','-.',':'};
tval = nt;
for www = 1:4
    loglog(resV, error_hlld(www,:,tval),wavestyle{www},'Color','b')
%     loglog(nx, timing_hlld(www,:),wavestyle{www},'Color','b')
    hold on
end
for www = 1:4
    loglog(resV, error_hlle(www,:,tval),wavestyle{www},'Color','g')
%     loglog(nx, timing_hlle(www,:),wavestyle{www},'Color','g')
end
for www = 1:3
    loglog(resV, error_iso(www,:,tval),wavestyle{www},'Color','r')
%     loglog(nx, timing_iso(www,:),wavestyle{www},'Color','r')
end
loglog(resV,80*resV.^(-2),'k:')
ylim([2e-4 3])
ylabel('Error')
xlabel('N_x')
legend({'Alfven','Slow','Fast','Entropy type 2'})


end

function out = makeSol2Compare(F,P)
if P.n(2)==1
    out = double(F.(P.var) -  mean(F.(P.var)) ).';
elseif P.n(3)==1
    out = double(F.(P.var)(:,P.yzind(1)) - ...
        mean(F.(P.var)(:,P.yzind(1))));
else
    out = double(F.(P.var)(:,P.yzind(1),P.yzind(2)) - ...
        mean(F.(P.var)(:,P.yzind(1),P.yzind(2))));
end
end

function amps = perturbationAmplitudes(V,P)
if P.n(2)==1
    dv = V.x(2)-V.x(1);
elseif P.n(3)==1
    dv = (V.x(2)-V.x(1))*(V.y(2)-V.y(1));
else
    dv = (V.x(2)-V.x(1))*(V.y(2)-V.y(1))*(V.z(2)-V.z(1));
end
sm = @(f) sum(sum(sum(f)));
if strcmp(P.type,'iso')
    % Measure amplitude of perturbations in each variable from rms
    amps = [sm( sqrt( (V.dens-mean(mean(mean((V.dens))))).^2) );...
        sm( sqrt( (V.mom1-mean(mean(mean((V.mom1))))).^2) );...
        sm( sqrt( (V.mom2-mean(mean(mean((V.mom2))))).^2) );...
        sm( sqrt( (V.mom3-mean(mean(mean((V.mom3))))).^2) );...
        0;0;...
        sm( sqrt( (V.Bcc1-mean(mean(mean((V.Bcc1))))).^2) );...
        sm( sqrt( (V.Bcc2-mean(mean(mean((V.Bcc2))))).^2) );...
        sm( sqrt( (V.Bcc3-mean(mean(mean((V.Bcc3))))).^2) )] / dv;
else 
    % Measure amplitude of perturbations in each variable from rms
    amps = [sm( sqrt( (V.dens-mean(mean(mean((V.dens))))).^2) );...
        sm( sqrt( (V.mom1-mean(mean(mean((V.mom1))))).^2) );...
        sm( sqrt( (V.mom2-mean(mean(mean((V.mom2))))).^2) );...
        sm( sqrt( (V.mom3-mean(mean(mean((V.mom3))))).^2) );...
        sm( sqrt( (V.Etot-mean(mean(mean((V.Etot))))).^2) );...
        sm( sqrt( (V.mu-mean(mean(mean((V.mu))))).^2) );...
        sm( sqrt( (V.Bcc1-mean(mean(mean((V.Bcc1))))).^2) );...
        sm( sqrt( (V.Bcc2-mean(mean(mean((V.Bcc2))))).^2) );...
        sm( sqrt( (V.Bcc3-mean(mean(mean((V.Bcc3))))).^2) )] / dv;
end
end

function en = perturbationEnergy(V,P)
if P.n(2)==1
    dv = V.x(2)-V.x(1);
elseif P.n(3)==1
    dv = (V.x(2)-V.x(1))*(V.y(2)-V.y(1));
else
    dv = (V.x(2)-V.x(1))*(V.y(2)-V.y(1))*(V.z(2)-V.z(1));
end
sm = @(f) sum(sum(sum(f)));
if strcmp(P.type,'iso')
    % Measure amplitude of perturbations in each variable from rms
    en = [sm( 0.5*(V.mom1.^2 + V.mom1.^2 + V.mom1.^2)./V.dens );...
        sm( 0.5*(V.Bcc1.^2 + V.Bcc2.^2 + V.Bcc3.^2) ); ...
        sm( V.dens ); 0 ] / dv;
else 
    % Measure amplitude of perturbations in each variable from rms
    en = [sm( 0.5*(V.mom1.^2 + V.mom1.^2 + V.mom1.^2)./V.dens );...
        sm( 0.5*(V.Bcc1.^2 + V.Bcc2.^2 + V.Bcc3.^2) ); ...
        sm( V.Etot );sm( V.dens )] / dv;
end
end

function decay = mesureDecayRate(t,amps,P)
% Measure decay rate of amps 
lamp = size(amps,1);
decay=NaN*zeros(lamp,1);
for nnn=1:lamp
    if amps(nnn,1)>1e-8
        fit = polyfit(t,log(amps(nnn,:)),1);
        decay(nnn)= fit(1);
    end
end

end

function [k, Xk]=workOutK(P, X,Y,Z)
% Work out the k that athena uses
if P.oblique==0 ||  P.n(3)==1 % 1D or 2D assume parallel to axis
    k=2*pi/P.l(1);
    Xk = X;
else
    % Just copied from athena here
    x1size = P.l(1);x2size = P.l(2);x3size = P.l(3);
    ang_3 = atan(x1size/x2size);
    sin_a3 = sin(ang_3);
    cos_a3 = cos(ang_3);
    ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
    sin_a2 = sin(ang_2);
    cos_a2 = cos(ang_2);
    x1 = x1size*cos_a2*cos_a3;
    x2 = x2size*cos_a2*sin_a3;
    x3 = x3size*sin_a2;
    lambda = x1;
    lambda = min(lambda,x2);
    lambda = min(lambda,x3);

    k = 2*pi/lambda;
    Xk = cos_a2*(X*cos_a3 + Y*sin_a3) + Z*sin_a2;
end

end

function [phase,amp] = findPhase(x,lx,f)
% Fit phase of initial sin wave
amp = (max(f)-min(f))/2;
f = f/amp;
f = f-mean(f);
[mf, pos] = max(f); pos = pos/length(x);
fun = @(p) sin(2*pi*lx*x + p) - f;
phase = lsqnonlin(fun,2*pi*(pos-0.25),0,2*pi);

disp('')
disp(['Initial wave found with amp = ' num2str(amp) ', phase = ' num2str(phase)])
disp('')
end

function [R, vals, L] = formCGLMatrix(pprp,pprl,Bx,By,Bz,rho,vx)
% Forms A matrix 
vy=0;vz=0;
% All the conserved variables, and other useful quantities
mx = vx/rho;my = vy/rho;mz = vz/rho;
m2 = mx^2 + my^2 + mz^2;
dp = pprp - pprl; B2 = Bx^2 + By^2 + Bz^2;
B = sqrt(B2);
bhx = Bx/B;bhy = By/B;bhz = Bz/B;
mu = pprp/B;
E = pprp + pprl/2 + m2/rho/2 + B2/2;
Bdm = Bx*mx + By*my + Bz*mz;
% Useful bits to make matrix seem simpler
% Derivatives of 1+dp/B2
fh = 1+dp/B2;
fhrho = -m2/rho^2/B2;
fhmx = 2*vx/B2;
fhmy = 2*vy/B2;
fhmz = 2*vz/B2;
fhE = -2/B2;
fhmu = 3/B;
fhBy = (3*mu*bhy + 2*By)/B2 - 2*By*dp/B2^2;
fhBz = (3*mu*bhz + 2*Bz)/B2 - 2*Bz*dp/B2^2;

A1 = [0,      1,       0,      0,      0,      0,      0,      0 ];...
A2 = [-vx^2-Bx^2*fhrho,  2*vx-Bx^2*fhmx,  -Bx^2*fhmy,  -Bx^2*fhmz,  -Bx^2*fhE, ...
            B-Bx^2*fhmu,   By+mu*bhy-Bx^2*fhBy,   Bz+mu*bhz-Bx^2*fhBz];
A3 = [-vx*vy-Bx*By*fhrho,   vy-Bx*By*fhmx,   vx-Bx*By*fhmy,  -Bx*By*fhmz,...
            -Bx*By*fhE,  -Bx*By*fhmu,   -Bx*fh-Bx*By*fhBy,   -Bx*By*fhBz];
A4 = [-vx*vz-Bx*Bz*fhrho,   vz-Bx*Bz*fhmx,   -Bx*Bz*fhmy,  vx-Bx*Bz*fhmz,...
            -Bx*Bz*fhE,  -Bx*Bz*fhmu,   -Bx*Bz*fhBy,   -Bx*fh-Bx*Bz*fhBz];
A5 = [-1/rho^2*(E*mx+mx*B2/2+mu*B*mx-Bx*Bdm*fh)-Bx*Bdm/rho*fhrho,...
            E/rho+B2/rho/2+mu*B/rho-Bx^2/rho*fh-Bx*Bdm/rho*fhmx,...
            -Bx*By/rho*fh-Bx*Bdm/rho*fhmy,-Bx*Bz/rho*fh-Bx*Bdm/rho*fhmz,...
            vx-Bx*Bdm/rho*fhE, B*vx-Bx*Bdm/rho*fhmu,...
            By*vx+bhy*vx*mu-vy*Bx*fh-Bx*Bdm/rho*fhBy,...
            Bz*vx+bhz*vx*mu-vz*Bx*fh-Bx*Bdm/rho*fhBz];
A6 = [-mu*mx/rho^2,   mu/rho,   0,   0,   0,   vx,   0,   0];
A7 = [-(mx*By-my*Bx)/rho^2,   By/rho,   -Bx/rho,   0,   0,   0,   vx,   0];
A8 = [-(mx*Bz-mz*Bx)/rho^2,   Bz/rho,   0,   -Bx/rho,   0,   0,   0,   vx];


A = [A1;A2;A3;A4;A5;A6;A7;A8];

[R, vals, L] = eig(A);
[~,ord] = sort(real(diag(vals)));
vals = diag(vals(ord,ord));
R = R(:,ord);
L=L(ord,:);
end

function [R, vals, L] = formIsoMHDMatrix(k,cs,Bx,By,Bz,rho,vx,nu)
X=0;Y=1; % Parameters from Stone
% Viscosity in athena seems to be div(rho*nu*(grad(u)+div(u)))
A1 = [0,      1,       0,      0,      0,      0,      ];
A2 = [vx^2+cs^2+X+2i*nu*k*vx,   2*vx-2i*nu*k,   0,   0,    By*Y,    Bz*Y];
A3 = [1i*nu*k*vx,     -1i*nu*k,     vx-1i*nu*k,     0,     -Bx,     0];
A4 = [1i*nu*k*vx,     -1i*nu*k,     0,     vx-1i*nu*k,     0,     -Bx];
A5 = [-By*vx/rho,    By/rho,   -Bx/rho,    0,   vx,    0];
A6 = [-Bz*vx/rho,    Bz/rho,   0,    -Bx/rho,   0,    vx];

A = [A1;A2;A3;A4;A5;A6];

[R, vals, L] = eig(A);
[~,ord] = sort(real(diag(vals)));
vals = diag(vals(ord,ord));
R = R(:,ord);
L=L(ord,:);
end


function tout = runAthena(P)

cd(P.tests)
if strcmp(P.type,'hlld')
    tic
    unix(['rm cgl/*;./athena-hlld -i athinput.linear_wave -d cgl']);
    tout = toc;
elseif strcmp(P.type,'hlle')
    tic
    unix(['rm cgl/*;./athena-hlle -i athinput.linear_wave -d cgl']);  
    tout = toc;
elseif strcmp(P.type,'iso')
    tic
    unix(['rm cgl/*;./athena-iso -i athinput.linear_wave -d cgl']);  
    tout = toc;
end
cd(P.home)

end


function writeInputFile(P)

%More loops could easily be added for various variables, care would need to
%be taken if looping over the viscosity coefficients

cd(P.tests)

fid = fopen('athinput.linear_wave','w');

fprintf(fid,'<comment>\n');
fprintf(fid,'problem   = linear wave\n');
fprintf(fid,'reference =\n');
fprintf(fid,'configure = --prob=linear_wave\n');
fprintf(fid,'\n');
fprintf(fid,'<job>\n');
fprintf(fid,'problem_id = LinWave  # problem ID: basename of output filenames\n');
fprintf(fid,'\n');
fprintf(fid,'<output1>\n');
fprintf(fid,'file_type  = hst      # History data dump\n');
fprintf(fid,'dt         = %1.2e    # time increment between outputs\n',P.dtout);
fprintf(fid,'\n');
fprintf(fid,'<output2>\n');
fprintf(fid,'file_type  = vtk      # VTK data dump\n');
fprintf(fid,'variable   = cons     # variables to be output\n');
fprintf(fid,'dt         = %1.2e    # time increment between outputs\n',P.dtout);
fprintf(fid,'\n');
fprintf(fid,'<time>\n');
fprintf(fid,'cfl_number = 0.3      # The Courant, Friedrichs, & Lewy (CFL) Number\n');
fprintf(fid,'nlim       = -1       # cycle limit\n');
fprintf(fid,'tlim       = %g       # time limit\n',P.tf);
fprintf(fid,'integrator = vl2      # time integration algorithm\n');
fprintf(fid,'xorder     = 2        # order of spatial reconstruction\n');
fprintf(fid,'ncycle_out = 1000000000  # interval for stdout summary info\n');
fprintf(fid,'\n');
fprintf(fid,'<mesh>\n');
fprintf(fid,'nx1        = %i       # Number of zones in X1-direction\n',P.n(1));
fprintf(fid,'x1min      = 0.0      # minimum value of X1\n');
fprintf(fid,'x1max      = %g       # maximum value of X1\n', P.l(1)); %(Change back to 2*xmaxfor 3d)
fprintf(fid,'ix1_bc     = periodic # inner-X1 boundary flag\n');
fprintf(fid,'ox1_bc     = periodic # outer-X1 boundary flag\n');
fprintf(fid,'\n');
fprintf(fid,'nx2        = %i       # Number of zones in X2-direction\n',P.n(2));
fprintf(fid,'x2min      = 0.0      # minimum value of X2\n');
fprintf(fid,'x2max      = %g       # maximum value of X2\n',P.l(2));
fprintf(fid,'ix2_bc     = periodic # inner-X2 boundary \n');
fprintf(fid,'ox2_bc     = periodic # outer-X2 boundary flag\n');
fprintf(fid,'\n');
fprintf(fid,'nx3        = %i       # Number of zones in X3-direction\n',P.n(3));
fprintf(fid,'x3min      = 0.0      # minimum value of X3\n');
fprintf(fid,'x3max      = %g       # maximum value of X3\n',P.l(3));
fprintf(fid,'ix3_bc     = periodic # inner-X3 boundary flag\n');
fprintf(fid,'ox3_bc     = periodic # outer-X3 boundary flag\n');
fprintf(fid,'\n');
fprintf(fid,'num_threads = 1       # maximum number of OMP threads\n');
fprintf(fid,'refinement  = none\n');
fprintf(fid,'\n');
fprintf(fid,'<meshblock>\n');
fprintf(fid,'nx1        = %i       # Number of zones in X1-direction\n',P.n(1));
fprintf(fid,'nx2        = %i       # Number of zones in X2-direction\n',P.n(2));
fprintf(fid,'nx3        = %i       # Number of zones in X3-direction\n',P.n(3));
fprintf(fid,'\n');
fprintf(fid,'<hydro>\n');
fprintf(fid,'gamma = 1.66666666667 # gamma = C_p/C_v\n');
fprintf(fid,'iso_sound_speed = %g      # isothermal sound speed\n', sqrt(P.pprl));
fprintf(fid,'\n');
fprintf(fid,'<problem>\n');
fprintf(fid,'compute_error = false # set value to "true" to compute L1 error compared to initial data\n');
fprintf(fid,'wave_flag  = %i       # Wave Type Identifier\n', P.waveflag);
fprintf(fid,'amp        = 1e-6     # Wave Amplitude\n');
fprintf(fid,'vflow      = %g      # Background Flow Velocity\n',P.vx);
if P.oblique==0 || P.n(3)==1; ang=0;else;ang=-999.9;end
fprintf(fid,'ang_2      = %g      # Angle from x-axis in xy-plane\n',ang);
fprintf(fid,'ang_3      = %g      # Angle from x-axis in xz-plane\n',ang);
%fprintf(fid,'ang_2_vert = false    # Set true to point wavevector along y-axis\n');
%fprintf(fid,'ang_3_vert = false    # Set true to point wavevector along z-axis\n');
fprintf(fid,'nu_iso     = %g       # Isotropic Viscosity Coefficient\n', P.nu_iso);
%fprintf(fid,'nu_aniso   = %g       # Anisotropic Viscosity Coefficient\n', nu_aniso);
fprintf(fid,'kn         = 1        # Number of wavelengths in box\n');
fprintf(fid,'bx0        = %1.16f      # Magnetic Field\n',P.Bx);
fprintf(fid,'by0        = %1.16f       \n',P.By);
fprintf(fid,'bz0        = %1.16f       \n',P.Bz);
fprintf(fid,'pp0      = %1.16f      \n',P.pprp);
fprintf(fid,'pl0      = %1.16f      \n',P.pprl);
fclose(fid);


cd(P.home)
end


