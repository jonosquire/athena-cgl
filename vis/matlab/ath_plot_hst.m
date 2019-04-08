function ath_plot_hst%( fname )
folder = '~/Research/athena/turb-tests/output'; % Folder with outputs
file = 'Turb324_vl_knl'; % Name of output
% Plots variables from hst file 
fulldata = importdata([folder '/' file '.hst']);
names = strsplit(fulldata.textdata{2},'  ');
names = names(1:end);


% nums = input(['Choose variables:  ' char(10) strjoin(names,char(10)) char(10)]);
dat = fulldata.data;
inds = find_restart_ind(dat(:,1)); % Remove intermediate times between a restart
t = dat(inds,1);
vol = 16;

figure;
if size(dat,2)==14;m1=13;m2=14;else;m1=12;m2=13;end
kinetic_energy = (dat(inds,7)+dat(inds,8)+dat(inds,9))/vol;
magnetic_energy = (dat(inds,m1)+dat(inds,m2))/vol;
uy = sqrt(2*(dat(inds,8))/vol);uz = sqrt(2*(dat(inds,9))/vol);
by = sqrt(2*(dat(inds,m1))/vol);bz = sqrt(2*(dat(inds,m2))/vol);
plot(t,uy,t,uz,t,by,t,bz)
norm = 2*dat(1,m1-1)/vol;
semilogy(t,kinetic_energy/norm,t,magnetic_energy/norm, ...
    t,(kinetic_energy+magnetic_energy)/norm,'k')
% legend(names(nums))


end


function inds = find_restart_ind(t)
i2d = find(diff(t)<0);
if ~isempty(i2d)
    if length(i2d)>1
        error('Upgrade to deal with multiple restarts')
    end
    tAboveRS = t(i2d+1);
    i2cut = find(t>tAboveRS,1);
    inds = 1:length(t);
    inds(inds>=i2cut & inds<=i2d)=[];
else 
    inds = 1:length(t);
end

end
