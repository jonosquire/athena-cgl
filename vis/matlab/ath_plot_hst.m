function ath_plot_hst%( fname )
folder = '~/Research/athena/turb-tests/decay'; % Folder with outputs
file = 'Turb'; % Name of output
% Plots variables from hst file 
fulldata = importdata([folder '/' file '.hst']);
names = strsplit(fulldata.textdata{2},'  ');
names = names(1:end);


% nums = input(['Choose variables:  ' char(10) strjoin(names,char(10)) char(10)]);
hold on
dat = fulldata.data;
inds = find_restart_ind(dat(:,1)); % Remove intermediate times between a restart
t = dat(inds,1);

kinetic_energy = dat(inds,7)+dat(inds,8)+dat(inds,9);
semilogy(t,kinetic_energy,'')
hold on 
legend(names(nums))


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
