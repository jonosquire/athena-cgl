% Plot_ath_hst.m

% General script for plotting of lots of ath_plot_hst simulations

folder = '../../../turb-tests/stampede-decay/';
files = {'MHD-N32-b128', 'visbrag-N32-b128', 'cglmhd-N32-b128-hlle',...
    'MHD-N64-b128', 'visbrag-N64-b128', 'cglmhd-N64-b128-hlle',...
     'visbrag-N128-b128', 'cglmhd-N128-b128-hlle'};
legnames = {'MHD $\beta = 128,\; N_{\perp}=32$','Viscous $\beta = 128,\; N_{\perp}=32$','CGL MHD $\beta = 128,\; N_{\perp}=32$',...
    'MHD $\beta = 128,\; N_{\perp}=64$','Viscous $\beta = 128,\; N_{\perp}=64$','CGL MHD $\beta = 128,\; N_{\perp}=64$',...
    'Viscous $\beta = 128,\; N_{\perp}=128$','CGL MHD $\beta = 128,\; N_{\perp}=128$'};
linesty = {':',':',':',...
    '--','--','--',...
    '-','-'};
cgl = [0 0 1 ...
    0 0 1 ...
    0 1];
col = [1 2 3 ...
    1 2 3 ...
    2 3];


origCO=get(gcf,'DefaultAxesColorOrder');
origCO = [[0 0 0];origCO];
figure
for kk=1:8
    fname = [folder files{kk} '/Turb'];
    fulldata = importdata([fname '.hst']);
    dat = fulldata.data;
    m2 = size(dat,2);
    nums=[8 9 m2-1 m2];
    semilogy(dat(1:end,1),sum(dat(1:end,nums),2)/(2*dat(1,m2-2)),linesty{kk},'Color',origCO(col(kk),:),'Markersize',2)
    drawnow
    hold on
end

legend(legnames,'interpreter','latex') 
xlabel('$t$','interpreter','latex') 
ylabel('$E_{\mathrm{turb}}$','interpreter','latex') 
xlim([0 260])
ylim([5e-5 1e-2])