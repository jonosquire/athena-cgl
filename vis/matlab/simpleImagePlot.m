function simpleImagePlot()


folder = '~/Research/athena/wave-tests/cgl';

file = 'LinWave'; % Name of output
output_id = 2; % Output id (set in input file)
nums = 0:100;

filename = @(n,oid) [folder '/' file '.block0.out' num2str(oid) '.'  sprintf('%05d',n) '.vtk'];


lims = 0.5*[-1 1];
dpstore = [];ts = [];
for nnn = nums
    V = readVTKpp(filename(nnn,output_id));
    ar = length(V.x)/length(V.y);
    
%     imagesc(V.x,V.y,V.Bcc2(:,:,1).',lims)
%     pbaspect([ar 1 1])
%     title(['t = ' num2str(V.t)]);    colorbar

    oneD = @(a) a(:,8,8);%@(a) mean(mean(a,2),3);
    subplot(211)
    plot(V.x,oneD(V.Bcc2),V.x,oneD(V.vel2))
    hold on;
    ax = gca;ax.ColorOrderIndex = 1;
    plot(V.x,oneD(V.Bcc3),'--',V.x,oneD(V.vel3),'--')
    hold off
    ylim(lims)
    title(['t = ' num2str(V.t)])
    subplot(212)
    plot(V.x,oneD((V.pprp-V.pprl)./(V.Bcc1.^2+ V.Bcc2.^2+V.Bcc3.^2)))
    ylim([-2 1])
    
    dpstore = [dpstore mean(mean(mean(V.pprp-V.pprl)))];
    ts = [ts V.t];
    
    drawnow
    ginput(1);
end
% semilogy(ts,dpstore,'r',ts,0.3*exp(-0.1*ts),':k','Linewidth',3)