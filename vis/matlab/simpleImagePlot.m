function simpleImagePlot()


folder = '~/Research/athena/wave-tests/cgl';
% folder = '~/Research/athena/turb-tests/decay';

file = 'LinWave'; % Name of output
output_id = 2; % Output id (set in input file)
nums = 0:250;

filename = @(n,oid) [folder '/' file '.block0.out' num2str(oid) '.'  sprintf('%05d',n) '.vtk'];readfunc = @(file) readVTKpp(file);
filename = @(n,oid) [folder '/' file '.out' num2str(oid) '.'  sprintf('%05d',n) '.athdf'];readfunc = @(file) readHDF5(file);


dpstore = [];estore=[];ts = [];
lims = 1*[-1 1];
for nnn = nums
    V = readfunc(filename(nnn,output_id));
    ar = length(V.x)/length(V.y);
    
%     tp = @(f) f.';%-mean(mean(mean(f)));
%     vars = {'vel2'};
%     for kkk=1:length(vars)
%         subplot(1,length(vars),kkk)
%         imagesc(V.x,V.y,tp(V.(vars{kkk})(:,:,1)))
%         pbaspect([ar 1 1])
%         title(['t = ' num2str(V.t)]);    colorbar
%     end
    
    if length(V.y)>1
        oneD = @(a) a(:,6,1);%@(a) mean(mean(a,2),3);
    else 
        oneD = @(a) a;%@(a) mean(mean(a,2),3);
    end
    subplot(211)
    plot(V.x,oneD(V.vel1),V.x,oneD(V.vel2))
    hold on;
    ax = gca;ax.ColorOrderIndex = 1;
    plot(V.x,oneD(V.Bcc2),'--',V.x,oneD(V.Bcc3),'--')
    hold off
    ylim(lims)
    title(['t = ' num2str(V.t)])
    subplot(212)
    if isfield(V,'pprp')
        plot(V.x,oneD((V.pprp-V.pprl)./(V.Bcc1.^2+ V.Bcc2.^2+V.Bcc3.^2)),'-')
    end
    ylim([-2 1])
    
%     dpstore = [dpstore mean(mean(mean(V.pprp-V.pprl)./(V.Bcc1.^2+ V.Bcc2.^2+V.Bcc3.^2)))];
%     ts = [ts V.t];
%     estore =[estore mean(mean(mean(0.5*V.Bcc2.^2 + 0.5*V.rho.*V.vel2.^2)))];
%     
    drawnow
    pause(0.0)
%      ginput(1);
end
% semilogy(ts,-dpstore,'r',ts,0.5*exp(-1*ts),':k','Linewidth',3)
subplot(212)
hold on
omA=2*pi;tp0onuc=3/10;
semilogy(ts, estore, ts, estore(1)./(1+estore(1)*2*tp0onuc*omA^2*ts))


end