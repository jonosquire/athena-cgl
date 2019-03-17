function testOutAnalyticFWithMirrorHoles()

% Try out different analytic distribution functions that lack trapped
% particles, to then calculate the stability

vprp = linspace(0,1.5,101);
vprl = linspace(-1.5,1.5,101);

[Vprl,Vprp]=ndgrid(vprl,vprp);

maxwell = Vprp.*exp(-Vprp.^2-Vprl.^2);

imagesc(vprp,vprl,maxwell.')
set(gca,'YDir','normal')