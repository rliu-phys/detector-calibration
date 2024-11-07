function f = findpowder(x,img_exp,plotflag)
% to load data
% ccd = double(h5read('h5/scan_5_028866.h5','/entry/data/data',[1,1,1],[1028,1062,1]));
% ccd = ccd.*(double(ccd<4e9));
% ccd = ccd.*(double(ccd>0));
% ccd = ccd';
% for ii = 2:20;ccd1 = double(h5read('h5/scan_5_028866.h5','/entry/data/data',[1,1,ii],[1028,1062,1]));ccd1 = ccd1.*(double(ccd1<4e9));ccd1 = ccd1.*(double(ccd1>0));ccd = ccd1'+ccd;end
% to view and fit
% findpowder(1e5*[1.8 -0.375 -0.30 .00037],ccd,1);
% x=fminsearch(@(x) findpowder(x,ccd,0),1e5*[1.8 -0.375 -0.30 .00037]);
% findpowder(x,ccd,1);
% x = [1.806360774795135e+05,-3.729935571548389e+04,-2.859329871988752e+04,36.970724978683790]

f1 = qbin_pilatus(img_exp, x(1), x(2), x(3), x(4), x(5), 0, 0, 0, plotflag);


if(plotflag)
    figure(701);imagesc(f1.powder);title('Thresholded experiment');axis image
    figure(702);imagesc(f1.theory);title('Theory');axis image
end
% two theta, R, (gam), X, Y
% input - R X Y 2th
if(plotflag)
    f=f1;
else
    f = f1.error;
end
end
