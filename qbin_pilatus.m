function dataout = qbin_pilatus(img_exp, Tthdet, Rdet, Gamdet, Xdet,Ydet,Samth, Samchi, Samphi, plotflag)

%---------------------------------------------
%PARAMETERS FOR MOMENTUM SPACE BINNING PIXIRAD IMAGES
% All vector notation in beam convention [+outboard +up +downstream]

% On/off CALCULATION FLAGS
Tthflag = 0; % output q axis in polar two theta (degrees), zero for 1/A
Binflag = 0; % do binning calculation
Powdflag = 1; % Calculate and display image with powder lines (DISABLED)


Ekev = 12.0; % beam energy in kev
kb = 2*pi*Ekev/12.39842; % beam momentum 1/A
%pxsz = 13; % pixel size in um (13.3mm x 1024 pixels)
hpxsz = 75; % pixel size in um EIGER
vpxsz = 75;
%hpxsz = 172; % pixel size in um PILATUS
%vpxsz = 172;
%hpxsz = 55; % pixel size in um (13.3mm x 1024 pixels)
%vpxsz = 55*sind(60); % pixel size in um (13.3mm x 1024 pixels)

%qmin = 0; % momentum space range start (zero for full ccd) 1/A
%qmax = 0; % momentum space range end (zero for full ccd)
qnum = 75; % number of momentum space bins

%gammin = 0; % gamma range start (zero for full ccd) rotation about beam axis
%gammax = 0; % gamma range end (zero for full ccd) (deg)
qlist=zeros(qnum,2);

Dzp = 150; % Zone plate diameter - um
drzp = 20; % Zone plate outermost zone - nm
Fzp = Dzp*drzp*Ekev/1.239842;  % Zone plate focal distance - um
Dstop = 60; % Zone plate central stop diameter - um

%lat = [5.6533 5.6533 5.6533];  %GaAs

lat = [5.4309 5.4309 5.4309]; % Si cubic lattice vectors a b c (A)
%lat = [3.6599 3.6599 3.6599]; % Si cubic lattice vectors a b c (A)
lat2 = [4.0782 4.0782 4.0782]; % Au cubic lattice vectors a b c (A)
refl2 = [0 0 4; 0 0 2; 0 2 2; 1 1 1; 1 1 3; 1 3 3; 3 3 3; 0 4 4]; % Au powder reflections to plot
%refl2 = [1 1 1;0 0 2];
refl = [0 0 4; 0 2 2; 1 1 1; 1 1 3; 1 3 3; 3 3 3; 0 4 4]; % Si powder reflections to plot
%refl =[];
%lat = [10.334 6.002  4.6915]; % cubic lattice vectors a b c (A)
%refl = [2 0 0; 0 2 0; 0 0 2; 2 2 0; 0 2 2; 2 2 0; 1 1 1; 1 0 1; 0 1 1; 1 1 0]; % which powder reflections to plot

%imdet = double(imread(imname));
imdet=img_exp;imname='Summed image';

if(nargin<8) 
    plotflag=1; % output images
end

if(nargin<7) 
    Samth =  24.9456; % Sample rotation about vertical axis (deg, right-hand-rule)
    Samchi = 0.0; % Sample rotation about beam axis (deg, RHR)
    Samphi = 0; % Sample rotation about outboard axis (deg, RHR)
end

if(nargin<4) 
% October experiment - run with '../../2017R3/20171017/Images/image_000011.tif'
%     Tthdet=35.2028; % **Cylindrical** two theta angle to the center of the detector (deg)
%     Rdet=171110;  % **Cylindrical** radius from sample to the center detector (um)
%     Gamdet=0;   % **Cylindrical** out of diffraction plane angle to the center of the detector (deg)
%     Xdet = 47364;
%     Ydet = -10046;
    Tthdet=35.2028; % **Cylindrical** two theta angle to the center of the detector (deg)
    Rdet=169510+5000;  % **Cylindrical** radius from sample to the center detector (um)
    Gamdet=0;   % **Cylindrical** out of diffraction plane angle to the center of the detector (deg)
    Xdet = 10500+4000;
    Ydet = -10046;
end

scrsz = get(0,'ScreenSize');  % Plot size parameters
imd = scrsz(3)/5;
imb = scrsz(4)-imd;

vNdet = size(imdet,1); 
vNcen = round(vNdet/2);
vimaxis = (1-vNcen):(vNdet-vNcen);
hNdet = size(imdet,2); 
hNcen = round(hNdet/2);
himaxis = (1-hNcen):(hNdet-hNcen);

imax = double(max(max(imdet)));
%Change scale to see weak peaks
%imtemp = imdet>(0.2*imax);
%imdet = imtemp*(0.2*imax)+(1-imtemp).*imdet;imax=0.2*imax;

if (plotflag==1)
    figure(101);set(101,'Position',[1 imb imd imd]);clf;imagesc(imdet);axis image;caxis([0 imax/100]);%colormap(gray);
    title(['Input image:' imname]);pause(0.5);
end

%initialize output images
imbin = imdet;

%Detector center and pixel plane vectors
tthd = pi * Tthdet / 180;
gamd = pi * Gamdet / 180;
Detcen = [Rdet*sin(tthd)*cos(gamd) Rdet*sin(gamd) Rdet*cos(tthd)*cos(gamd)];
Detunit = Detcen./(sqrt(dot(Detcen,Detcen)));
jvec = cross([0 1 0],Detunit); % horizontal pixels unit vector (low pixel number, low tth)
jvec = jvec./(sqrt(dot(jvec,jvec)));
ivec = cross(jvec,Detunit); % vertical pixels unit vector (low pixel number, more up)
%display(Detcen)
%display(ivec)
%display(jvec)
%Sample 00L rod unit vector
%samvec=[cos(sch)*cos(sth) sin(sch) (-cos(sch)*sin(sth))];
%Sample rotation from phi --> theta, [00L] (sqz) starts outboard at theta 0
sth=Samth*pi/180; sph=Samphi*pi/180; sch=Samchi*pi/180;
sROT1 = [cos(sth) 0 sin(sth); 0 1 0; -sin(sth) 0 cos(sth)];
sROT2 = [cos(sch) -sin(sch) 0; sin(sch) cos(sch) 0; 0 0 1];
sROT3 = [1 0 0; 0 cos(sph) -sin(sph); 0 sin(sph) cos(sph)];
sROT = sROT1*sROT2*sROT3; %nanoprobe stage stack
%sROT = sROT2*sROT1*sROT3; %LCLS stage stack
sqx = sROT*[0 0 1]'; sqy = sROT*[0 1 0]'; sqz = sROT*[1 0 0]';

kfmat=zeros(vNdet,hNdet,4);
kimat=zeros(vNdet,hNdet,4);
prefac=zeros(vNdet,hNdet);
qmat=zeros(vNdet,hNdet,4);

% Assign magnitude and vector q to active detector pixels within q and gamma range
%display('assigning q vectors...');
%Generate pixel arrays (um)
pix_x = repmat((himaxis.*(hpxsz)),vNdet,1)+Xdet;
pix_y = repmat((vimaxis.*(vpxsz))',1,hNdet)+Ydet;
%Calculate real space pixel positions in mm as a vector sum
kfmat(:,:,1) = pix_x(:,:).*jvec(1)+pix_y(:,:).*ivec(1)+Detcen(1);
kfmat(:,:,2) = pix_x(:,:).*jvec(2)+pix_y(:,:).*ivec(2)+Detcen(2);
kfmat(:,:,3) = pix_x(:,:).*jvec(3)+pix_y(:,:).*ivec(3)+Detcen(3);
% if(plotflag)
%        figure(11);clf;set(11,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,1));axis square;
%        title('X real space position'); pause(0.5);
%        figure(16);set(16,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,2));axis square;
%        title('Y real space position'); 
%        figure(17);set(17,'Position',[imd imb imd imd]);clf;imagesc(kfmat(:,:,3));axis square;
%        title('Z real space position');
% end
kfmat(:,:,4) = sqrt(kfmat(:,:,1).^2+kfmat(:,:,2).^2+kfmat(:,:,3).^2);
%Assign kfinal unit vector to each pixel based on position
kfmat(:,:,1) = kfmat(:,:,1)./kfmat(:,:,4);
kfmat(:,:,2) = kfmat(:,:,2)./kfmat(:,:,4);
kfmat(:,:,3) = kfmat(:,:,3)./kfmat(:,:,4);

%Find kinitial unit vector symmetric about 00L rod
prefac(:,:) = -2*(kfmat(:,:,1).*sqz(1)+kfmat(:,:,2).*sqz(2)+kfmat(:,:,3).*sqz(3));
kimat(:,:,1) = kfmat(:,:,1)+prefac(:,:).*sqz(1);
kimat(:,:,2) = kfmat(:,:,2)+prefac(:,:).*sqz(2);
kimat(:,:,3) = kfmat(:,:,3)+prefac(:,:).*sqz(3);
%Scale kinitial back to optic plane, find radial intersection in zp
prefac(:,:) = abs(Fzp./kimat(:,:,3)).*sqrt(kimat(:,:,1).*kimat(:,:,1)+kimat(:,:,2).*kimat(:,:,2));
kimat(:,:,4) = (prefac>(Dstop/2)).*(prefac<(Dzp/2));
if(Powdflag)
    kimat(:,:,1)=zeros(vNdet,hNdet);
    kimat(:,:,2)=zeros(vNdet,hNdet);
    kimat(:,:,3)=ones(vNdet,hNdet);
    kimat(:,:,4)=ones(vNdet,hNdet);
end
%Calculate q vectors
qmat(:,:,1) = kb*(kfmat(:,:,1)-kimat(:,:,1));
qmat(:,:,2) = kb*(kfmat(:,:,2)-kimat(:,:,2));
qmat(:,:,3) = kb*(kfmat(:,:,3)-kimat(:,:,3));
qmat(:,:,4) = sqrt(qmat(:,:,1).^2+qmat(:,:,2).^2+qmat(:,:,3).^2);
if(plotflag)
       figure(181);set(181,'Position',[imd imb imd imd]);clf;imagesc(acosd(kfmat(:,:,3)));axis image;
       title('Two theta per pixel map');colorbar;
       figure(19);set(19,'Position',[imd imb imd imd]);clf;imagesc(90-acosd(kfmat(:,:,2)));axis image;
       title('Gamma per pixel map');colorbar
       figure(20);set(19,'Position',[imd imb imd imd]);clf;imagesc(qmat(:,:,4));axis image;
       title('Momentum transfer (A^-1) per pixel map');colorbar
end
%Select active detector pixels in desired gamma range and q range
%dmask=(repmat(imaxis.^2,Ndet,1)+repmat((imaxis.^2)',1,Ndet))<Ncen^2; 
%if(qmax>0)
%    dmask=dmask.*(qmat(:,:,4)>=qmin);
%    dmask=dmask.*(qmat(:,:,4)<=qmax);
%end
%if(gammax>0)
%    dmask=dmask.*((atan(qmat(:,:,2)./qmat(:,:,1))*(180/pi))>=gammin);
%    dmask=dmask.*((atan(qmat(:,:,2)./qmat(:,:,1))*(180/pi))<=gammin);
  %  tempm=(atan(qmat(:,:,2)./qmat(:,:,1))*(180/pi))>=gammin;
  %  qmat(:,:,4)=qmat(:,:,4).*tempm;
  %  tempm=(atan(qmat(:,:,2)./qmat(:,:,1))*(180/pi))<=gammax;
  %  qmat(:,:,4)=qmat(:,:,4).*tempm;
%end
%qmat(:,:,4)=qmat(:,:,4).*dmask;

%imq=qmat(:,:,4); imqx=qmat(:,:,1);imqy=qmat(:,:,2);imqz=qmat(:,:,3);

if (plotflag==1)
%     figure(11);set(11,'Position',[imd imb imd imd]);clf;imagesc(imq);axis square;
%     title('Magnitude q per pixel'); pause(0.5);
%     figure(16);set(16,'Position',[imd imb imd imd]);clf;imagesc(imqx);axis square;
%     title('Magnitude qx per pixel'); 
%     figure(17);set(17,'Position',[imd imb imd imd]);clf;imagesc(imqy);axis square;
%     title('Magnitude qy per pixel'); 
%     figure(18);set(18,'Position',[imd imb imd imd]);clf;imagesc(imqz);axis square;
%     title('Magnitude qz per pixel'); 

end

% Create a list of average intensities binned in magnitude q
if(Binflag)
%    display('binning average intensity in magnitude q...');
    %qmat(:,:,4) = qmat(:,:,4).*kimat(:,:,4);
    tempm=qmat(:,:,4);
    qmin=min(tempm(find(tempm)));
    %qmin=round(1000*qmin)/1000;
    %qmax=qmin+0.035;
    qmax=max(tempm(find(tempm)));
    %qmax=round(1000*qmax)/1000;
    for i=1:qnum
        q1=qmin+(i-1)*(qmax-qmin)/qnum;
        q2=qmin+i*(qmax-qmin)/qnum;
        qlist(i,1)=q1;
        tempm1=qmat(:,:,4)>=q1;
        tempm2=qmat(:,:,4)<q2;
        tempm=tempm1.*tempm2;
        qlist(i,2)=(sum(sum(tempm))~=0)*(sum(sum(imdet.*tempm))/sum(sum(tempm)));
        tempm=(abs(qmat(:,:,4)-q1))<0.0001; 
        imbin=imbin.*(1-tempm)+imax*tempm;   %put bin lines on image
        if (Tthflag)
           qlist(i,1)=2*asin((12.39842/Ekev)*(q1+q2)/(2*4*pi))*180/pi; 
        end
    end
% Create unbinned complete list magnitude q and intensity for output
[r,c,v]=find(qmat(:,:,4));
numref=max(size(v));
outlist = zeros(numref,2);
outlist(:,1)=v(:);
for i=1:numref
    outlist(i,2)=imdet(r(i),c(i));
end

    if (plotflag==1)
        figure(12);set(12,'Position',[2*imd imb-imd imd 0.5*imd]);clf reset;plot(qlist(:,1),log(qlist(:,2)+1)); 
        ylabel('Avg Intensity (ADU)');
        if (Tthflag) xlabel('Two theta (deg)');
        else xlabel('Momentum transfer (1/A)');
        end
        figure(13);set(13,'Position',[2*imd imb imd imd]);clf;imagesc(imbin);axis square;
        title('Image with bin lines');pause(0.5);
    end
end

% Display powder lines
%%{
if(Powdflag)
    %display('generating powder lines...');
    numref=size(refl,1);
    iref=zeros(numref,1);
    jref=zeros(numref,1);
    numref2=size(refl2,1);
    iref2=zeros(numref2,1);
    jref2=zeros(numref2,1);
    tempm1=zeros(vNdet,hNdet);
    for ii=1:numref
        q1=2*pi*norm((refl(ii,:)./lat));
        tempm=(abs(qmat(:,:,4)-q1))<0.007;
        tempm1=tempm1+tempm;    
        [r,c,v] = find(tempm);
        if(min(size(r))>0)
            iref(ii)=r(1);jref(ii)=c(1);
        else
            iref(ii)=-1;jref(ii)=-1;
        end
    end
    for ii=1:numref2
        q1=2*pi*norm((refl2(ii,:)./lat2));
        tempm=(abs(qmat(:,:,4)-q1))<0.014;
        tempm1=tempm1+tempm;    
        [r,c,v] = find(tempm);
        if(min(size(r))>0)
            iref2(ii)=r(1);jref2(ii)=c(1);
        else
            iref2(ii)=-1;jref2(ii)=-1;
        end
    end
    tempm1=tempm1>=1;
    impow=imdet.*(1-tempm1)+imax*tempm1;
    if (plotflag==1)
 figure(142);set(142,'Position',[3*imd imb imd imd]);clf;imagesc(impow);axis image;caxis([0 imax/200]);
 title('Image with powder lines');
    for ii=1:numref
        if(iref(ii)>=0)
text(jref(ii),iref(ii),['Si[' num2str(refl(ii,1)) num2str(refl(ii,2)) num2str(refl(ii,3)) ']'],'FontSize',14,'Color',[1,1,1]); 
        end
    end
    for ii=1:numref2
        if(iref2(ii)>=0)
text(jref2(ii),iref2(ii),['Au[' num2str(refl2(ii,1)) num2str(refl2(ii,2)) num2str(refl2(ii,3)) ']'],'FontSize',14,'Color',[1,1,1]); 
        end
    end 
    pause(0.5);
    end
end
%}
%{
imdet2=zeros(size(imdet,1),size(imdet,2));
for ii=1:487
    if((ii>10)&&(ii<30))
        imdet2(:,ii)=imdet(:,ii)>320;
    end
    if((ii>150)&&(ii<210)) 
        imdet2(:,ii)=imdet(:,ii)>545;
    end
     if((ii>250)&&(ii<320)) 
        imdet2(:,ii)=imdet(:,ii)>795;
    end   
end
%}
imdet = double(imdet>40);
for ii=1:size(imdet,1)
    for jj=1:size(imdet,2)
        rrr=sqrt((ii-size(imdet,1))^2+jj^2);
        if(rrr<550)imdet(ii,jj)=0.0;
        end
        %if(rrr>1200)imdet(ii,jj)=0.0;end
    end
end
%imdet=imdet2;        
%imdet=imdet(:,40:100)>1800;
tempm2 = (tempm1-(imdet)).^2;
%dataout=outlist;
%dataout=qlist;
%dataout = sum(sum(tempm2(:,25:450)));
dataout.error=sum(sum(tempm2(:,:)));
%dataout.error=sum(sum(tempm2(:,[5:45,75:116,219:270,328:384])));
dataout.twotheta = acosd(kfmat(:,:,3));
dataout.qmat = qmat;
dataout.gamma = 90-acosd(kfmat(:,:,2));
dataout.twotheta_mot = atand(kfmat(:,:,1)./kfmat(:,:,3));
dataout.gamma_mot = asind(kfmat(:,:,2));
dataout.powder = imdet;
dataout.theory = tempm1;
if(plotflag)
    csvwrite("gamma.csv", dataout.gamma);
    csvwrite("twotheta.csv", dataout.twotheta);
    csvwrite("twotheta_mot.csv", dataout.twotheta_mot);
    csvwrite("gamma_mot.csv", dataout.gamma_mot);
end
end
