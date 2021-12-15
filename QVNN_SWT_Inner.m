clear
load Inner
Gh=v;
gh=Gh(end-51940+1:end);
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db1');
%SWT
delay=100;
for ii=1:length(gh)-delay
xs=gh(ii:delay+ii-1);
xl1=conv(xs,Lo_D,'valid')';
xh1=conv(xs,Hi_D,'valid')';
xl2=conv(xl1,Lo_D,'valid')';
xh2=conv(xl1,Hi_D,'valid')';
mg(:,ii)=[xl1(end), xh1(end), xl2(end), xh2(end)];
end
mg=mg';
%Normalization
for i=1:4, xN(:,i)=normalis(mg(:,i),0.1,0.9); end
% To quaternion
yy=quaternion(xN(1:end,1),xN(1:end,2),xN(1:end,3),xN(1:end,4));
yN2=yy;
% Initialization
Nb_Iu = 1;Nb_Hu = 40; Nb_Ou = 1;
mt_1Z=0; mt_1W=0; mt_1z0=0; mt_1w0=0;
vt_1hatZ=0; vt_1hatz0=0; vt_1hatW=0; vt_1hatw0=0;
ep=quaternion(eps, eps, eps ,eps);
alpha=0.001;
beta1=0;
beta2=0.999;

W=1*randq(Nb_Iu,Nb_Hu);w0=1*randq(1,Nb_Hu);Z=1*randq(Nb_Hu,Nb_Ou);z0=1*randq(1,Nb_Ou);
k=0;
vt_1Z=zerosq(size(Z));
vt_1z0=zerosq(size(z0));
vt_1W=zerosq(size(W));
vt_1w0=zerosq(size(w0));
N_iter=100; 
delay1=delay-3;
mu=5;
nn= 38880;
%Training
while k <= N_iter
    k=k+1;
    mtz1=zeros(size(Z));
    mtw1=zeros(size(W));
    mtz0=zeros(size(z0));
    mtw0=zeros(size(w0));
    tic
    for i=2:nn
        XX=yN2(i-1,:);
        U=XX*W+w0;
        hh=[ReLU(realq1(U))' ReLU(x(U))' ReLU(y(U))' ReLU(z(U))'];
        h=quaternion(hh(:,1),hh(:,2),hh(:,3),hh(:,4));
        hc=conj(h, 'hamilton');
        S=hc'*Z+z0;
        St=conj(S, 'hamilton')';
        oo=[logsig(realq1(St)) logsig(x(St)) logsig(y(St)) logsig(z(St))];
        o=quaternion(oo(:,1),oo(:,2),oo(:,3),oo(:,4));
        
        or(i,:,:)=[realq1(o) x(o) y(o) z(o)]';
        ot=conj(o, 'hamilton')';
        e=yN2(i,:)-ot;
         et=conj(e, 'hamilton')';
        dz0=[(realq1(et).*(1-realq1(o)).*realq1(o)), ...
            x(et).*(1-x(o)).*x(o),...
            y(et).*(1-y(o)).*y(o),...
            z(et).*(1-z(o)).*z(o)];
        deltaz0=quaternion(dz0(:,1),dz0(:,2),dz0(:,3),dz0(:,4));
        deltaZ =  hc*conj(deltaz0, 'hamilton')';
        mtZ=beta1*mt_1Z+(1-beta1)*(-conj(deltaZ, 'hamilton')');
        vtZhat=beta2*vt_1hatZ+(1-beta2)*(-deltaZ.*deltaZ);
        vtZr=max(realq1(vt_1Z),realq1(vtZhat));
        vtZx=max(x(vt_1Z),x(vtZhat));
        vtZy=max(y(vt_1Z),y(vtZhat));
        vtZz=max(z(vt_1Z),z(vtZhat));
        vtZ=quaternion(vtZr, vtZx, vtZy, vtZz);
        vhZ=vtZ;
        mhZ=mtZ;
        mt_1Z=mtZ; vt_1Z=vtZ;vt_1hatZ=vtZhat;
        mtz0=beta1*mt_1z0+(1-beta1)*(-deltaz0);
        vtz0hat=beta2*vt_1hatz0+(1-beta2)*(-deltaz0.*deltaz0);
        vtz0r=max(realq1(vt_1z0),realq1(vtz0hat));
        vtz0x=max(x(vt_1z0),x(vtz0hat));
        vtz0y=max(y(vt_1z0),y(vtz0hat));
        vtz0z=max(z(vt_1z0),z(vtz0hat));
        vtz0=quaternion(vtz0r, vtz0x, vtz0y, vtz0z);
        vhz0=vtz0;
        mhz0=mtz0;
        mt_1z0=mtz0; vt_1z0=vtz0;vt_1hatz0=vtz0hat;
        QZ=[realq1(vhZ), x(vhZ) , y(vhZ) , z(vhZ)];
        QZa=sqrt(sum(QZ'.*QZ'));
        [r1Z, r2Z, r3Z] = quat2angle(QZ);
        RZ=angle2quat(r1Z/2, r2Z/2, r3Z/2);
        RZq=(sqrt(QZa')).*quaternion(RZ(:,1), RZ(:,2), RZ(:,3), RZ(:,4));
        QZrootq=mu*RZq;
        QPZ=(1/mu)*qqlog(1+exp(QZrootq));
        Z=Z- alpha*conj(mhZ,'hamilton')'./QPZ;
        Qz0=[realq1(vhz0), x(vhz0), y(vhz0) , z(vhz0) ];
        Qz0a=sqrt(sum(Qz0'.*Qz0'));
        [r1z0, r2z0, r3z0] = quat2angle(Qz0);
        Rz0=angle2quat(r1z0/2, r2z0/2, r3z0/2);
        Rz0q=(sqrt(Qz0a')).*quaternion(Rz0(:,1), Rz0(:,2), Rz0(:,3), Rz0(:,4));
        Qz0rootq=mu*Rz0q;
        QPz0=(1/mu)*qqlog(1+exp(Qz0rootq));
        z0=z0- alpha*mhz0./QPz0;
        Ze=conj(conj(deltaz0, 'hamilton')'*conj(conj(Z, 'hamilton'), 'hamilton')', 'hamilton')';
        dw0= [(realq1(h)>0).*realq1(Ze),...
            +(x(h)>0).*x(Ze),...
            +(y(h)>0).*y(Ze),...
            +(z(h)>0).*z(Ze)];
        deltaw0=quaternion(dw0(:,1),dw0(:,2),dw0(:,3),dw0(:,4));
        deltaw0=conj(deltaw0, 'hamilton')';
        deltaW=conj(conj(XX, 'hamilton'),'hamilton')'*deltaw0;
        mtW=beta1*mt_1W+(1-beta1)*(-deltaW);
        vtWhat=beta2*vt_1hatW+(1-beta2)*(-deltaW.*deltaW);
        vtWr=max(realq1(vt_1W),realq1(vtWhat));
        vtWx=max(x(vt_1W),x(vtWhat));
        vtWy=max(y(vt_1W),y(vtWhat));
        vtWz=max(z(vt_1W),z(vtWhat));
        vtW=quaternion(vtWr, vtWx, vtWy, vtWz);
        vhW=vtW;
        mhW=mtW;
        mt_1W=mtW; vt_1W=vtW;vt_1hatW=vtWhat;
        mtw0=beta1*mt_1w0+(1-beta1)*(-deltaw0);
        vtw0hat=beta2*vt_1hatw0+(1-beta2)*(-deltaw0.*deltaw0);
        
        vtw0r=max(realq1(vt_1w0),realq1(vtw0hat));
        vtw0x=max(x(vt_1w0),x(vtw0hat));
        vtw0y=max(y(vt_1w0),y(vtw0hat));
        vtw0z=max(z(vt_1w0),z(vtw0hat));
        vtw0=quaternion(vtw0r, vtw0x, vtw0y, vtw0z);
        
        vhw0=vtw0;
        mhw0=mtw0;
        mt_1w0=mtw0; vt_1w0=vtw0;vt_1hatw0=vtw0hat;
        QW=[realq1(vhW)', x(vhW)' , y(vhW)' , z(vhW)'];
        QWa=sqrt(sum(QW'.*QW'));
        [r1W, r2W, r3W] = quat2angle(QW);
        RW=angle2quat(r1W/2, r2W/2, r3W/2);
        
        RWq=(sqrt(QWa')).*quaternion(RW(:,1), RW(:,2), RW(:,3), RW(:,4));
        QWrootq=mu*RWq;
        QPW=(1/mu)*qqlog(1+ exp(QWrootq));
        W=W- alpha*mhW./conj(QPW, 'hamilton')';
        Qw0=[realq1(vhw0)', x(vhw0)' , y(vhw0)' , z(vhw0)' ];
        Qw0a=sqrt(sum(Qw0'.*Qw0'));
        [r1w0, r2w0, r3w0] = quat2angle(Qw0);
        Rw0=angle2quat(r1w0/2, r2w0/2, r3w0/2);
        Rw0q=(sqrt(Qw0a')).*quaternion(Rw0(:,1), Rw0(:,2), Rw0(:,3), Rw0(:,4));
        Qw0rootq=mu*Rw0q;
        QPw0=(1/mu)*qqlog(1+ exp(Qw0rootq));
        w0=w0- alpha*mhw0./conj(QPw0, 'hamilton')';
        V1=((realq1(o) - 0.1)*(max(mg(:,1)) - min(mg(:,1))))/0.8+min(mg(:,1)); Q1(i)=mg(i,1)-V1; RE1(i)=V1;
        V2=((x(o) - 0.1)*(max(mg(:,2)) - min(mg(:,2))))/0.8+min(mg(:,2));Q2(i)=mg(i,2)-V2; RE2(i)=V2;
        V3=((y(o) - 0.1)*(max(mg(:,3)) - min(mg(:,3))))/0.8+min(mg(:,3));Q3(i)=mg(i,3)-V3; RE3(i)=V3;
        V4=((z(o) - 0.1)*(max(mg(:,4)) - min(mg(:,4))))/0.8+min(mg(:,4));Q4(i)=mg(i,4)-V4; RE4(i)=V4;
    end
    time(k)=toc;
    for ii=1:length(RE1)
        oh2=conv(RE4(ii),Hi_D)';
        ol2=conv(RE3(ii),Lo_D)';
        olt1=(oh2+ol2);
        ol1=conv(olt1,Lo_D)';
        oh1=conv(RE2(ii),Hi_D)';
        Ye=(oh1+ol1');
        recon1(ii)=Ye(end);
    end
    E=gh(1+delay1:nn+delay1)-recon1(1:nn)';
    moy(k)=sqrt(sum(E.^2)/nn);
    fprintf('%d          %f\n',k,moy(k));
    E=0;
end
% SWT coefficient prediction 
for i=2:51840-1
    XX=yN2(i-1,:);
    U=XX*W+w0;
    hh=[ReLU(realq1(U))' ReLU(x(U))' ReLU(y(U))' ReLU(z(U))'];
    h=quaternion(hh(:,1),hh(:,2),hh(:,3),hh(:,4));
    hc=conj(h, 'hamilton');
    S=hc'*Z+z0;
    St=conj(S, 'hamilton')';
    oo=[logsig(realq1(St)) logsig(x(St)) logsig(y(St)) logsig(z(St))];
    o=quaternion(oo(:,1),oo(:,2),oo(:,3),oo(:,4));
    or(i,:,:)=[realq1(o) x(o) y(o) z(o)]';
    ot=conj(o, 'hamilton')';
    e=yN2(i,:)-ot;
    et=conj(e, 'hamilton')';
    V1=((realq1(o) - 0.1)*(max(mg(:,1)) - min(mg(:,1))))/0.8+min(mg(:,1)); Q1(i)=mg(i,1)-V1; RE1(i)=V1;
    V2=((x(o) - 0.1)*(max(mg(:,2)) - min(mg(:,2))))/0.8+min(mg(:,2));Q2(i)=mg(i,2)-V2; RE2(i)=V2;
    V3=((y(o) - 0.1)*(max(mg(:,3)) - min(mg(:,3))))/0.8+min(mg(:,3));Q3(i)=mg(i,3)-V3; RE3(i)=V3;
    V4=((z(o) - 0.1)*(max(mg(:,4)) - min(mg(:,4))))/0.8+min(mg(:,4));Q4(i)=mg(i,4)-V4; RE4(i)=V4;
    
end
NE= 38880;
NE1=51840-1;
% Using the predicted SWT to have the predited wind speed
for ii=1:NE1
    xs=gh(ii:delay+ii-2);
    xl1=conv(xs,Lo_D,'valid')';
    xh1=conv(xs,Hi_D,'valid')';
    xl2=conv(xl1,Lo_D,'valid')';
    xh2=conv(xl1,Hi_D,'valid')';
    xh2(end)=RE4(ii);
    xl2(end)=RE3(ii);
    xh1(end)=RE2(ii);
    oh2=conv(xh2,Hi_D)';
    ol2=conv(xl2,Lo_D)';
    olt1=(oh2+ol2);
    ol1=conv(olt1,Lo_D)';
    oh1=conv(xh1,Hi_D)';
    Ye=(oh1+ol1');
    r(ii)=Ye(end);
    t(ii)=gh(ii+delay1);
   
end
% Compute the mertics
NE= 38880;
NE1=51840-1;
Re=r;
n=NE1-NE; % Three months = 24*30*3*6 = 12959
% Compute the metrics
y1=t(NE+1:NE1)';
yhat=Re(NE+1:NE1)';
E=y1-yhat;
MSEv=sum((abs(E)).^2/n)
MAE=sum((abs(E))/n)
MAPE=100*sum(abs(E)./abs(y1))/n
RMSE=sqrt(MSEv)