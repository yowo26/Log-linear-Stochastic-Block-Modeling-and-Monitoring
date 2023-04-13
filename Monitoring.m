tic
d1=6;d2=8;d3=11;d4=15;d=d1+d2+d3+d4;M=4;
%============================community 1================================
s=500;
rng(s);
thta1 = rand(d1,1);thta11 = rand(d1,1);
thta2 = rand(d2,1);thta22 = rand(d2,1);
thta3 = rand(d3,1);thta33 = rand(d3,1);alph0=thta3;
thta4 = rand(d4,1);thta44 = rand(d4,1);
thta1=thta1/sum(thta1(:,1));thta11=thta11/sum(thta11(:,1));
thta2=thta2/sum(thta2(:,1));thta22=thta22/sum(thta22(:,1));
thta3=thta3/sum(thta3(:,1));thta33=thta33/sum(thta33(:,1));
thta4=thta4/sum(thta4(:,1));thta44=thta44/sum(thta44(:,1));
%============================community parameters========================
gama = zeros(M,M); 
s=502;
rng(s);
nc0=[d1;d2;d3;d4];
for i=1:M
    for j=1:M
        gama(i,j)=200+(nc0(i,1)+nc0(j,1))*round(rand(1,1)*10)*10; 
    end
end

%=========================================================================
thta0=[thta1(1:d1,1);thta2(1:d2,1);thta3(1:d3,1);thta4(1:d4,1)];
thta00=[thta11(1:d1,1);thta22(1:d2,1);thta33(1:d3,1);thta44(1:d4,1)];
c0=[ones(d1,1);2*ones(d2,1);3*ones(d3,1);4*ones(d4,1)];

s=504;
rng(s);
cc=0;
rho=zeros(d,d);
pi=zeros(d,d);
ldaa=zeros(d,d);
pi0=zeros(d,d);
ldaa0=zeros(d,d);
k=0;
for i = 1:(d-1)
    for j = (i+1):d
        cc=cc +1;
        pi(i,j)=2*rand(1,1)-1;
        pi(j,i)=2*rand(1,1)-1;
        ldaa(i,j)=thta0(i,1)*thta00(j,1)*gama(c0(i,1),c0(j,1));
        ldaa(j,i)=thta0(j,1)*thta00(i,1)*gama(c0(j,1),c0(i,1));
        if c0(i,1)==1 && c0(j,1)==1 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<20 && ldaa(j,i)<20
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==2 && c0(j,1)==2 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<10 && ldaa(j,i)<10
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==3 && c0(j,1)==3 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<20 && ldaa(j,i)<20
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==4 && c0(j,1)==4 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<18 && ldaa(j,i)<18
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif abs(c0(i,1)-c0(j,1))>0 && ldaa(i,j)>2 && ldaa(j,i)>2 && ldaa(i,j)<15 && ldaa(j,i)<15
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        end
    end
end

mup=zeros(d,d);
for i=1:d
    for j=1:d
        if pi0(i,j)>0
            mup(i,j)=(1-pi0(i,j))*ldaa0(i,j);
            mup(j,i)=(1-pi0(j,i))*ldaa0(j,i);
        end
    end
end

%======================================================================

ee=[1 1;0 0];
e2=ones(2,1);
nn1=zeros(nc0(1,1),nc0(1,1)-1);
nn2=zeros(nc0(2,1),nc0(2,1)-1);
nn3=zeros(nc0(3,1),nc0(3,1)-1);
nn4=zeros(nc0(4,1),nc0(4,1)-1);

for i=1:M
    for j=1:(nc0(i,1)-1)
        if i==1
           nn1(j,j)=1;
        elseif i==2
           nn2(j,j)=1;
        elseif i==3
           nn3(j,j)=1;
        elseif i==4
           nn4(j,j)=1;
        end
    end
end
nn1(nc0(1,1),:)=-1; nn2(nc0(2,1),:)=-1; nn3(nc0(3,1),:)=-1;nn4(nc0(4,1),:)=-1; 

nn5=zeros(M^2,M^2);
for i=1:(M^2)
    nn5(i,i)=1;
end
xx1=blkdiag(nn1,nn1);xx2=blkdiag(nn2,nn2);xx3=blkdiag(nn3,nn3);xx4=blkdiag(nn4,nn4);
xx=blkdiag(nn1,nn2,nn3,nn4);
nn=blkdiag(nn1,nn2,nn3,nn4,nn1,nn2,nn3,nn4,nn5);

%==========================================================================
s=504;
rng(s);
N=2000;
mu=[0 0];
ii=0;

          ob=cell(d,d);
          for i = 1:(d-1)
              for j = (i+1):d 
                  if pi0(i,j)>0
                      lda=zeros(2,1);lda(1,1)=ldaa0(i,j);lda(2,1)=ldaa0(j,i);
                      ob{i,j}=zeros(3,N);
                      for k=1:N
                          %sig=[1 rho(i,j);rho(i,j) 1];r=mvnrnd(mu,sig,1);
                          ob{i,j}(1,k)=binornd(1,1-pi0(i,j));
                          if ob{i,j}(1,k)==1
                             ob{i,j}(2,k)=random('Poisson',lda(1,1));
                             ob{i,j}(3,k)=random('Poisson',lda(2,1));
                          else
                             ob{i,j}(2,k)=0; 
                             ob{i,j}(3,k)=0;
                          end
                      end
                  end
              end
          end
      MaxA=zeros(25,1);
      MaxB=zeros(25,1);
      pi1=zeros(d,d);
      Alph=0.02*ones(d-M,1);
      Thta=0.001*ones(2*d+M^2-2*M,1);
      b=0;Maxm=1;
      while b<25 && Maxm>0.02
          b=b+1;
          n01=zeros(d,d);
          af0=zeros(d-M,1);
          af1=zeros(d-M,d-M);
          f0=zeros(2*d+M^2-2*M,1);
          f1=zeros(2*d+M^2-2*M,2*d+M^2-2*M);
          Ldaa=zeros(d,d);
          for i = 1:(d-1)
              for j = (i+1):d 
                  if pi0(i,j)>0
                     bb=zeros(2,1);
                     for k=1:N
                         obb=zeros(2,1);obb(1,1)=ob{i,j}(2,k);obb(2,1)=ob{i,j}(3,k);
                         if sum(obb(:,1))==0
                            bb(1,1)=bb(1,1)+1;
                         else
                            bb(2,1)=bb(2,1)+1;
                            n01(i,j)=n01(i,j)+1;
                            n01(j,i)=n01(j,i)+1;
                            mm=zeros(2,2*d+M^2);
                            mm(1,i)=1;mm(1,d+j)=1;mm(2,d+i)=1;mm(2,j)=1;
                            mm(1,2*d+(c0(j,1)-1)*M+c0(i,1))=1;mm(2,2*d+(c0(i,1)-1)*M+c0(j,1))=1; 
                            mn=mm*nn;Lda=exp(mn*Thta);
                            Ldaa(i,j)=Lda(1,1); Ldaa(j,i)=Lda(2,1);
                            if exp(e2'*Lda)-1>0
                               f0 = f0 + mn'*obb - mn'*Lda/(1-exp(-e2'*Lda));
                               f1 = f1 - mn'*diag(Lda)* mn/(1-exp(-e2'*Lda)) + mn'*(Lda * Lda')*mn*exp(-e2'*Lda)/(1-exp(-e2'*Lda))^2;
                            end
                         end
                     end
                     ww=zeros(2,d);ww(1,i)=1;ww(2,j)=1;
                     wx=ee*ww*xx;
                     pp=exp(wx*Alph-log(e2'*exp(wx*Alph)));
                     diap=diag(pp)-pp*pp';
                     af0 = af0 + wx'*(bb(:,1)-N*pp(:,1));
                     af1 = af1 - N*wx'* diap * wx;
                     pi1(i,j)=pp(1,1);
                     pi1(j,i)=pi1(i,j);
                  end
              end
          end
          Thta1=Thta;
          Thta=Thta-f1\f0;
          MaxA(b,1)=sum((Thta-Thta1).^2);
          Maxm=MaxA(b,1);
          if b<8
             Alph1=Alph;
             Alph=Alph-af1\af0;
             MaxB(b,1)=sum((Alph-Alph1).^2);
          end
      end

%=========================================================================

shft=[0; -0.05; -0.10; -0.15; -0.20; -0.25; 0.05; 0.10; 0.15; 0.20; 0.25];

s=510;
rng(s);
nSimu=1000;

lmd=0.1;

for aa=1:11
    dRL=zeros(nSimu,1);
    %pi11=pi0;
    %for i=1:d
    %    pi11(25,i)=pi11(25,i)*(1+shft(aa,1));
    %    pi11(i,25)=pi11(i,25)*(1+shft(aa,1));
    %end
    alph1=alph0;
    alph1(11,1)=alph1(11,1)+shft(aa,1);
    alph1=alph1/sum(alph1(:,1));
    alph11=[thta1(1:d1,1);thta2(1:d2,1);alph1(1:d3,1);thta4(1:d4,1)];
    %gama1=gama;
    %gama1(3,2)=gama1(3,2)*(1+shft(aa,1));
    ldaa1=zeros(d,d);
    for i = 1:(d-1)
        for j = (i+1):d
            if pi0(i,j)>0
               ldaa1(i,j)=thta00(j,1)*alph11(i,1)*gama(c0(i,1),c0(j,1));
               ldaa1(j,i)=thta00(i,1)*alph11(j,1)*gama(c0(j,1),c0(i,1));
            end
        end
    end
    nIndx=0;
for a=1:nSimu
    b=0;
    MaxA1=0;
    zp=zeros(d-M,1);
    zta=zeros(2*d+M^2-2*M,1);
    zff=ff;
    while MaxA1<lmt && b<2000
          b=b+1;
          ob=zeros(d*(d-1)/2,3);
          cc=0;
          f0=zeros(d-M,1);
          f1=zeros(2*d+M^2-2*M,1);
          f2=zeros(2*d+M^2-2*M,2*d+M^2-2*M);
          for i = 1:(d-1)
              for j = (i+1):d 
                  cc=cc+1;
                  if pi0(i,j)>0
                     mm=zeros(2,2*d+M^2);
                     mm(1,i)=1;mm(1,d+j)=1;mm(2,d+i)=1;mm(2,j)=1;
                     mm(1,2*d+(c0(j,1)-1)*M+c0(i,1))=1;mm(2,2*d+(c0(i,1)-1)*M+c0(j,1))=1; 
                     mn=mm*nn;
                     lda=exp(mn*Thta);
                     ww=zeros(2,d);ww(1,i)=1;ww(2,j)=1;
                     wx=ee*ww*xx;
                     
                     oo=zeros(2,1);obb=zeros(2,1);
                     for k=1:1
                         %if b>10.5
                         %   ob(cc,1)=binornd(1,1-pi11(i,j));
                         %else
                            ob(cc,1)=binornd(1,1-pi0(i,j));
                         %end
                          %ob(cc,2)=binornd(1,1-pi0(j,i));
                         if ob(cc,1)==1 && b>10.5
                            ob(cc,2)=random('Poisson',ldaa1(i,j));
                            ob(cc,3)=random('Poisson',ldaa1(j,i));
                         elseif ob(cc,1)==1 && b<10.5
                            ob(cc,2)=random('Poisson',ldaa0(i,j));
                            ob(cc,3)=random('Poisson',ldaa0(j,i));
                         else
                            ob(cc,2)=0;
                            ob(cc,3)=0;
                         end
                         if sum(ob(cc,2:3))>0
                             oo(2,1)=oo(2,1)+1;
                             obb(:,1)=obb(:,1)+(ob(cc,2:3))';
                         else
                             oo(1,1)=oo(1,1)+1;
                         end
                     end
                     if sum(obb(:,1))>0
                        lda0=zeros(2,1);lda0(1,1)=Ldaa(i,j);lda0(2,1)=Ldaa(j,i);
                        f1 = f1 + mn'*obb - oo(2,1)*mn'*lda0/(1-exp(-e2'*lda0));
                        f2 = f2 + oo(2,1)*(mn'*diag(lda0)* mn/(1-exp(-e2'*lda0)) - mn'*(lda0 * lda0')*mn*exp(-e2'*lda0)/(1-exp(-e2'*lda0))^2); 
                     end
                     pp=zeros(2,1);pp(1,1)=pi1(i,j);pp(2,1)=1-pi1(i,j); 
                     f0=f0 + wx'*(oo(:,1)-pp(:,1));
                  end
              end
          end
          zp=(1-lmd)*zp + lmd*f0;
          zta=(1-lmd)*zta + lmd*f1;
          zff=(1-lmd)*zff + lmd*f2;
          MaxA1=(zp')*(gg\zp)+(zta')*(zff\zta);
          MaxA1=MaxA1/10;
    end
    if b>10
        nIndx=nIndx+1;
        dRL(nIndx,1)=b-10;
    end
end
    ARL1=sum(dRL(1:nIndx,1))/nIndx;
    rslt=[aa ARL1 std(dRL(1:nIndx,1))];
    fprintf(fd,'%12.4f\r',rslt);
end
      
toc