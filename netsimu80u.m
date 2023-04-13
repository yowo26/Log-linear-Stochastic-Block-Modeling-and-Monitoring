tic
d1=6;d2=9;d3=12;d4=15;d5=18;d6=20;d=d1+d2+d3+d4+d5+d6;M=6;
%============================community 1================================
s=1;
rng(s);
thta1 = 2*rand(d1,1)+1;thta11 = 2*rand(d1,1)+1;
thta2 = 2*rand(d2,1)+1;thta22 = 2*rand(d2,1)+1;
thta3 = 2*rand(d3,1)+1;thta33 = 2*rand(d3,1)+1;
thta4 = 2*rand(d4,1)+1;thta44 = 2*rand(d4,1)+1;
thta5 = 2*rand(d5,1)+1;thta55 = 2*rand(d5,1)+1;
thta6 = 2*rand(d6,1)+1;thta66 = 2*rand(d6,1)+1;
thta1=thta1/sum(thta1(:,1));thta11=thta11/sum(thta11(:,1));
thta2=thta2/sum(thta2(:,1));thta22=thta22/sum(thta22(:,1));
thta3=thta3/sum(thta3(:,1));thta33=thta33/sum(thta33(:,1));
thta4=thta4/sum(thta4(:,1));thta44=thta44/sum(thta44(:,1));
thta5=thta5/sum(thta5(:,1));thta55=thta55/sum(thta55(:,1));
thta6=thta6/sum(thta6(:,1));thta66=thta66/sum(thta66(:,1));
%============================community parameters========================
gama = zeros(M,M); 
s=502;
rng(s);
nc0=[d1;d2;d3;d4;d5;d6];
for i=1:M
    for j=1:M
        %if i==j
        if (i-4)<1 && (j-4)<1
           gama(i,j)=200+(nc0(i,1)+nc0(j,1))*round(rand(1,1)*10)*10;
        else
           gama(i,j)=200+(nc0(i,1)+nc0(j,1))*round(rand(1,1)*10)*15;
        end
           
         %  gama(i,j)=rand(1,1)+1;
        %else
        %   gama(i,j)=rand(1,1)+1;
        %end
    end
end
%==========================================================================
nn1=zeros(nc0(1,1),nc0(1,1)-1);
nn2=zeros(nc0(2,1),nc0(2,1)-1);
nn3=zeros(nc0(3,1),nc0(3,1)-1);
nn4=zeros(nc0(4,1),nc0(4,1)-1);
nn5=zeros(nc0(5,1),nc0(5,1)-1);
nn6=zeros(nc0(6,1),nc0(6,1)-1);

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
        elseif i==5
           nn5(j,j)=1;
        elseif i==6
           nn6(j,j)=1;
        end
    end
end
nn1(nc0(1,1),:)=-1; nn2(nc0(2,1),:)=-1; nn3(nc0(3,1),:)=-1;nn4(nc0(4,1),:)=-1; nn5(nc0(5,1),:)=-1;nn6(nc0(6,1),:)=-1;

nn7=zeros(M^2,M^2);
for i=1:(M^2)
    nn7(i,i)=1;
end
xx1=blkdiag(nn1,nn1);xx2=blkdiag(nn2,nn2);xx3=blkdiag(nn3,nn3);xx4=blkdiag(nn4,nn4);
xx=blkdiag(nn1,nn2,nn3,nn4,nn5,nn6);
nn=blkdiag(nn1,nn2,nn3,nn4,nn5,nn6,nn1,nn2,nn3,nn4,nn5,nn6,nn7);

%=========================================================================
thta0=[thta1(1:d1,1);thta2(1:d2,1);thta3(1:d3,1);thta4(1:d4,1);thta5(1:d5,1);thta6(1:d6,1)];
thta00=[thta11(1:d1,1);thta22(1:d2,1);thta33(1:d3,1);thta44(1:d4,1);thta55(1:d5,1);thta66(1:d6,1)];
c0=[ones(d1,1);2*ones(d2,1);3*ones(d3,1);4*ones(d4,1);5*ones(d5,1);6*ones(d6,1)];

e2=ones(2,1);

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
        if c0(i,1)==1 && c0(j,1)==1 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<15 && ldaa(j,i)<15
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==2 && c0(j,1)==2 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<21 && ldaa(j,i)<21
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==3 && c0(j,1)==3 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<17 && ldaa(j,i)<17
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==4 && c0(j,1)==4 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<16 && ldaa(j,i)<16
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif c0(i,1)==5 && c0(j,1)==5 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<25 && ldaa(j,i)<25
                k=k+1;
                rho(i,j)=rand(1,1);
                if rho(i,j)>0.3
                   ldaa0(i,j)=ldaa(i,j);
                   ldaa0(j,i)=ldaa(j,i);
                   pi0(i,j)=normcdf(pi(i,j));
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif c0(i,1)==6 && c0(j,1)==6 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<16 && ldaa(j,i)<16
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                pi0(j,i)=normcdf(pi(j,i));
        elseif abs(c0(i,1)-c0(j,1))>0 && ldaa(i,j)>2 && ldaa(j,i)>2 && ldaa(i,j)<15 && ldaa(j,i)<15
                k=k+1;
                rho(i,j)=rand(1,1);
                %if rho(i,j)>0.5
                   ldaa0(i,j)=ldaa(i,j);
                   ldaa0(j,i)=ldaa(j,i);
                   pi0(i,j)=normcdf(pi(i,j));
                   pi0(j,i)=normcdf(pi(j,i));
                %end
        end
    end
end

mup=zeros(d,d);
for i=1:d
    for j=1:d
        if pi0(i,j)>0
            mup(i,j)=(1-pi0(i,j))*(ldaa0(i,j)+ldaa0(j,i));
            %mup(i,j)=ldaa0(i,j)+ldaa0(j,i);
        end
    end
end
%==========================================================================
s=504;
rng(s);
N=500;
NN=100;
mu=[0 0];
pmup=zeros(d,d);
for a=1:NN
          ii=0;
          ob=cell(d,d);
          for i = 1:(d-1)
              for j = (i+1):d 
                  if pi0(i,j)>0
                      ii=ii+1;
                      lda=zeros(2,1);lda(1,1)=ldaa0(i,j);lda(2,1)=ldaa0(j,i);
                      ob{i,j}=zeros(4,N);
                      for k=1:N
                          %sig=[1 rho(i,j);rho(i,j) 1];r=mvnrnd(mu,sig,1);
                          ob{i,j}(1,k)=binornd(1,1-pi0(i,j));
                          %ob{i,j}(2,k)=binornd(1,1-pi0(j,i));
                          if ob{i,j}(1,k)==1 %r(1,1)>pi(i,j)%
                             %ob{i,j}(1,k)=1;
                             ob{i,j}(3,k)=random('Poisson',lda(1,1));
                             ob{i,j}(4,k)=random('Poisson',lda(2,1));
                          else
                             %ob{i,j}(1,k)=0; 
                             ob{i,j}(3,k)=0;
                             ob{i,j}(4,k)=0;
                          end
                          %if ob{i,j}(2,k)==1%r(1,2)>pi(j,i)%ob{i,j}(2,k)==1
                             %ob{i,j}(2,k)=1;
                           %  ob{i,j}(4,k)=random('Poisson',lda(2,1));
                          %else
                             %ob{i,j}(2,k)=0; 
                            % ob{i,j}(4,k)=0;
                          %end
                      end
                  end
              end
          end
      MaxA=zeros(25,1);
      MaxB=zeros(25,1);
      pi1=zeros(d,d);
      Alph=0.02*ones(d-M,1);
      Thta=0.002*ones(2*d+M^2-2*M,1)-0.001;
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
                     bb=zeros(2,1);bb(1,1)=length(find(ob{i,j}(3,1:N)==0 & ob{i,j}(4,1:N)==0));
                     bb(2,1)=N-bb(1,1);
                     obb=zeros(2,1);obb(1,1)=sum(ob{i,j}(3,1:N));obb(2,1)=sum(ob{i,j}(4,1:N));
                     mm=zeros(2,2*d+M^2);
                     mm(1,i)=1;mm(1,d+j)=1;mm(2,d+i)=1;mm(2,j)=1;
                     mm(1,2*d+(c0(j,1)-1)*M+c0(i,1))=1;mm(2,2*d+(c0(i,1)-1)*M+c0(j,1))=1; 
                     mn=mm*nn;Lda=exp(mn*Thta);
                     Ldaa(i,j)=Lda(1,1); Ldaa(j,i)=Lda(2,1);
                     if exp(e2'*Lda)-1>0
                        f0 = f0 + mn'*(obb - bb(2,1)*Lda/(1-exp(-e2'*Lda)));
                        f1 = f1 - bb(2,1)*((mn'*diag(Lda)* mn)/(1-exp(-e2'*Lda)) - (mn'*(Lda * Lda')*mn)*exp(-e2'*Lda)/(1-exp(-e2'*Lda))^2);
                     end
                     ww=zeros(2,d);ww(1,i)=1;ww(2,j)=1;eee=[1 1;0 0];
                     wx=eee*ww*xx;
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
      p1=n01/N;
      pmup0=zeros(d,d);
      for i=1:(d-1)
          for j=(i+1):d
             if pi0(i,j)>0 && Ldaa(i,j)+Ldaa(j,i)>0 
                pmup0(i,j)=(1-pi1(i,j))*(Ldaa(i,j)+Ldaa(j,i))/(1-exp(-Ldaa(i,j)-Ldaa(j,i)));
             end
          end
     end
      pmup=pmup+(pmup0-mup).^2;
end

pMaxA=zeros(d*d,1);
cc=0;
for i=1:(d-1)
    for j=(i+1):d
        if pi0(i,j)>0
            cc=cc+1;
            pMaxA(cc,1)=sqrt(pmup(i,j)/NN);
        end 
    end
end



toc