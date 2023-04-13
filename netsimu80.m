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
        if (i-4)<1 && (j-4)<1
           gama(i,j)=200+(nc0(i,1)+nc0(j,1))*round(rand(1,1)*10)*10;
        else
           gama(i,j)=200+(nc0(i,1)+nc0(j,1))*round(rand(1,1)*10)*15;
        end
    end
end
%==========================================================================
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
kc=0;
for i = 1:(d-1)
    for j = (i+1):d
        cc=cc +1;
        kk=rand(1,1);
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
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif c0(i,1)==2 && c0(j,1)==2 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<21 && ldaa(j,i)<21
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif c0(i,1)==3 && c0(j,1)==3 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<17 && ldaa(j,i)<17
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif c0(i,1)==4 && c0(j,1)==4 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<16 && ldaa(j,i)<16
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif c0(i,1)==5 && c0(j,1)==5 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<25 && ldaa(j,i)<25
                k=k+1;
                rho(i,j)=rand(1,1);
                if rho(i,j)>0.3
                   ldaa0(i,j)=ldaa(i,j);
                   ldaa0(j,i)=ldaa(j,i);
                   pi0(i,j)=normcdf(pi(i,j));
                   %kk=rand(1,1);
                   if kk>0.4
                      kc=kc+1;
                      pi0(j,i)=normcdf(pi(j,i));
                   end
                end
        elseif c0(i,1)==6 && c0(j,1)==6 && ldaa(i,j)>1 && ldaa(j,i)>1 && ldaa(i,j)<16 && ldaa(j,i)<16
                k=k+1;
                rho(i,j)=rand(1,1);
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
        elseif abs(c0(i,1)-c0(j,1))>0 && ldaa(i,j)>2 && ldaa(j,i)>2 && ldaa(i,j)<15 && ldaa(j,i)<15
                k=k+1;
                rho(i,j)=rand(1,1);
                %if rho(i,j)>0.5
                ldaa0(i,j)=ldaa(i,j);
                ldaa0(j,i)=ldaa(j,i);
                pi0(i,j)=normcdf(pi(i,j));
                %kk=rand(1,1);
                if kk>0.4
                   kc=kc+1;
                   pi0(j,i)=normcdf(pi(j,i));
                end
                %end
        end
    end
end

mup=zeros(d,d);
for i=1:d
    for j=1:d
        if pi0(i,j)>0
            if pi0(j,i)==0
               mup(i,j)=(1-pi0(i,j))*(ldaa0(i,j)+ldaa0(j,i)); 
            else
               mup(i,j)=(1-pi0(i,j))*ldaa0(i,j)+(1-pi0(j,i))*ldaa0(j,i);
            end
        end
    end
end