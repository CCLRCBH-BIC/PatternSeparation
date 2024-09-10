function [Cor,Beta,P_Fstat,Const,Res,ResT]=Univaranalysisv2(X,Y,mask,C,tstat_flag,fstat_flag)
%% X and Y should have been standarized already, each column in C is a contrast
% Y is the dataset after mask applied
% compare with v1: F statistics has replaced t statistics
% in v2: the constant regressor is added.
if exist('tstat_flag','var')== 0 | isempty(tstat_flag)
    tstat_flag = 0;
end
if exist('fstat_flag','var')== 0 | isempty(fstat_flag)
    fstat_flag = 0;
end
X = [ones(size(X,1),1),X];
if exist('C','var')==0 || isempty(C)
    C = [1;zeros(size(X,2)-1,1)];
else
    C = [zeros(1,size(C,2));C];
end
[Nreg,Ncon]=size(C);
[tdim,N]=size(Y);

beta=pinv(X)*Y;
Yest=X*beta;
ResT = zeros(tdim,N);

cor=zeros(N,1);
res = zeros(N,1);
const=zeros(N,Ncon);
p_Fstat = zeros(N,1);
for i=1:N
    cor(i)=corr(Y(:,i),Yest(:,i));
    res(i) = mean((Y(:,i)-Yest(:,i)).^2);
    ResT(:,i) = Y(:,i)-Yest(:,i);
    if tstat_flag==0
        const(i,:)=lambdastatistic(X,Y(:,i),1,C,1);
    else
        const(i,:)=Tstatistic(X,Y(:,i),C,beta(:,i));
    end
    if fstat_flag>0
        try
            mdl = fitlm(X(:,2:end),Y(:,i));
            p_Fstat(i) = coefTest(mdl);
        catch
            disp('sth is wrong.')
        end
    end
end
res = sqrt(res);
if exist('mask','var')&& isempty(mask)==0
    sz=size(mask);
    Cor=zeros(prod(sz),1);
    P_Fstat = zeros(size(mask));
    Res = zeros(prod(sz),1);
    Const=zeros([prod(sz),Ncon]);
    Beta=zeros([prod(sz),Nreg-1]);
    Cor(mask(:)==1)=cor;
    Res(mask(:)==1) = res;
    Const(mask(:)==1,:)=const;
    Beta(mask(:)==1,:)=beta(2:end,:)';
    Cor=reshape(Cor,sz);
    Res = reshape(Res,sz);
    Const=reshape(Const,[sz,Ncon]);
    Beta=reshape(Beta,[sz,Nreg-1]);
    P_Fstat(mask>0) = p_Fstat;
else
    Cor=cor;Const=const;Beta=beta(2:end,:)';Res = res;P_Fstat = p_Fstat;
end
end
    
    
