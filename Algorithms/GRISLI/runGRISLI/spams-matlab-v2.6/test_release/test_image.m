

Io=double(imread('data/lena.png'))/255;
sigma=20/255;
n=10;
I=Io+(sigma)*randn(size(Io));

X=mexExtractPatches(I,n);
me=mean(X);
X=bsxfun(@minus,X,me);

paramdict.iter=1000;
paramdict.batchsize=1000;
paramdict.mode=1;
paramdict.K=256;
paramdict.verbose=false;
paramdict.lambda=n*n*(1.1)*sigma*sigma;
D=mexTrainDL(X,paramdict);

paramomp.eps=paramdict.lambda;
paramomp.L=paramdict.K;
alpha=mexOMP(X,D,paramomp);

Xhat=D*alpha+ones(size(D,1),1)*me;
Ihat=mexCombinePatches(Xhat,I,n);

subplot(1,3,1);
imagesc(Io); colormap(gray);
subplot(1,3,2);
imagesc(I); colormap(gray);
subplot(1,3,3);
imagesc(Ihat); colormap(gray);

paramfista.lambda=0.1;
paramfista.max_it=200;
paramfista.it0=1;
paramfista.verbose=true;
paramfista.regul='l1';
paramfista.loss='square';
paramfista.L0=0.0001;
sizeMapX=size(I,2)-n+1;
sizeMapY=size(I,1)-n+1;
p=size(D,2);
sizeMap=sizeMapX*sizeMapY;
alpha0=full(alpha);
alpha0=reshape(alpha0,[p sizeMapX sizeMapY]);
alpha=mexConvFista(I,D,alpha0,paramfista);
alpha=reshape(alpha,[p sizeMap]);
