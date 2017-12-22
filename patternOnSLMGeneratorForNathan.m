function patternOnSLMGeneratorForNathan
    dr=13.58;
    filePath='C:\medowlark\Nathan';
    if ~exist(filePath,'dir')
        mkdir(filePath)
    end
    piValueAll=100:350;
    for ii=1:length(piValueAll)
        phasePattern=patternOnSLMGenerator(2*dr,512,512,piValueAll(ii));
        file=['dr',num2str(dr,'%.2f'),'_pi_',num2str(piValueAll(ii),'%03d'),'.bmp'];
        imwrite(phasePattern,fullfile(filePath,file));
    end
    
    
function phasePattern=patternOnSLMGenerator(S,pixelNumberX,pixelNumberY,piValue)
% patternOnSLMGenerator is to generate binary pattern used for SLM;
% S is the period of binary pattern on SLM;
% pixelNumberX is the dimension of SLM on x axis;
% pixelNumberY is the dimension of SLM on y axis;
%---------example------------------
% S=24.2;
% pixelNumberX=512;
% pixelNumberY=512;
% phasePattern=patternOnSLMGenerator(S,pixelNumberX,pixelNumberY);   

if nargin==0
    S=24.2;
    pixelNumberX=512;
    pixelNumberY=512;
end
m=1:pixelNumberY;
n=1:pixelNumberX;
[nn,mm]=meshgrid(n,m);
yy=1*(mm-round(pixelNumberY/2));
xx=1*(nn-round(pixelNumberX/2));
rr=yy.^2+xx.^2;rr=sqrt(rr);
beta=1*pi/4;
phasePattern_tmp=cos(2*pi*rr/S-beta)+eps;
phasePattern=zeros(pixelNumberY,pixelNumberX);
phasePattern(phasePattern_tmp<0)=0;
phasePattern(phasePattern_tmp>0)=1;
phasePattern=uint8(phasePattern*piValue);
% figure(4);clf;imshow(phasePattern);title('binary phase pattern on SLM')
