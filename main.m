clear;clc;close all;

pic1 = imread('1st.jpg');
pic2 = imread('2nd.jpg');

point1  = [734,1653;2416,43;3620,1290;1907,2889];
point2  = [1112,455;3647,476;3560,2359;1104,2241];
hLen    = 20;
tLen    = 2*hLen+1;
dAngle  = 45;
nCorner = 4;
gaussianSigma = 4;
gaussianFilter = fspecial('gaussian', [tLen tLen], 1.5*gaussianSigma);

for cnt = 1:nCorner
    pIdx = ['p',num2str(cnt)];
    
    % Make block and convert to grayscale
    % Do the gaussian blurring
    sub1.(pIdx) = rgb2gray(pic1(point1(cnt,2)+(-hLen:hLen), point1(cnt,1)+(-hLen:hLen), :));
    sub1.(pIdx) = imgaussfilt(sub1.(pIdx),gaussianSigma);
    
    sub2.(pIdx) = rgb2gray(pic2(point2(cnt,2)+(-hLen:hLen), point2(cnt,1)+(-hLen:hLen), :));
    sub2.(pIdx) = imgaussfilt(sub2.(pIdx),gaussianSigma);
    
    hist1.(pIdx) = zeros(1,360/dAngle);
    hist2.(pIdx) = zeros(1,360/dAngle);
end


for cnt = 1:nCorner
    pIdx = ['p',num2str(cnt)];
    % 1.5 times of scale sigma
    
    [Gmag1,Gdir1] = imgradient(sub1.(pIdx),'sobel');
    hValue1 = int16(ceil((Gdir1+180.0001)/dAngle));
    Gmag1 = Gmag1.*gaussianFilter;
    [Gmag2,Gdir2] = imgradient(sub2.(pIdx),'sobel');
    hValue2 = int16(ceil((Gdir2+180.0001)/dAngle));
    Gmag2 = Gmag2.*gaussianFilter;
    
    for cnt2 = 1:tLen
        for cnt3 = 1:tLen
            hist1.(pIdx)(hValue1(cnt2,cnt3)) = hist1.(pIdx)(hValue1(cnt2,cnt3)) + Gmag1(cnt2,cnt3);
            hist2.(pIdx)(hValue2(cnt2,cnt3)) = hist2.(pIdx)(hValue2(cnt2,cnt3)) + Gmag2(cnt2,cnt3);
        end
    end

end
2;
for cnt = 1:nCorner
    pIdx = ['p',num2str(cnt)];
    [~,temp]=sort(hist1.(pIdx));
    hIdx1(cnt) = temp(end)*dAngle-dAngle/2;
    [~,temp]=sort(hist2.(pIdx));
    hIdx2(cnt) = temp(end)*dAngle-dAngle/2;
end
2;
startPoint = hLen-8+1;
spSize = 4;
spNum = 16;
spAngle = 45;
nAngle = 360/spAngle;
gaussianFilter = fspecial('gaussian', [spSize spSize], 1.5*gaussianSigma);
for cnt = 1:nCorner
    pIdx = ['p',num2str(cnt)];
    fHist1.(pIdx) = zeros(1,spNum*nAngle);
    fHist2.(pIdx) = zeros(1,spNum*nAngle);
    
    % do every sub patch
    counter = 0;
    for cnt2 = startPoint : spSize : startPoint + spSize*3
        for cnt3 = startPoint : spSize : startPoint + spSize*3
            counter = counter+1;
            xRange = cnt2 - 1*spSize + (1:spSize);
            yRange = cnt3 - 1 + (1:spSize);
            temp1 = sub1.(pIdx)(xRange,yRange);
            temp2 = sub2.(pIdx)(xRange,yRange);
            
            % get each sub patch gradient
            [Gmag1,Gdir1] = imgradient(temp1,'sobel');
            [Gmag2,Gdir2] = imgradient(temp2,'sobel');
            Gdir1 = Gdir1 + 180.0001 - hIdx1(cnt);
            Gdir2 = Gdir2 + 180.0001 - hIdx2(cnt);
            Gmag1 = Gmag1.*gaussianFilter;
            Gmag2 = Gmag2.*gaussianFilter;
            % to make every angle positive
            for cnt4 = 1:spSize
                for cnt5 = 1:spSize
                    while(Gdir1(cnt4,cnt5)<=0)
                        Gdir1(cnt4,cnt5) = Gdir1(cnt4,cnt5) + 360.0;
                    end
                    while(Gdir2(cnt4,cnt5)<=0)
                        Gdir2(cnt4,cnt5) = Gdir2(cnt4,cnt5) + 360.0;
                    end
                end
            end
            
            % convert to angle index
            hValue1 = int16(ceil(Gdir1/spAngle));
            hValue2 = int16(ceil(Gdir2/spAngle));
            tHist1 = zeros(1,nAngle);
            tHist2 = zeros(1,nAngle);
            for cnt4 = 1:spSize
                for cnt5 = 1:spSize
                    tHist1(hValue1(cnt4,cnt5)) = tHist1(hValue1(cnt4,cnt5)) + Gmag1(cnt4,cnt5);
                    tHist2(hValue2(cnt4,cnt5)) = tHist2(hValue2(cnt4,cnt5)) + Gmag2(cnt4,cnt5);
                end
            end
            3;
            fHist1.(pIdx)((counter-1)*nAngle + (1:nAngle)) = tHist1;
            fHist2.(pIdx)((counter-1)*nAngle + (1:nAngle)) = tHist2;
        end
    end
    figure();
    plot(fHist1.(pIdx));
    figure();
    plot(fHist2.(pIdx));
end

for cnt = 1:nCorner
    pIdx1 = ['p',num2str(cnt)];
    tempo1 = fHist1.(pIdx1);
    for cnt2 = 1:nCorner
        pIdx2 = ['p',num2str(cnt2)];
        tempo2 = fHist2.(pIdx2);
        result(cnt,cnt2) = sqrt(sum((tempo1-tempo2).^2));
    end
end
disp(result);

