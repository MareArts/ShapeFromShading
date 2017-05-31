clc;
clear all;

%%광원 설정
light1=[0 0 40];
light2=[20 5 40];
light3=[-5 15 40];
light4=[-10 -10 40];
light5=[5 -10 40];
light6=[5 20 40];



%%라이트 단위 벡터 만들기
light1 = light1 / norm(light1);
light2 = light2 / norm(light2);
light3 = light3 / norm(light3);
light4 = light4 / norm(light4);
light5 = light5 / norm(light5);
light6 = light6 / norm(light6);

%%주전자
img1 = imread('s_ss1.bmp');
img2 = imread('s_ss2.bmp');
img3 = imread('s_ss3.bmp');
img4 = imread('s_ss4.bmp');
img5 = imread('s_ss5.bmp');
img6 = imread('s_ss6.bmp');


%{
%%주전자
img1 = imread('s_ss1.bmp');
img2 = imread('s_ss2.bmp');
img3 = imread('s_ss3.bmp');
img4 = imread('s_ss4.bmp');
img5 = imread('s_ss5.bmp');
img6 = imread('s_ss6.bmp');

%%구 
img1 = imread('s_gg1.bmp');
img2 = imread('s_gg2.bmp');
img3 = imread('s_gg3.bmp');
img4 = imread('s_gg4.bmp');
img5 = imread('s_gg5.bmp');
img6 = imread('s_gg6.bmp');

%%원뿔
img1 = imread('s_cc1.bmp');
img2 = imread('s_cc2.bmp');
img3 = imread('s_cc3.bmp');
img4 = imread('s_cc4.bmp');
img5 = imread('s_cc5.bmp');
img6 = imread('s_cc6.bmp');

%%교수님 타원
img1 = imread('ellipse1.bmp');
img2 = imread('ellipse2.bmp');
img3 = imread('ellipse3.bmp');
img4 = imread('ellipse4.bmp');
img5 = imread('ellipse5.bmp');
img6 = imread('ellipse6.bmp');
%}

%광원 벡터 6x3
S= [light1; light2 ;light3; light4; light5; light6];

%노멀 벡터가 들어갈 곳 
b=ones(240, 320,3);
b=double(b);

%z를 구하기 위한 p,q 
p=ones(240, 320);
p=double(p);
q=p;
%실제 깊이는 Z에 들어감
Z=ones(240, 320);
Z=double(Z);


for i=1:240
    for j=1:320
        %normal vector를 위한 계산 
        
        %i,j에서의 밝기값 6x1
        E=[img1(i,j) img2(i,j) img3(i,j) img4(i,j) img5(i,j) img6(i,j)];
        E=double(E');
        %의사 코드 풀기
        tb= (inv(S'*S))*S'*E;
        %tb는 최종 3x1이 나온다.
        nbm = norm(tb);
        if( nbm == 0)
            b(i,j,:) = 0; 
        else
            b(i,j,:) = tb / nbm; %tb를 albedo -> norm(tb) 로 나누면 노멀벡터가 구해짐
        end
        
        %Z를 위한 계산 p q를 구함
        tM = [b(i,j,1) b(i,j,2) b(i,j,3)];
        nbm = norm(tM);
        if( nbm == 0)
            tM = [0 0 0];
        else
            tM = tM / nbm; 
        end        
        p(i,j)=tM(1,1);
        q(i,j)=tM(1,2);
    end
end


%pq를 이용 Z 값을 저장함
for i=1:240
    for j=1:320
        Z(i,j) = sum(q(1:i, 1)) + sum(p(i,1:j));
    end
end
Z = Z*-1; %깊이가 역전되어 보임..

figure(1);
hold on;
for i=1:2:240
    for j=1:2:320
        %i,j에 대한 노멀 벡터 뿌림 
        plot3(j+b(i,j,1),i+b(i,j,2),b(i,j,3),'b' );
    end
end
hold off;

figure(2);
mesh(Z);


