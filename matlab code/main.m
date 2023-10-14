clear;
clc;

global hght;
global wdth;
global wsp;
global hsp;


l_eye=[202 233 235 238 236 282 244 207 203 241 240 234 239 237 283 245 360 359 138 137 133 134 124 126 122 131 112 116 397 148 149 430 136 135 132 123 125 121 130 129 127 128 248 249 250 451 247 114 36 37 38 350 109 105 107 98 102 103 100 119 117 354 355 358 454 113 115 110 111 106 108 99 104 101 120 118 356 357 402 403 452];
r_eye=[641 679 633 635 632 630 599 642 680 634 636 631 637 638 600 644 834 647 525 528 519 523 521 531 530 534 535 754 604 837 835 797 796 646 645 524 526 527 518 522 520 529 532 533 536 547 546 791 755 753 750 749 514 516 497 500 499 495 504 502 506 745 484 483 482 511 513 752 751 515 517 498 501 496 505 503 508 507 512 510 509];
nose=[219 365 366 475 199 197 24 5 6 26 594 596 850 761 760 616 217 220 261 260 433 303 304 305 306 477 478 388 198 279 281 261 260 259 258 262 90 91 92 94 479 419 387 296 263 264 89 93 339 334 336 337 11 83 7 8 87 86 84 85 88 9 10 595 852 702 701 700 821 781 813 854 491 782 853 490 852 488 15 733 732 489 485 486 656 655 660 661 735 730 487 658 659 657 693 676 678 614 617 379 30 32 434 380 31 33 34 35 774 440 825 442 441 21 20 19 457 459 481 826 827 461 841 460 400 398 404 232 18 17 28 29 458 480 629 798 840 225 224 223 280 278 222 221 215 216 218 364 399 201 200 25 22 23 27 597 598 793 759 615 613 792 794 622 621 620 612 618 619 677 675];
lips=[267 266 265 298 269 268 297 179 295 272 271 270 301 186 189 273 302 300 299 173 178 188 307 319 317 321 185 182 183 187 320 318 323 322 325 324 167 175 184 170 171 172 315 316 313 314 311 312 1 55 53 59 57 61 65 63 69 67 71 76 2 56 54 60 58 62 66 64 70 68 72 78 565 572 581 567 568 569 711 712 709 710 707 708 570 575 582 579 580 584 716 714 719 718 721 720 576 692 583 586 585 703 715 713 717 697 696 670 699 695 669 668 667 664 663 666 665 694 664 663 662 698];
edge=[95 96 97 294 392 150 151 155 333 332 331 401 255 257 256 412 411 405 406 408 427 425 426 417 416 415 424 423 429 444 809 810 811 818 817 819 802 800 799 805 806 653 808 654 816 652 814 778 795 727 777 762 728 729 763 549 553 548 786 691 494 493 492 414 422 420 384 383 367 368 369 154 153 152 293 292 276 274 275 242 243 288 764 552 551 550 690 689 673 671 672 639 640 685];
B=[l_eye r_eye nose lips edge];

ii=30;
jj=65;

File=['F:\2023.05.04\1\calculate\01\',num2str(ii),'.bmp'];

figure;imshow(File);hold on;
set(gcf,'Name','calculation area')
axis ij
rect=getrect;
rect=floor(rect);

startx=rect(1);
starty=rect(2);
endx=rect(1)+rect(3);
endy=rect(2)+rect(4);

startpt=[starty,startx];
endpt=[endy,endx];
    
N=endpt-startpt+1;

hght=35;
wdth=35;  

wsp=35;
hsp=35;

xBD=[startx startx endx endx];
yBD=[starty endy endy starty];
plot(xBD,yBD,'r*');hold on;

load('F:\district\triMesh.mat');
A=load(['F:\2023.05.04\1\calculate\01\landmark-',num2str(ii),'.txt']);
region_num=size(triMesh,1);

File1=['F:\2023.05.04\1\calculate\01\',num2str(ii),'.bmp'];
File2=['F:\2023.05.04\1\calculate\02\',num2str(ii),'.bmp'];

beimage=double(imread(File1)); 
afimage=double(imread(File2));

[image_x,image_y]=size(beimage);

point=zeros(image_x,image_y);

for xx=1:wsp:N(2)
    for yy=1:hsp:N(1)
        point(startpt(1)+yy-1,startpt(2)+xx-1)=1;
    end
end

figure;
imshow(uint8(beimage));
hold on;

BW=zeros(image_x,image_y);
for i=1:region_num

    if ismember(i,B)
        continue
    end
    im=beimage;
    f=triMesh(i,:);

    c=[A(f(1)+1,2) A(f(2)+1,2) A(f(3)+1,2)];
    r=[A(f(1)+1,3) A(f(2)+1,3) A(f(3)+1,3)];

       BW=BW+double(roipoly(im,c,r));
end

point=point.*BW;
[p(:,1),p(:,2)]=find(point);


plot(p(:,2),p(:,1),'+','Color','w')
    
    [acu,acv,coef,pointx,pointy]=DIC_cxs(p,beimage,afimage);
    er=find(coef>-0.7);
    acu(er)=[];
    acv(er)=[];
    pointx(er)=[];
    pointy(er)=[];

    TimePoint=[pointy(:) pointx(:)];
    Ry=pointy+acu(:);
    Rx=pointx+acv(:);
    figure; imshow(File1);hold on
    plot(pointy(:),pointx(:),'r+');hold off;
    
    figure;imshow(File2);hold on;
    plot(Ry, Rx,'r+');hold off;
    PtOI=stereo([pointy(:) pointx(:)],[Ry  Rx]);
    [m,n]=size(pointx);

    p=[pointx,pointy];
    p1=[Rx,Ry];
    
    p_n=length(pointx);
    
l=jj-ii;
lcoef=zeros(p_n,l);
rcoef=zeros(p_n,l);
smile_left=zeros(p_n,3,l);
px=zeros(p_n,l);
py=zeros(p_n,l);

    
parfor j=ii+1:jj
    
        p_=p;
        p1_=p1;
        pointx_=pointx;
        pointy_=pointy;
        PtOI_=PtOI;
        acu_=acu;
        acv_=acv;
        
        FileL1=imread(['F:\2023.05.04\1\calculate\01\',num2str(j),'.bmp']);

        FileR1=imread(['F:\2023.05.04\1\calculate\02\',num2str(j),'.bmp']);

        IMGL1=double(FileL1);
        IMGR1=double(FileR1);
        [acuL,acvL,coefL,pointxL,pointyL]=DIC_cxs2(p_,beimage,IMGL1);
        erL=find(coefL>-0.8);

        [acuR0,acvR0,coefR,pointxR,pointyR]=DIC_cxs2(p1_,afimage,IMGR1);
        erR=find(coefR>-0.8);

        acuR=acuR0+acu_;
        acvR=acvR0+acv_;

        Lx1=pointx_(:)+acvL(:);
        Ly1=pointy_(:)+acuL(:);
        Rx1=pointx_(:)+acvR(:);
        Ry1=pointy_(:)+acuR(:);
        output=stereo([Ly1 Lx1],[Ry1 Rx1]);

        
        stereotrans=output-PtOI_;
        lcoef(:,j-ii)=coefL';
        rcoef(:,j-ii)=coefR';
        smile_left(:,:,j-ii)=stereotrans;
        px(:,j-ii)=pointx_;
        py(:,j-ii)=pointy_;
        

end

strain;

clear uu;
clear vv;
acu=smile_left(:,1,2);
acv=smile_left(:,2,2);
pointx=px;
pointy=py;
    acuu=acu;
    acvv=acv;
    pointxx=pointx(:,2);
    pointyy=pointy(:,2);
    for i=1:1285
        uu(pointxx(i),pointyy(i))=acu(i);
        vv(pointxx(i),pointyy(i))=acv(i);
    end
    
    u_x=sum(uu,2);
    u_y=sum(uu);
    v_x=sum(vv,2);
    v_y=sum(vv);
    uu(find(~u_x),:)=[];
    uu(:,find(~u_y))=[];
    vv(find(~v_x),:)=[];
    vv(:,find(~v_y))=[];
    
    mu=mean(uu(:));
    su=std(uu(:));
    Tu=[mu-3*su mu+3*su];
    uu(find(uu<Tu(1)))=nan;
    uu(find(uu>Tu(2)))=nan;
    
    mv=mean(vv(:));
    sv=std(vv(:));
    Tv=[mv-4*sv mv+4*sv];
    vv(find(vv<Tv(1)))=nan;
    vv(find(vv>Tv(2)))=nan;

