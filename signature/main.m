im=imread("sig.jpg");
% imshow(im)
im=im(:,:,1); b=im>60;
im(b)=255;

% imwrite(im,"sig2.png",Alpha=uint8(b))
imwrite(im,"sig2.png")
%%
system("convert sig2.png -transparent white sig3.png")

%%
[im, map, alpha] =imread("sig.jpg");