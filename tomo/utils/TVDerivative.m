function Vst = TVDerivative(img)

    fxy = padarray(img, [1,1], 0, 'both');
    fxnegy = circshift(fxy, [-1, 0]);
    fxposy = circshift(fxy, [1, 0]);
    fnegxy = circshift(fxy, [0, -1]);
    fposxy = circshift(fxy, [0, 1]);
    fposxnegy = circshift(fxy, [-1 1]);
    fnegxposy = circshift(fxy, [1 -1]);
    eps = realmin;

    vst1 = (2*(fxy - fnegxy) + 2*(fxy - fxnegy))./sqrt(eps + (fxy - fnegxy).^2 ...
    + (fxy - fxnegy).^2);

    vst2 = (2*(fposxy - fxy))./sqrt(eps + (fposxy - fxy).^2 ...
    + (fposxy - fposxnegy).^2);

    vst3 = (2*(fxposy - fxy))./sqrt(eps + (fxposy - fxy).^2 ...
    + (fxposy - fnegxposy).^2);

    vst = vst1 - vst2 - vst3;
    Vst = vst(2:end-1,2:end-1);

end 
