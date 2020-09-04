function [ data ] = interpolateMap( map ,RA, DECL, dirs )
%
% dirs is a Np x 2 matrix containing
%    ra1   decl1
%    ra2   decl2
%    ...   .....
%    raNp  declNp

    if length(RA)<2 || length(DECL)<2
        error('Need full RA and DECL vectors');
    end


    %Np=size(dirs,1);
    Nra  =length(unique(RA));
    Ndecl=length(unique(DECL));
    sRA   =reshape(RA    ,Ndecl,Nra);
    sDECL =reshape(DECL  ,Ndecl,Nra);
    smap  =reshape(map   ,Ndecl,Nra);
    data=interp2(sRA,sDECL,smap,dirs(:,1),dirs(:,2));

end

