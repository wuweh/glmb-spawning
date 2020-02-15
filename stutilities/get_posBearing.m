function [ beta, vmag ] = get_posBearing( pos1,pos2,dt )

d = pos2-pos1;
dmag = norm(d);
vmag = dmag/dt;
beta = atan2(d(2)/dmag,d(1)/dmag); betaDeg = rad2deg(beta); %#ok<NASGU> % for debugging
while beta > 2*pi; beta = beta-2*pi; end
while beta < 0;    beta = beta+2*pi; end

end

