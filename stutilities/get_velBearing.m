function [ beta ] = get_velBearing( vx,vy )

v = norm([vx,vy]);
beta = atan2(vy/v,vx/v);
while beta > 2*pi; beta = beta-2*pi; end
while beta < 0;    beta = beta+2*pi; end

end

