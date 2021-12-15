function yq=qqlog(quat)
yq= quaternion(log(normq(quat))./2, axis(quat) .* angle(quat));
end