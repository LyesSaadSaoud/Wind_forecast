function yy=sqrtq(x)

yy=sqrt(abs(x)) .* exp(axis(x) .* angle(x) ./ 2);
end