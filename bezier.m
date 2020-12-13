function pts = bezier(p0, p1, p0c, p1c, qtd)

i = 1;
pts = zeros(qtd, 2);
for t=0:1/(qtd-1):1
    b0 = (1.0 - t) * (1.0 - t) * (1.0 - t);
	b0c = 3.0 * (1.0 - t) * (1.0 - t) * t;
	b1c = 3 * (1.0 - t) * t * t;
	b1 = t * t * t;
   
	pts(i,:) = p0 * b0 + p0c * b0c + p1c * b1c + p1 * b1;
    i = i + 1;
end
   
end

