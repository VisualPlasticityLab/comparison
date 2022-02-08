function ODIfinal=convertODI(ODIi)

ODIfinal=ODIi;
for i=1:numel(ODIi)
    ODIfinal(i)= sqrt(abs(ODIi(i)))*sign(ODIi(i));
end
