

function y=hh(x)
if norm(x)~=0
    y=x/norm(x);
else if norm(x)==0
y=0;
  end
end
end