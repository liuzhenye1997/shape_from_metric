function result=Dt(iteration)
if iteration<10
    result=2;
elseif iteration<20
    result=1;
elseif iteration<50
    result=1;
else
    result=0.2;
end