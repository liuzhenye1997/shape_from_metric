function result=batch_normalize(vector)
result=vector./sqrt(sum(vector.^2,2));