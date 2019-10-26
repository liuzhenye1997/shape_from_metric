%将向量vector变成单位向量
function result=batch_normalize(vector)
result=vector./sqrt(sum(vector.^2,2));