#######################################################
#Inputs
#######################################################
#depth is the depth of the tree
#product(s)=ds(s)*ds(s+1)*...*ds(depth),s=1,...,depth
#index: index of a scenario 
#######################################################
#######################################################
#Output
#######################################################
#etas(1),...,etas(depth) is the scenario with index
#index.
#######################################################

function index_to_scenario(index,product,depth)

    etas = zeros(Int,depth,1)
    for i=1:(depth-1)
        etas[i] = ceil(index/product[i+1])
        index = index-(etas[i]-1)*product[i+1]
    end
    etas[depth] = index

    return etas
end