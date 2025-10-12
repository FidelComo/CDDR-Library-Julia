#######################################################
#Inputs.
#######################################################
#depth: depth of the tree
#product(s)=ds(s)*ds(s+1)*...*ds(depth),s=1,...,depth
#etas: scenario etas(1),etas(2),...,etas(depth)
#######################################################
#Outputs
#######################################################
#index is the index of scenario etas in tree given by 
#product
#######################################################

function scenario_to_index(etas,product,depth)

    index = etas[depth]+(etas[1:depth-1].-1)'*product[2:depth]
    
    return index
end