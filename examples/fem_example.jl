topology = [1 2 3  5  6  7;
            2 3 4  6  7  8;
            6 7 8 10 11 12;
            5 6 7  9 10 11]
Le = zeros(Int64,size(topology,1),maximum(topology))
element = 3
element_topology = topology[:,element]
for (idx,element) in enumerate(element_topology)
    Le[idx,element] = 1
end
Le
p = collect(1:maximum(topology))
pe = Le*p

Ke = rand(4,4)
K1 = Le'*Ke*Le
K2 = zeros(size(K1))
K2[element_topology,element_topology] = Ke
K1 - K2
Le*K2*Le'

Le*Le'
M = Le'*Le


T = collect(1:maximum(topology))'
Te = T*Le'

Te*Le

T*(Le'*Le)


Te*Le


y = rand(size(T,2))

Te*Le*y
