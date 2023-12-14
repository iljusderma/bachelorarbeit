using BenchmarkTools, Random

function BitInject(state)
    L = length(state)
    I = BitArray(undef, L)
    I[1] = 1
    if state .& I == BitArray(undef, L)
        state[1] = 1
    end
 end

 function Inject(state)
    L = length(state)
    if state[1] == 0
        state[1] == 1
    end
 end

 @btime BitInject(bitrand(1e10))
 @btime Inject(zeros(1e10))