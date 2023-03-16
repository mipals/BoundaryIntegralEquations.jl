#===========================================================================================
                                Scaling of model
===========================================================================================#
nDOF = n
nDOF = 154_000
#* A single BEM matrix
2(nDOF)^2*2*8/(2^30)
#* Full dense system (100 times a single BEM Matrix)
(10nDOF)^2*2*8/(2^30)
#* 6 BEM Systems, 2 inverse products (Gv^{-1}Hv, Gh^{-1}Hv) and 1 collection
(9)*(nDOF)^2*2*8/(2^30) #! Computing inverse would also kill performance...

Base.summarysize(LGM)/(2^30)
#===========================================================================================
                                Memory tracking
===========================================================================================#
Ga = FMMGOperator(mesh,kₐ;n_gauss=1,tol=1e-6,offset=0.2,nearfield=true,depth=1)
@time Ga*pa;
Ha = FMMFOperator(mesh,kₐ;n_gauss=3,tol=1e-6,offset=0.2,nearfield=true,depth=1) + 0.5I
@time Ha*pa;

mem_matrix = Base.summarysize(LGM)/(2^20)
mem_multiply = @allocated begin
    LGM*vn0
end
@time LGM*vn0;
@time LGM.Ga*vn0;
@time LGM.Ha*vn0;

@time LGM.inner*vn0;

Base.summarysize(LGM.Gv)/(2^20)
Base.summarysize(LGM.Hv)/(2^20)
Base.summarysize(LGM.Gh)/(2^20)
Base.summarysize(LGM.Hh)/(2^20)
Base.summarysize(LGM.Ga)/(2^20)
Base.summarysize(LGM.Ha)/(2^20)
Base.summarysize(LGM.Nd)/(2^20)
# Note that Dc requires less memory as we use a CSC format (and Dc has less columns)
Base.summarysize(LGM.Dr)/(2^20)
Base.summarysize(LGM.Dc)/(2^20)
Base.summarysize(LGM.inner)/(2^20)

(Base.summarysize(LGM.Gv) + Base.summarysize(LGM.Hv) +
 Base.summarysize(LGM.Gh) + Base.summarysize(LGM.Hh) +
 Base.summarysize(LGM.Ga) + Base.summarysize(LGM.Ha) +
 Base.summarysize(LGM.Dr) + Base.summarysize(LGM.Dc) +
 Base.summarysize(LGM.Nd))/(2^20)
Base.summarysize(LGM)/(2^20) # Why is this less?

LGM.Dr[1:n,0n+1:1n] - LGM.Dc[0n+1:1n,1:n]
LGM.Dr[1:n,1n+1:2n] - LGM.Dc[1n+1:2n,1:n]
LGM.Dr[1:n,2n+1:3n] - LGM.Dc[2n+1:3n,1:n]
