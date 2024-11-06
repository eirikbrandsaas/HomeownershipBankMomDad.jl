include("../src/HomeownershipBankMomDad.jl")
using Test

##
@testset "Test that solution didn't change and code runs" begin
    nstate = 10 # Grid points in state
    nchoice = 30 # Grid points in choice
    @testset "Model solution" begin 

        mp, np = benchpar(nstate,nchoice)
        np.Ndyn = 100
        np.Ngen = 3
        Mb = Model(mp,np)
        solvemodel2nd!(Mb,np,mp)

        @test Mb.moms[Mb.moms.type.=="Model",:owneryoung][1] == 0.413
        @test Mb.simpan.xp[15] == 78.40585097145825
        @test Mb.g.gp.t[1,10,3,2,2,1,14,1]  == 357.96712339010395
    end

    ##
    @testset "Estimation" begin
        Nest = 1 # How many iterations in each step of the estimation
        pars = [:η]
        targ_moms = [:medwealthyoung]
        par_lims = par_lims_ini()
        GlobalSearch(par_lims, nstate,nchoice,targ_moms;nsim = 1)

        df_boot, globest,dat = LoadEstimation(nstate,nchoice,targ_moms,pars);

        @test globest.smm_obj[1]== 0.7092005784698485 
    end
end
