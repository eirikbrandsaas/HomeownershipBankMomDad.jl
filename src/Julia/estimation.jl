function EstPar(nstate,nchoice)
    np = NumPar(nxk=nstate,nxp=nstate,nbp=nstate,
        ncp_choice=nchoice,nck_choice=nchoice,ntp_choice=nchoice,
        nxphat = nstate, ntp_choice2 = nchoice, ncp_choice2 = nchoice)
    mp = ModPar(np)

    return np,mp
end

## Global
function GlobalSearch(par_lims,nstate,nchoice,targ_moms::Vector{Symbol};nsim = 1,start=1)
    @assert start<=nsim

    npar = ncol(par_lims)
    sobolseq = SobolSeq(npar)

    np,mp = EstPar(nstate,nchoice)
    M = Model(mp,np)

    @assert np.nhr==1
    @assert np.nho==1

    pars = DataFrame(iter=[0])
    for col in (names(par_lims))
        pars[!,col] = [0.0]
    end
    if start > 1
        for i = 1:start-1
            sobolveq =next!(sobolseq)
        end
    end

    for i = start:nsim
        M.numdet.cnvrgd = false
        sobolveq =next!(sobolseq)
        isobol = 0
        for col in (names(par_lims))
            isobol += 1
            pars[!,col] .=  (par_lims[1,col][2] - par_lims[1,col][1])*sobolveq[isobol] + par_lims[1,col][1]
        end
        pars[!,:iter] .= i

        mp = ModPar(np,
            η = pars.η[1],
            χ = pars.χ[1],
            ho = pars.ho[1]
        )

        println(@sprintf("Iteration %4i",i))
        display(pars)
        M = Model(mp,np)
        @time solvemodel2nd!(M,np,mp)

        smmobj = smm_obj(M.moms[M.moms.type .== "Model",:],M.moms[M.moms.type .== "Data",:],targ_moms)
        df = DataFrame(iter = i,
                        smm_obj = smmobj,
                        converged = M.numdet.cnvrgd,
                        distance = M.numdet.dist,
                        simdistance = M.numdet.simdist,
                        cyclelength = M.numdet.cycle_length
                        )
        df = leftjoin(pars,df,on=:iter)
        M.moms[!,:iter] .= i # Add so you can joine the two files
        output=leftjoin(df,M.moms[M.moms.type .== "Model",:],on = :iter)
        path = joinpath(@__DIR__,"../../data/model/est/")
        filename = path*"est_global_nstate$(nstate)_nchoice$(nchoice)"

        if i == 1
            CSV.write(filename*".csv",output)
            CSV.write(filename*"_data.csv",M.moms[M.moms.type .== "Data",:])
        else
            CSV.write(filename*".csv",output,append=true)
        end


    end
end

## Function to load best parameters
function BestPars(nstate,nchoice)

    par = (
            ho = 2.65,
            η = 0.2,
            χ = 1.145,
        )
    try
        est = DataFrame(CSV.File("data/model/est/estpar_nstate$(nstate)_nchoice$(nchoice).csv"))
        par = (
            ho = est.ho[1],
            η = est.η[1],
            χ = est.χ[1],
        )
    catch
        @warn "No valid parameter combination given, using benchmark values!"
    end

    return par
end



function shrink_searchspace(firstiter,lastiter,nstate,nchoice,par_lims)
    ilo = firstiter
    ihi = lastiter
    globest = LoadGlobest(nstate,nchoice)
    @assert ilo >= 1
    @assert ihi <= maximum(globest.iter)
    globest = globest[ilo .<= globest.iter .<= ihi,:]
    sort!(globest,:smm_obj)


    perc = 5
    imax = floor(Int,perc/100*nrow(globest))
    _df = globest[1:imax,:] # keep the `perc`th best objectives

    par_lims_n = DataFrame(
        η = (minimum(_df[!,:η]),maximum(_df[!,:η])),
        χ = (minimum(_df[!,:χ]),maximum(_df[!,:χ])),
        ho = (minimum(_df[!,:ho]),maximum(_df[!,:ho])))

    println("Shrinks the search space")
    ηr = ((par_lims_n.η[1][2] - par_lims_n.η[1][1])/(par_lims.η[1][2]-par_lims.η[1][1]))
    hor = ((par_lims_n.ho[1][2] - par_lims_n.ho[1][1])/(par_lims.ho[1][2]-par_lims.ho[1][1]))
    χr = ((par_lims_n.χ[1][2] - par_lims_n.χ[1][1])/(par_lims.χ[1][2]-par_lims.χ[1][1])) 
    @printf("η : (%1.3f,%1.3f) vs (%1.3f,%1.3f). Size ratio of new over old %1.3f \n",par_lims_n.η[1][1],par_lims_n.η[1][2],par_lims.η[1][1],par_lims.η[1][2], ηr)
    @printf("ho: (%1.3f,%1.3f) vs (%1.3f,%1.3f). Size ratio of new over old %1.3f \n",par_lims_n.ho[1][1],par_lims_n.ho[1][2],par_lims.ho[1][1],par_lims.ho[1][2], hor)
    @printf("χ : (%1.3f,%1.3f) vs (%1.3f,%1.3f). Size ratio of new over old %1.3f \n",par_lims_n.χ[1][1],par_lims_n.χ[1][2],par_lims.χ[1][1],par_lims.χ[1][2],χr)
    @printf("Total reduction in search area: %1.3f percent  \n",(1 - ηr*hor*χr)*100)

    return par_lims_n
end



function par_lims_ini()
    par_lims = DataFrame(
            η = (0.11,0.15),
            χ = (1.2,1.45),
            ho = (2.45,2.9),
        )
end