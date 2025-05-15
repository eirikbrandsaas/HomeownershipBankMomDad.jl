function varnames()
    varnames = Dict(:owneryoung =>"Owner (25-44)",
    :medwealthyoung =>"Median Wealth (25-44)",
    :medwealthold =>"Median Wealth (55-74)",
    :medwealthfyoung =>"Median Family Wealth (25-44)",
    :wealthatpurchase => "Wealth at Purchase (25-44)",
    :transfrateyoung => "Transfer Receipt Rate (25-44)",
    :h2wyoung =>"House Value / Wealth (25-44)",
    :r2incyoung =>"Rent / Income (25-44)",
    :wealthatpurchaseyoung => "Wealth at Purchase (25-44)",
    :first_ownyoung => "Age First Own (25-44)",
    :LTVatpurchaseyoung => "LTV at Purchase (25-44)",
    :tp2wealthpyoung => "Transfers / Parental Wealth (55-74) (\\%)",
    :transferbuyersyoung => "Transfers Around Purchase (25-44)",
    :hand2mouthyoung  => "Hand to Mouth (25-44)",
    :mortg2incyoung => "Mortgage / Income (25-44)",
    :transfsizeyoung => "Transfer Size (55-74)",
    :ownerall => "Owner (25-74)",
    :ownerold => "Owner (55-74)",
    :price => "House Price",
    :α1 => "Supply Elasticity (\$\\alpha_1\$)",
    :α0 => "Supply El. Intercept",
    :supply => "Housing Supply",
    :wealthp_owneryoung => "Own: Parent Wealth",
    :wealthp_renteryoung => "Rent: Parent Wealth",
    :wealthpgradyoung => "Parent Wealth Gradient (avg)",
    :medwealthpgradyoung => "Parent Wealth Gradient (med)",
    :medwealthp_owneryoung => "Own: Parent Wealth",
    :medwealthp_renteryoung => "Rent: Parent Wealth",
    :medconsp2conskyoung => "Consumption Ratio (P/K)",
    :utilk => "Lifetime Utils Kid",
    :utilp => "Lifetime Utils Parent",
    :mortgyoung => "Mortgage (25-44)",
    :owneryoung_prich => "Owner (25-44), Top 33\\%",
    :owneryoung_pmiddle => "Owner (25-44), Middle 33\\%",
    :owneryoung_ppoor => "Owner (25-44), Bottom 33\\%",
    :transfrateyoung_prich => "Transfers (25-44), Parent Top 50\\%",
    :transfrateyoung_ppoor => "Transfers (25-44), Parent Bot 50\\%",
    :transfraterenteryoung => "Transfers (25-44), Renters",
    :transfrateowneryoung => "Transfers (25-44), Owners",
    :smm_obj => "SMM Objective Function",
    )
end

function parnames_s()
    parnames_s = Dict(:β => "beta",
    :γ => "gamma",
    :χ => "chi",
    :η => "eta",
    :rmk => "rmk",
    :ho => "ho",
    :ho2 => "ho2",
    :price => "price" )
end

function parnames()
    parnames = Dict(:β => "Discount Factor (β)",
    :γ => "Risk Aversion (γ)",
    :χ => "Ownership Pref. (χ)",
    :η => "Altruism (η)",
    :rmk => "Mortg. Prem Kid (r_m)",
    :ho => "Size Ratio (h_1/h_r)",
    :ho2 => "Size Ratio (h_2/h_1)",
    :price => "House Price (p)" )
end

function parnames_t()
    parnames_t = Dict(:β => "Discount Factor",
    :γ => "Risk Aversion",
    :χ => "Ownership Pref.",
    :η => "Altruism",
    :rmk => "Mortg. Prem Kid",
    :ho => "Size Ratio",
    :ho2 => "Size Ratio 2",
    :price => "House Price" )
end



function LoadGlobest(nstate,nchoice)
    globest = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/model/est/est_global_nstate$(nstate)_nchoice$(nchoice).csv")))
end

function LoadEstimation(nstate,nchoice,targ_moms,pars)
    globest = LoadGlobest(nstate,nchoice)
    Md = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/model/est/est_global_nstate$(nstate)_nchoice$(nchoice)_data.csv")))

    ## Import data
    # Import Bootstrap sample
    dfb = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/moments/moments_race_age_boot.csv")))
    dfb = dfb[dfb.race .== "All",:]
    dfb = select!(dfb,Not(:race))
    dfb = dfb[:,unique(targ_moms)]

    for var in [:transfsizeyoung, :transferbuyersamntyoung, :transferbuyersamntyoung,:LTVatpurchaseyoung, :transferbuyersyoung]
        globest[isequal.(globest[!,var],NaN),var] .= Inf # Relace NaN's with Inf
    end

    # Import main data
    dat1 = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/moments/moments_race_age.csv")))
    dat = dat1[dat1.race .== "All",:]
    dat[!,:type] .= "Data"
    Md = collect(dat[1,unique(targ_moms)])

    ## Standard Errors
    df_boot, _est_b = BootstrapStandardErrors(globest,targ_moms,pars,dfb)
    println("There were $(length(unique(_est_b.iter))) different parameter vectors that were the best fit for the $(nrow(_est_b)) bootstrap sample. The first number should be 'high'. ")
    est=sort(globest,:smm_obj)[1:1,:]

    filename = joinpath(@__DIR__,"../../data/model/est/estpar_nstate$(nstate)_nchoice$(nchoice)")
    CSV.write(filename*".csv",est[1:1,pars])

    ## Returns:
    # df_boot: Dataframe with estimated parameters and the standard errors
    # globest: Dataframe with results from global estimation.
    # dat: Dataframe with the data used in the estimation.

    return df_boot, globest, dat
end


function BootstrapStandardErrors(est,targ_moms,pars,dfb)
    optb = fill(0,nrow(dfb))
    smmb = fill(-Inf64,nrow(dfb))

    Nb = 100
    println("Bootstrapping standard errors for global estimation. $Nb samples")
    for ib = 1:Nb # nrow(dfb)
        smmtmp = fill(Inf64,nrow(est))
        Mdb = dfb[ib:ib,:]
        for ie = 1:nrow(est)
            smmtmp[ie] = smm_obj(est[ie:ie,targ_moms],Mdb,targ_moms)
        end
        optb[ib] = est.iter[argmin(smmtmp)]
        smmb[ib] = smmtmp[optb[ib]]
        print("$ib,")
    end

    est_b = DataFrame(iter = optb,smm_obj_boot=smmb)
    est_b = leftjoin(est_b,est,on=:iter)[:,unique(union([:iter,:smm_obj,:smm_obj_boot],pars))]

    df_boot = DataFrame(Parameter=pars,
        Value = collect(sort(est,:smm_obj)[1,pars]),
        StdError = collect(combine(est_b, pars .=> std)[1,:]),
        ValueBoot = collect(combine(est_b, pars .=> mean)[1,:]),)
    return df_boot, est_b
end


function PlotGlobalEstimation(dat,est,pars,targ_moms;store=false)
    ## For dictionaries stored in Model()
    _parnames = parnames()
    _varnames = varnames()
    _parnames_s = parnames_s()    

    pylims = DataFrame(
        medwealthyoung = (10,50),
        medwealthold = (100,300),
        owneryoung = (0.4,0.6),
        r2incyoung = (0.18,0.27),
        first_ownyoung = (29,37),
        LTVatpurchaseyoung = (0.4,0.8),
        transfrateyoung = (0.17,0.35),
        transferbuyersyoung = (0,1.0),
        wealthatpurchaseyoung = (25.0,45.0),
        mortg2incyoung = (0.5,2.),
        h2wyoung = (0.0,10.0),
        ownerall = (0.6,0.9),
        mortgyoung = (80,160),
        smm_obj = (0.00005,0.4),
        )

    replace!(est.smm_obj,NaN=>Inf)
    nqs = 20
    dfq, dfp = calculate_global_est(nqs,est,pars,[targ_moms;:smm_obj])
    plots = [plot() for npar = 1:length(pars), nmom = 1:length(targ_moms) ]
    for (ipar,par) in enumerate(pars)
        for (imom,mom) in enumerate(targ_moms)

            plot(dfq[!,par],[dfp[!,Symbol(par,"lo_",mom)] dfp[!,Symbol(par,"hi_",mom)]],color=:gray,linestyle=:dash)
            plot!(dfq[!,par],[dfp[!,Symbol(par,"med_",mom)]],color=:black,)
            if mom != :smm_obj
                hline!([dat[!,mom]],color="orange",style=:dot)
            end

            plot!(legend=false)
            ptemp = plot!(title=string(par)*" "*string(mom))
            plots[ipar,imom] = ptemp
        end
    end
    plots_obj = [plot() for npar = 1:length(pars) ]
    for (ipar,par) in enumerate(pars)
        mom = :smm_obj

            plot(dfq[!,par],[dfp[!,Symbol(par,"lo_",mom)] dfp[!,Symbol(par,"hi_",mom)]],color=:gray,linestyle=:dash)
            plot!(dfq[!,par],[dfp[!,Symbol(par,"med_",mom)]],color=:black,)
            plot!(dfq[!,par],[dfp[!,Symbol(par,"mnm_",mom)]],color=:red,)
            ptemp  = plot!(legend=false,ylims=pylims[1,mom],title=_varnames[mom],yscale=:log)            
            plots_obj[ipar] = ptemp

    end


    p1 =deepcopy(plots)
    p2 =deepcopy(plots)
    size1 = (700,400).*0.35

    ### Plots each plot alone
    for (ipar,par) in enumerate(pars)
        for (imom,mom) in enumerate(targ_moms)
        p1[ipar,imom] = plot(plots[ipar,imom],ylims=pylims[1,mom],title="",xlabel=_parnames[par])
        p2tmp = deepcopy(p1[ipar,imom])
        p2[ipar,imom] = plot!(p2tmp,title=_varnames[mom],xlabel="")
        punitemp = deepcopy(p1[ipar,imom])
        puni = deepcopy(plot!(punitemp,title="",size=size1))
        if store == true
            Plots.pdf(puni,"tabfig/est/identification/"*_parnames_s[par]*"_"*String(mom))
        end
        end
    end

    ## Plot parameter - objective function
    Plots.pdf(plot(plots_obj...),"tabfig/est/identification/smmobj")
    ### Plots each moment for all parameter
    for (imom,mom) in enumerate(targ_moms)
    ptmp = deepcopy(plot(p1[:,imom]...))
    if store == true
        Plots.pdf(ptmp,"tabfig/est/identification/"*String(mom))
    end
    end

    ### Plots the effect of each parameter
    for (ipar,par) in enumerate(pars)
        ptmp = deepcopy(plot(p2[ipar,:]...,plots_obj[ipar],titlefontsize=9,layout=(2,3),size=(800,300).*0.8))
        display(ptmp)
        if store == true
            Plots.pdf(ptmp,"tabfig/est/identification/"*_parnames_s[pars[ipar]])
        end
    end
end


function filename_wdate(filename::String)
    filename*"_"*@sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
end


function output_esttable(dfb,pars,moms,M::Model;store=false)

    parlatex = Dict(
        :β => L"\beta",
        :χ => L"\chi",
        :η => L"\eta",
        :rmk => L"r^m",
        :ho => L"\frac{h_o}{h_r}",
        :price => L"p"
    )

    
    Model(ModPar(NumPar()),NumPar());

    _parnames_t = parnames_t()  # Dictionary for pretty printing
    dfp = DataFrame(
        Parameter = pars,
        Parameter_str = [_parnames_t[par]*" ("*parlatex[par]*")" for par in pars],
        )

    dfp = leftjoin(dfp,dfb,on=:Parameter)
    select!(dfp,[:Parameter_str,:Value,:StdError])
    rename!(dfp,:StdError=>"Standard Error")
    display(dfp)
    par_latex = latexify(dfp; env=:table, latex=false,fmt="%1.3f",escape_underscore=true,booktabs = true)
    par_latex = replace(par_latex,"ccc" => "lrr") # Change Alignement
    par_latex = replace(par_latex,"Mortg. Prem Kid" => "Mortg. Prem.") # Change Alignement
    par_latex = replace(par_latex,"Parameter_str" => " \\textbf{Parameter}") # Change column headers
    par_latex = replace(par_latex,"}\n\\toprule" => "}\n\\toprule   ") # Change column headers
    par_latex = replace(par_latex,"\\end{tabular}\n" => "")
    par_latex = replace(par_latex,"\\bottomrule" => "\\cmidrule(lr){1-3} ")

    dfm = copy(M.moms)
    dfm.smm_obj .= NaN
    select!(dfm,union([:smm_obj,:type],moms))
    
    
    out = DataFrame(Moment = names(dfm[!,moms]),
            Data = Matrix(dfm[dfm.type .== "Data",moms])[1,:],
            Model = Matrix(dfm[dfm.type .== "Model",moms])[1,:]
            )

    _varnames = varnames()
    out[!,:Moment] = [_varnames[Symbol(var)] for var in out[!,:Moment]]
    smmfactor = 100 # Increase scale
    push!(out,["\\midrule \n Sum Squared Distances (\$ \\times $smmfactor\$)", M.moms[M.moms.type.=="Data",:smm_obj][1], M.moms[M.moms.type.=="Model",:smm_obj][1]*smmfactor])
    display(out)
    out_latex = latexify(out,; env=:table, latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,
            head= ["\\textbf{Moment}", "Data", "Model"])
    out_latex = replace(out_latex,"ccc" => "lrr")
    out_latex = replace(out_latex,"NaN" => "")
    out_latex = replace(out_latex,"}\n\\toprule" => "}\n ") # Remove toprule
    out_latex = replace(out_latex,"\\begin{tabular}{lrr}" => "")

    filename = "tabfig/est/esttable"
    filename2 = filename_wdate(filename)
    table = par_latex*out_latex
    println(table)

    display(string(table))

    if store == true
            write(filename*".tex", table)
            write(filename2*".tex", table)
    end
end

## Extra Moments
function output_extramoments(Mb::Model,mom_extra::Vector{Symbol},;store=false)
    
    dfe = Mb.moms[:,[mom_extra;:type]]

    try
        dfe.tp2wealthpyoung *= 100
    catch
    end

    out = DataFrame(Moment = names(dfe[!,mom_extra]),
            Data = Matrix(dfe[dfe.type .== "Data",mom_extra])[1,:],
            Model = Matrix(dfe[dfe.type .== "Model",mom_extra])[1,:]
            )

    _varnames = varnames()
    out[!,:Moment] = [_varnames[Symbol(var)] for var in out[!,:Moment]]

    display(out)
    out_latex = latexify(out,; env=:table, latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,
    head= ["\\multicolumn{1}{l}{Moment}", "\\multicolumn{1}{c}{Data}", "\\multicolumn{1}{c}{Model}"])
    out_latex = replace(out_latex,"-Inf" => "0.07")
    out_latex = replace(out_latex,"ccc" => "lrr")

    ## Code that prints some distributional moments to show that the model matches these patterns okay
    wa_label = "Wealth percentiles, Age 35 (10,25,50,75,90)" 
    wa_dat = string(round.(Int,collect(Mb.moms[Mb.moms.type.=="Data","wealth_age3536_p" .* string.([10, 25, 50, 75, 90]) .* "young"][1,:])))
    wa_mod = string(round.(Int,collect(Mb.moms[Mb.moms.type.=="Model","wealth_age3536_p" .* string.([10, 25, 50, 75, 90]) .* "young"][1,:])))
    wa_row  = wa_label * " & " * wa_dat * " & "* wa_mod * "  \\\\ "

    oi_label = "Ownership, by income tertile, Age 35"
    oi_dat = string(round.(collect(Mb.moms[Mb.moms.type.=="Data","owner_age3536_inc" .* string.([1,2,3]) .* "young"][1,:]);digits=2))
    oi_mod = string(round.(collect(Mb.moms[Mb.moms.type.=="Model","owner_age3536_inc" .* string.([1,2,3]) .* "young"][1,:]);digits=2))
    oi_row  = oi_label * " & " * oi_dat * " & "* oi_mod * "  \\\\ "

    out_latex = replace(out_latex,"\\bottomrule" => wa_row * "\n\\bottomrule")
    out_latex = replace(out_latex,"\\bottomrule" => oi_row * "\n\\bottomrule")
    out_latex = replace(out_latex,"[" => "")
    out_latex = replace(out_latex,"]" => "")
    
    println(out_latex)
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/est/extramoments_"*append_info*".tex"
    filename2 = "tabfig/est/extramoments_mostrecent.tex"
    if store == true
        write(filename, out_latex)
        write(filename2, out_latex)
        println("Output stored as ",filename)
    end
end


"""
    function calculate_global_est(data,simdata,parameters,moments)

divides the parameters into quartiles for identification plots
"""
function calculate_global_est(nqs,dfs,parameters::Vector{Symbol},moments::Vector{Symbol})
    step = 1/(nqs)

    dfq = DataFrame() # Quantile levels

    for ipar in 1:length(parameters)
        symb = parameters[ipar]
        dfq[!,symb] = quantile(dfs[!,symb],collect(step:step:1))
    end

    insert!.(eachcol(dfq), 1, [-Inf]) # Adds a row of -Infs to have lower bounds for quartiles
    dfp = DataFrame() # Plotting data
    lo = zeros(nqs)
    hi = zeros(nqs)
    med = zeros(nqs)
    avg = zeros(nqs)
    mnm = zeros(nqs)
    dfq2 = copy(dfq)
    for mom in moments
        for par in parameters
            for iq = 1:nqs
                iqconv = iq + 1 # Because you need to start at -inf, so the vectors are nq+1 long
                subgroup = dfs[!,mom][(dfs[!,par].>dfq[!,par][iqconv-1]) .&  (dfs[!,par] .< dfq[!,par][iqconv]) ,:][:,1]
                lo[iq] = percentile(skipmissing(subgroup),25)
                hi[iq] = percentile(skipmissing(subgroup),75)
                med[iq] = median(skipmissing(subgroup))
                avg[iq] = mean(skipmissing(subgroup))
                mnm[iq] = minimum(skipmissing(subgroup))
                dfq2[!,par][iqconv] = mean(dfs[!,par][(dfs[!,par].>dfq[!,par][iqconv-1]) .&  (dfs[!,par] .< dfq[!,par][iqconv]) ,:][:,1])
            end
            dfp[!,Symbol(par,"lo_",mom)] = copy(lo)
            dfp[!,Symbol(par,"hi_",mom)] = copy(hi)
            dfp[!,Symbol(par,"avg_",mom)] = copy(avg)
            dfp[!,Symbol(par,"med_",mom)] = copy(med)
            dfp[!,Symbol(par,"mnm_",mom)] = copy(mnm)

        end
    end
    delete!(dfq,1)
    delete!(dfq2,1)

    return dfq2,dfp

end