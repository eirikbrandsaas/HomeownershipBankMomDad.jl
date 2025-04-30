"""
    function altrcontrhousing(Mbench::Model)

Performs the main quantiative exercise by solving the model without altruism
"""
function altrcontrhousing(mpin,npin)

    ## initialize
    mp = deepcopy(mpin)
    np = deepcopy(npin)
    
    ## Solve model without altruism

    mp.η   = 0.0
    mp.τt   = 1.0
    mp.τb   = 1.0
    Mn = Model(mp,np)
    solvemodel2nd!(Mn,np,mp)

    return Mn
end

function output_noalt(Mb::Model,Mn::Model,targ_moms::Vector{Symbol},nontarg_moms::Vector{Symbol};store=false)
    @assert Mb.moms[Mb.moms.type.=="Model",:transfrateyoung][1] > 0.0
    @assert Mn.moms[Mn.moms.type.=="Model",:transfrateyoung][1] == 0.0

    out_moms =  union(targ_moms,nontarg_moms)
    # out_moms = [:owneryoung, :medwealthyoung, :medwealthold, :wealthatpurchase, :transfrate, :r2incyoung]
    df = DataFrame(title=["Data", "Altruism", "Without"])
    temp = vcat(Mb.moms[:,out_moms],Mn.moms[Mn.moms.type .== "Model",out_moms])
    df = hcat(df,temp)

    ##

    out = DataFrame(varname = names(df[:,2:end]),
        Data = vec(Matrix(df[df.title .=="Data",out_moms])),
        Altruism = vec(Matrix(df[df.title .=="Altruism",out_moms])),
        Without = vec(Matrix(df[df.title .=="Without",out_moms])),
        )

    allowmissing!(out,[:Data,:Altruism,:Without])
    _varnames = varnames()
    out[!,:varname] = ["\\;"*_varnames[Symbol(var)] for var in out[!,:varname]]
    insert!.(eachcol(out), length(targ_moms)+1, ["\\textit{Non-Targeted Moments}", missing, missing,missing])
    insert!.(eachcol(out), 1, ["\\textit{Targeted Moments}", missing, missing,missing])
    display(out)
    df_latex = latexify(out; env=:table, head=["Moment","Data","\\multicolumn{1}{p{2.5cm}}{\\centering Altruism \\\\ \$ \\eta=$(round(Mb.mp.η;digits=3)) \$}","\\multicolumn{1}{p{2.5cm}}{\\centering No Altruism \\\\ \$ \\eta=0 \$}"],latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true)
    df_latex = replace(df_latex,"cccc" => "lrrr")
    df_latex = replace(df_latex,"missing" => "")
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/postest/quant_pres_"*append_info*".tex"
    filename2 = "tabfig/postest/quant_mostrecent.tex"
    println(df_latex)
    write(filename, df_latex)
    if store == true
        write(filename2, df_latex)
    end
    println("Output stored as ",filename)

end

function output_decomp(df_decomp,targ_moms::Vector{Symbol},nontarg_moms::Vector{Symbol};store=false)
    df = df_decomp[df_decomp.type.== "Model",:]
    out_moms =  union(targ_moms,nontarg_moms)
    
    keepvars = [:friction,:altr]
    append!(keepvars,out_moms)

    df = df[!,keepvars]
    out = DataFrame(varname = names(df[:,out_moms]))

    for istrue in [true, false]
        for frictions in unique(df.friction)
            display([istrue,frictions])
            out[!,Symbol(frictions*"_"*string(istrue))] = collect(vec(Matrix(df[(df.friction .== frictions) .& (df.altr .== istrue),out_moms])))
        end
    end

    _varnames = varnames()
    out[!,:varname] = ["\\;"*_varnames[Symbol(var)] for var in out[!,:varname]]

    insert!.(eachcol(out), length(targ_moms)+1, ["\\textit{Non-Targeted Moments}",missing,missing,missing,missing,missing,missing,missing,missing])
    insert!.(eachcol(out), 1, ["\\textit{Targeted Moments}",missing,missing,missing,missing,missing,missing,missing,missing])

    out_narrow = out[!,Symbol.(["varname", "Benchmark_true" , "LTV_true", "Adjustment Costs_true","Income Risk_true",
            "Benchmark_false" , "LTV_false", "Adjustment Costs_false","Income Risk_false"])]

    display(out_narrow)
    df_latex = latexify(out_narrow; env=:table, head=["Moment",
        "\\multicolumn{1}{c}{Bench}", "\\multicolumn{1}{c}{LTV 0.9}","\\multicolumn{1}{c}{No Cost}","\\multicolumn{1}{c}{Inc Cer.}",
        "\\multicolumn{1}{c}{Bench}", "\\multicolumn{1}{c}{LTV 0.9}","\\multicolumn{1}{c}{No Cost}","\\multicolumn{1}{c}{Inc Cer.}"],latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true)
    df_latex = replace(df_latex,"ccccccccc" => "l rrrr rrrr")
    df_latex = replace(df_latex,"missing" => "")
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/postest/decomp.tex"
    filename2 = filename_wdate(filename)
    println("Output stored as ",filename)
    write(filename, df_latex)
    if store == true
        write(filename2, df_latex)
    end

    out_narrow = out[!,Symbol.(["varname", "Benchmark_true" , "Benchmark_false" ,
            "LTV_true", "LTV_false",
             "Adjustment Costs_true", "Adjustment Costs_false",
             "Income Risk_true",           "Income Risk_false"])]
     display(out_narrow)
     df_latex = latexify(out_narrow; env=:table, head=["Moment",
         "\\multicolumn{1}{c}{Altr}", "\\multicolumn{1}{c}{No Altr}","\\multicolumn{1}{c}{Altr}", "\\multicolumn{1}{c}{No Altr}",
         "\\multicolumn{1}{c}{Altr}", "\\multicolumn{1}{c}{No Altr}","\\multicolumn{1}{c}{Altr}", "\\multicolumn{1}{c}{No Altr}"],latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true)
     df_latex = replace(df_latex,"ccccccccc" => "l rrrr rrrr")
     df_latex = replace(df_latex,"\\toprule\n" => "\\toprule\n & \\multicolumn{2}{c}{Benchmark} & \\multicolumn{2}{c}{No LTV} & \\multicolumn{2}{c}{Liquid Housing} & \\multicolumn{2}{c}{No Income Risk}   \\\\  \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9} \n" )
     df_latex = replace(df_latex,"missing" => "")

     append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
     filename = "tabfig/postest/decomp_sort.tex"
     filename2 = filename_wdate(filename)
     write(filename, df_latex)
     if store == true
         write(filename2, df_latex)
     end

     println("The following lines will tell you how much of the decrease is driven by what")
     i = 2
     bench = (out_narrow[4,i+1] - out_narrow[4,i])/out_narrow[4,i]
     decr = (out_narrow[4,i+1] - out_narrow[4,i])
     tmp = decr/out_narrow[4,i]
     explain = 1-tmp/bench
     println(round.([bench,decr,tmp,explain];digits=4))
     i = 4
     decr = (out_narrow[4,i+1] - out_narrow[4,i])
     tmp = decr/out_narrow[4,i]
     explain = 1-tmp/bench
     println(round.([bench,decr,tmp,explain];digits=4))
     i = 6
     decr = (out_narrow[4,i+1] - out_narrow[4,i])
     tmp = decr/out_narrow[4,i]
     explain = 1-tmp/bench
     println(round.([bench,decr,tmp,explain];digits=4))
     i = 8
     decr = (out_narrow[4,i+1] - out_narrow[4,i])
     tmp = decr/out_narrow[4,i]
     explain = 1-tmp/bench
     println(round.([bench,decr,tmp,explain];digits=4))


end
function output_policy(df_decomp,targ_moms::Vector{Symbol},bypwealth_moms::Vector{Symbol};store=false)
    df = df_decomp[df_decomp.type.== "Model",:]
    out_moms =  union(targ_moms,bypwealth_moms)
    @assert minimum(df.transfrateyoung) > 0
    keepvars = [:friction,:oldval,:newval]
    append!(keepvars,out_moms)

    df = df[!,keepvars]
    dfout = DataFrame(varname = names(df[:,out_moms]))

    for frictions in unique(df.friction)
        dfout[!,Symbol(frictions)] = collect(vec(Matrix(df[(df.friction .== frictions),out_moms])))
    end


    dfout[!,:varname] = ["\\;"*varnames()[Symbol(var)] for var in dfout[!,:varname]]
    
    insert!.(eachcol(dfout), length(targ_moms)+1, ["\\textit{By Initial Parent Wealth}",missing,missing,missing,missing,missing,missing])
    insert!.(eachcol(dfout), 1, ["\\textit{Targeted Moments}",missing,missing,missing,missing,missing,missing,])

    
    out_narrow = dfout[!,Symbol.(["varname", "Benchmark" , "LTV", "PMI", "Price", "Purchase Costs","Sales Cost"])]

    display(out_narrow)
    df_latex = latexify(out_narrow; env=:table, head=["Moment",
        "\\multicolumn{1}{c}{Bench.}", "\\multicolumn{1}{c}{LTV}","\\multicolumn{1}{c}{\$PMI\$}","\\multicolumn{1}{c}{\$Price\$}","\\multicolumn{1}{c}{\$m_b\$ }","\\multicolumn{1}{c}{\$m_s\$}"]
        ,latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true)

    df_latex = replace(df_latex,"ccccccc" => "l rrrrrr")
    df_latex = replace(df_latex,"missing" => "")

    ## Next, add top rows of old and new parameter values
    oldvals = "\\; Old Parameter value"
    for val in df.oldval
        oldvals *= " & $(round(val; digits=3))"
    end
    oldvals *= " \\\\ \n"
    newvals = ("\\; New Parameter value")
    for val in df.newval
        newvals *= (" & $(string(round(val;digits=3)))" ) 
    end
    newvals *= (" \\\\ \n")


    df_latex = replace(df_latex,"\\midrule" => oldvals * newvals * "\\midrule")
    df_latex = replace(df_latex,"NaN" => "")


    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/postest/policy_pres_"*append_info*".tex"
    filename2 = "tabfig/postest/policy_mostrecent.tex"
    write(filename, df_latex)
    if store == true
        write(filename2, df_latex)
    end
    println("Output stored as ",filename)

end

##

function policyeffects(Mbench::Model)
    _mp = deepcopy(Mbench.mp)
    nstate=Mbench.np.nxk
    nchoice=Mbench.np.nck_choice
    df = Mbench.moms
    df = leftjoin(df,DataFrame(type=["Data", "Model"], friction =[missing,"Benchmark"], oldval = [missing,NaN], newval =[missing,NaN]),on = :type)

    ## Reduce price
    mp, np = benchpar(nstate,nchoice)
    cval =  _mp.price*0.98 - _mp.price
    update_prices!(mp.price+cval,np)
    MbPrice = Model(mp,np)
    solvemodel2nd!(MbPrice,np,mp)
    # Append moments
    moms = deepcopy(MbPrice.moms)
    moms[!,:friction] .= "Price"
    moms[!,:newval] .= _mp.price + cval 
    moms[!,:oldval] .= _mp.price 
    append!(df,moms[moms.type.=="Model",:])

    ## Reduce LTV
    mp, np = benchpar(nstate,nchoice)
    cval = -0.03
    mp.LTV .+= cval
    MbLTV = Model(mp,np)
    solvemodel2nd!(MbLTV,np,mp)
    # Append moments
    moms = deepcopy(MbLTV.moms)
    moms[!,:friction] .= "LTV"
    moms[!,:newval] .= _mp.LTV[1] + cval
    moms[!,:oldval] .= _mp.LTV[1] 
    append!(df,moms[moms.type.=="Model",:])

    ## Reduce PMI 
    mp, np = benchpar(nstate,nchoice)
    cval = +0.07
    mp.PMIlim += cval
    MbPMI = Model(mp,np)
    solvemodel2nd!(MbPMI,np,mp)
    # Append moments
    moms = deepcopy(MbPMI.moms)
    moms[!,:friction] .= "PMI"
    moms[!,:newval] .= _mp.PMIlim + cval
    moms[!,:oldval] .= _mp.PMIlim 
    append!(df,moms[moms.type.=="Model",:])

    ## Remove Purhcase costs
    mp, np = benchpar(nstate,nchoice)
    cval = -0.02
    mp.mb .+= cval
    Mbadj = Model(mp,np)
    solvemodel2nd!(Mbadj,np,mp)
    # Append moments
    moms = deepcopy(Mbadj.moms)
    moms[!,:friction] .= "Purchase Costs"
    moms[!,:newval] .= _mp.mb[1] + cval
    moms[!,:oldval] .= _mp.mb[1] 
    append!(df,moms[moms.type.=="Model",:])

    ## Remove Sales costs
    mp, np = benchpar(nstate,nchoice)
    cval = -0.02
    mp.ms .+= cval
    Mbsal = Model(mp,np)
    solvemodel2nd!(Mbsal,np,mp)
    # Append moments
    moms = deepcopy(Mbsal.moms)
    moms[!,:friction] .= "Sales Cost"
    moms[!,:newval] .= _mp.ms[1] + cval
    moms[!,:oldval] .= _mp.ms[1] 
    append!(df,moms[moms.type.=="Model",:])

    sort!(df,[:friction,:type])
    return df
end
