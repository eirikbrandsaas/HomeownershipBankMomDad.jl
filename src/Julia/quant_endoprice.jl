function SolveModelDifferentPrices(Min,pgrd)
    mp = Min.mp
    np = Min.np
    df = DataFrame()
    for price in pgrd
        update_prices!(price,np)
        _M = altrcontrhousing(mp,np)
        demand = _M.moms[_M.moms.type.=="Model",:].ownerall[1]
        append!(df, [_M.moms[_M.moms.type .== "Model",:] DataFrame(price=price,demand = demand)])
        @assert np.p_grd[1]==price # Check that the correct prices are used
    end

    return df
end



function output_endoprice(dfp,Mb,moms,agg_moms,α;store=false)
    keep = union(agg_moms,moms)
    Bmoms= Mb.moms[Mb.moms.type.=="Model",:]

    # Add variables to model
    Bmoms[!,:price] .= Mb.mp.price
    Bmoms[!,:α1] .= NaN
    Bmoms[!,:supply] .= Mb.moms[Mb.moms.type.=="Model",:ownerall]
    df= copy(Bmoms)
    df = vcat(df,dfp,cols=:union)
    # df[!,:α1] = [NaN;α]

    df = df[!,keep]
    df = select(df,keep)
    @assert issorted(α) "Supply elasticities not sorted. Next block assumes that they are sorted"
    out = DataFrame(Moment = names(df),
        Benchmark = Matrix(df[isnan.(df.α1),:])[1,:],
        Elastic = Matrix(df[df.α1.==α[3],:])[1,:],
        Middle = Matrix(df[df.α1.==α[2],:])[1,:],
        Low = Matrix(df[df.α1.==α[1],:])[1,:],
        )
    out[!,:Moment] = ["\\;"*varnames()[Symbol(var)] for var in out[!,:Moment]]
    insert!.(eachcol(out), 1, ["\\textit{Aggregate Moments}"; fill(NaN,length(α)+1)])
    insert!.(eachcol(out), length(aggmoms)+2, ["\\textit{Targeted Moments}"; fill(NaN,length(α)+1)])
    display(out)
    df_latex = latexify(out; env=:table,latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true)

    df_latex = replace(df_latex,"ccccc}\n\\toprule" => "lrrrr}\n \\toprule & \\multicolumn{1}{c}{Altruism} & \\multicolumn{$(length(α))}{c}{Without Altruism}  \\\\  \\cmidrule(lr){2-2} \\cmidrule(lr){3-$(length(α)+2)} \n")
    df_latex = replace(df_latex,"NaN" => "")
    df_latex = replace(df_latex,"& Inf &" => "& \$ \\infty \$ &")
    # df_latex = replace(df_latex,"\\begin{tabular}{lrrrr}\n" =>"test")
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/postest/quant_endogenprice_"*append_info*".tex"
    filename2 = "tabfig/postest/quant_endogenprice_mostrecent.tex"

    write(filename, df_latex)
    if store == true
        write(filename2, df_latex)
    end
end
##



function housedemand(Min)
    M = altrcontrhousing(Min.mp,Min.np)
    # agg = aggregate(M.simpan,[sum, mean])
    # HD = (agg.ownk_mean + agg.ownp_mean)[1]/2
    HD = M.moms[M.moms.type.=="Model",:].ownerall[1]

    return HD, M
end

function supply_func(α0,α1,price)
    exp(α0 + α1*log(price))
end

function Hsupply_intercept(Hbar,Pbar,α1)
    α0 = log(Hbar) - α1*log(Pbar)
end
