function robustness(nstate,nchoice)
    # Parental wealth shocks
    mpν, npν = benchpar(nstate,nchoice;robust_ν=true,robust_s=false)
    Mbν = Model(mpν,npν)
    solvemodel2nd!(Mbν,npν,mpν)
    Mnν = altrcontrhousing(Mbν.mp,Mbν.np)
 
    # Income shocks
    mps, nps = benchpar(nstate,nchoice;robust_ν=false,robust_s=true)
    Mbs = Model(mps,nps)
    solvemodel2nd!(Mbs,nps,mps)
    Mns = altrcontrhousing(Mbs.mp,Mbs.np)

    return Mbν,Mnν,Mbs,Mns
end



function output_robust(Mb::Model,Mn::Model,Mbν::Model,Mnν::Model,Mbs::Model,Mns::Model,
        targ_moms::Vector{Symbol},nontarg_moms::Vector{Symbol};store=false)
    @assert Mb.moms[Mb.moms.type.=="Model",:transfrateyoung][1] > 0.0
    @assert Mn.moms[Mn.moms.type.=="Model",:transfrateyoung][1] == 0.0
    @assert Mnν.moms[Mnν.moms.type.=="Model",:transfrateyoung][1] == 0.0
    @assert Mns.moms[Mns.moms.type.=="Model",:transfrateyoung][1] == 0.0

    out_moms =  union(targ_moms,nontarg_moms)
    # out_moms = [:owneryoung, :medwealthyoung, :medwealthold, :wealthatpurchase, :transfrate, :r2incyoung]
    df = DataFrame(title=["Data", "Altruism", "Without"])
    temp = vcat(Mb.moms,Mn.moms[Mn.moms.type .== "Model",:])[:,out_moms]
    df = hcat(df,temp)

    ##

    out = DataFrame(varname = names(df[:,2:end]),
        Data = vec(Matrix(df[df.title .=="Data",out_moms])),
        Altruism = vec(Matrix(df[df.title .=="Altruism",out_moms])),
        Without = vec(Matrix(df[df.title .=="Without",out_moms])),
        Altruism_oldrisk = vec(Matrix(Mbν.moms[Mbν.moms.type .=="Model",out_moms])),
        Without_oldrisk = vec(Matrix(Mnν.moms[Mnν.moms.type .=="Model",out_moms])),
        Altruism_pricerisk = vec(Matrix(Mbs.moms[Mbs.moms.type .=="Model",out_moms])),
        Without_pricerisk = vec(Matrix(Mns.moms[Mns.moms.type .=="Model",out_moms])),
        )
    allowmissing!(out,Not(:varname)) # Allowmissing in all columns

    out[!,:varname] = ["\\;"*varnames()[Symbol(var)] for var in out[!,:varname]]
    insert!.(eachcol(out), length(targ_moms)+1,  ["\\textit{Non-Targeted Moments}"; fill(missing,ncol(out)-1)])
    insert!.(eachcol(out), 1, ["\\textit{Targeted Moments}"; fill(missing,ncol(out)-1)])
    display(out)
    df_latex = latexify(out; env=:table,latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true,
        head=["Moment","Data","\\multicolumn{1}{c}{\$ \\eta>0 \$}","\\multicolumn{1}{c}{\$ \\eta=0 \$}","\\multicolumn{1}{c}{\$ \\eta>0 \$}","\\multicolumn{1}{c}{\$ \\eta=0 \$}","\\multicolumn{1}{c}{\$ \\eta>0 \$}","\\multicolumn{1}{c}{\$ \\eta=0 \$}",])
    df_latex = replace(df_latex,"cccccccc" => "lr rr rr rr")
    df_latex = replace(df_latex,"\n\\toprule" => "\n \\toprule & & \\multicolumn{2}{c}{Benchmark} & \\multicolumn{2}{c}{Old Risk} & \\multicolumn{2}{c}{Price Risk}\n \\\\  \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}")
    df_latex = replace(df_latex,"missing" => "")
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/postest/robust_"*append_info*".tex"
    filename2 = "tabfig/postest/robust_mostrecent.tex"
    write(filename, df_latex)
    println(df_latex)
    if store == true
        write(filename2, df_latex)
    end
    println("Output stored as ",filename)

end
