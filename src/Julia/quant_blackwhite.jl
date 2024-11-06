
function blackwhitegap(Mbench::Model)
    gaps = DataFrame(CSV.File("data/model/calibration/racialincomegaps.csv"))
    racdiff = ( # element of all racial differenctes
        bincome = gaps[gaps.race.=="Black",:gap][1],
        wincome = gaps[gaps.race.=="White",:gap][1],
        bhsize = 0.75,
        whsize = 1.0
    )    
    display("Racial income gap: White=$(round(racdiff.wincome;digits=3)). Black=$(round(racdiff.bincome;digits=3)). ")
    
    ## initialize
    moms = deepcopy(Mbench.moms)
    moms[!,:race] = ["bench", "bench"]
    
    
    # Black
    mp = deepcopy(Mbench.mp)
    np = deepcopy(Mbench.np)
    Mblack = Model(mp,np)
    np.h_grd .*=  racdiff.bhsize
    np.y_grd .*=  racdiff.bincome
    display(np.h_grd)
    display(np.y_grd)
    solvemodel2nd!(Mblack,np,mp)
    Mblack.moms[!, :race] .= ["Black w altr"]

    Mblack_n  = altrcontrhousing(mp,np)
    Mblack_n.moms[!, :race] .= ["Black no altr"]

    # White    
    mp = deepcopy(Mbench.mp)
    np = deepcopy(Mbench.np)
    Mwhite = Model(mp,np)
    np.h_grd .*= racdiff.whsize
    np.y_grd .*= racdiff.wincome
    display(np.h_grd)
    display(np.y_grd)
    solvemodel2nd!(Mwhite,np,mp)
    Mwhite.moms[!, :race] .= ["White w altr"]

    Mwhite_n  = altrcontrhousing(mp,np)
    Mwhite_n.moms[!, :race] .= ["White no altr"]

    Racemoms = vcat(moms,Mwhite.moms,Mwhite_n.moms,Mblack.moms,Mblack_n.moms)
    Racemoms = Racemoms[(Racemoms.type .== "Model") .| (Racemoms.race .== "bench") ,:]

    return Racemoms


end

function output_BW(BW::DataFrame,targ_moms::Vector{Symbol},nontarg_moms::Vector;store = false)
    BW_data = DataFrame(CSV.File("data/moments/moments_race_age.csv"))
    out_moms =  union(targ_moms,nontarg_moms)
    BW_out = BW[BW.race.!= "bench",:]
    BW_data = BW_data[BW_data.race .!= "All",:]
    BW_data.race = replace(BW_data.race, "Black" => "Black data")
    BW_data.race = replace(BW_data.race, "White" => "White data")



    BW_out = vcat(BW_out,BW_data,cols=:union)

    sort!(BW_out,:race,rev=true)
    filter!(x->x≠:race,out_moms)
    # push!(out_moms,:race)
    # Transpose the data

    out = DataFrame(varname = names(BW_out[:,out_moms]),
                    WhiteData = vec(Matrix(BW_out[BW_out.race .== "White data",out_moms])),
                    WhiteAltr = vec(Matrix(BW_out[BW_out.race .== "White w altr",out_moms])),
                    WhiteNoAltr = vec(Matrix(BW_out[BW_out.race .== "White no altr",out_moms])),
                    BlackData = vec(Matrix(BW_out[BW_out.race .== "Black data",out_moms])),
                    BlackAltr = vec(Matrix(BW_out[BW_out.race .== "Black w altr",out_moms])),
                    BlackNoAltr = vec(Matrix(BW_out[BW_out.race .== "Black no altr",out_moms])),
                    )
    allowmissing!(out,Not(:varname)) # Allowmissing in all columns

    _varnames = varnames()
    out[!,:varname] = ["\\;"*_varnames[Symbol(var)] for var in out[!,:varname]]
    if isempty(nontarg_moms) == false # Only insert these if you are keeping other moments
        insert!.(eachcol(out), 1, ["\\textit{Targeted Moments}"; fill(missing,nrow(BW_out))])
        insert!.(eachcol(out), length(targ_moms)+2, ["\\textit{Non-Targeted Moments}"; fill(missing,nrow(BW_out))])
    end

    display(out)
    BW_latex = latexify(out; env=:table, latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,
        head= ["\\multicolumn{1}{c}{Moment}", "\\multicolumn{1}{c}{Data}", "\\multicolumn{1}{c}{Altr.}", "\\multicolumn{1}{c}{No Altr.}", "\\multicolumn{1}{c}{Data}", "\\multicolumn{1}{c}{Altr.}", "\\multicolumn{1}{c}{No Altr.}"])
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    BW_latex = replace(BW_latex,"\\toprule" => "\\toprule & \\multicolumn{3}{c}{White} & \\multicolumn{3}{c}{Black} \\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}")
    BW_latex = replace(BW_latex,"ccccccc" => "l rrr rrr")
    BW_latex = replace(BW_latex,"missing" => "")
    println(BW_latex)
    filename = "tabfig/postest/BW_"*append_info*".tex"
    filename2 = "tabfig/postest/BW_mostrecent.tex"

    write(filename, BW_latex)
    println("Output stored as ",filename)
    if store == true
        write(filename2, BW_latex)
    end

end
