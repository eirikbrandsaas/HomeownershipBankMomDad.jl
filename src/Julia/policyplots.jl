function plot_policyfuncs(Mb,Mn;store=false)
    mp = Mb.mp
    np = Mb.np

    ## Setup for which state to plot
    bph = 400
    ibph = argmin(abs.(np.bp_grd.-bph))

    bpl = 100
    ibpl = argmin(abs.(np.bp_grd.-bpl))

    ihp = 1
    ivk = 2
    iok = 1
    is = 1
    iak = 1
    ihk = 2
    _h2size = Dict(1=>L"h_r",2=>L"h_o")

    plot(xlabel="Child Cash-on-Hand $(L"\tilde x_c")",ylabel="Child Bond Choice $(L"b'_c")")

    iown = findfirst(x->x==true,Mn.g.gk.disc[:,ibpl,ihp,ivk,iok,is,iak,ihk].==1)
    plot!(np.xk_grd[1:iown],Mn.g.gk.b[1:iown,ibpl,ihp,ivk,iok,is,iak,1] ,label="No Altruism",marker=:circle,color=:black,linewidth=2)
    plot!(np.xk_grd[iown:end],Mn.g.gk.b[iown:end,ibpl,ihp,ivk,iok,is,iak,2] ,label="",marker=:circle,color=:black,linewidth=2)
    vline!([np.xk_grd[iown]],line=(:black,:dash),label="",linewidth=2)

    iown = findfirst(x->x==true,Mb.g.gk.disc[:,ibpl,ihp,ivk,iok,is,iak,ihk].==1)
    plot!(np.xk_grd[1:iown],Mb.g.gk.b[1:iown,ibpl,ihp,ivk,iok,is,iak,1] ,label="Altruism, Non-Rich Parent",marker=:xcross,color=:brown,linewidth=2)
    plot!(np.xk_grd[iown:end],Mb.g.gk.b[iown:end,ibpl,ihp,ivk,iok,is,iak,2] ,label="",marker=:xcross,color=:brown,linewidth=2)
    vline!([np.xk_grd[iown]],line=(:brown,:dash),label="",linewidth=2)

    iown = findfirst(x->x==true,Mb.g.gk.disc[:,ibph,ihp,ivk,iok,is,iak,ihk].==1)
    plot!(np.xk_grd[1:iown],Mb.g.gk.b[1:iown,ibph,ihp,ivk,iok,is,iak,1] ,label="Altruism, Rich Parent",marker=:cross,color=:orange,linewidth=2)
    plot!(np.xk_grd[iown:end],Mb.g.gk.b[iown:end,ibph,ihp,ivk,iok,is,iak,2] ,label="",marker=:cross,color=:orange,linewidth=2)
    vline!([np.xk_grd[iown]],line=(:orange,:dash),label="",linewidth=2)

    xlims!(minimum(np.xk_grd)-5,350)
    p2 = plot!(ylims=(-220,200),size=(600,350).*0.8)
   
    if store == true
        Plots.pdf(p2,"tabfig/model/pol_bond")
        open("tabfig/model/pol_bond_state.tex", "w") do file
            write(file, "$(L"h'_p")=$(_h2size[ihp]), $(L"y_c")=$(round(np.v_grd[iak,ivk];digits=1)), $(L"h_c")=$(_h2size[ihk]), and $(L"a_c")=$(round(Int,np.age_grd[iak])). The non-rich and rich parent have $(L"b'_p")=$(round(Int,bpl)) and $(L"b'_p")=$(round(Int,bph))" )
        end

    
    end
    
    display(p2)
end

